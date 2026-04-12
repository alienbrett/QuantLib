/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/pricingengines/vanilla/vnbinomialengine.hpp>
#include <ql/exercise.hpp>
#include <ql/math/distributions/binomialdistribution.hpp>
#include <ql/pricingengines/greeks.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    VNBinomialVanillaEngine::VNBinomialVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size timeSteps)
    : process_(std::move(process)),
      dividends_(std::move(dividends)),
      timeSteps_(timeSteps) {
        QL_REQUIRE(timeSteps >= 2,
                   "at least 2 time steps required, "
                   << timeSteps << " provided");
        registerWith(process_);
    }

    void VNBinomialVanillaEngine::calculate() const {

        DayCounter rfdc  = process_->riskFreeRate()->dayCounter();
        DayCounter divdc = process_->dividendYield()->dayCounter();

        Real S0 = process_->stateVariable()->value();
        QL_REQUIRE(S0 > 0.0, "negative or null underlying given");

        Date referenceDate = process_->riskFreeRate()->referenceDate();
        Date maturityDate  = arguments_.exercise->lastDate();
        Time T = rfdc.yearFraction(referenceDate, maturityDate);

        Rate r = process_->riskFreeRate()->zeroRate(
            maturityDate, rfdc, Continuous, NoFrequency);
        Rate q = process_->dividendYield()->zeroRate(
            maturityDate, divdc, Continuous, NoFrequency);
        Volatility sigma = process_->blackVolatility()->blackVol(
            maturityDate, S0);

        auto payoff =
            ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");
        Real K = payoff->strike();
        QL_REQUIRE(K > 0.0, "strike must be positive");

        bool isAmerican =
            (arguments_.exercise->type() == Exercise::American);

        // Earliest exercise time (for American options with a window)
        Time earliestExercise = 0.0;
        if (arguments_.exercise->dates().size() > 1)
            earliestExercise = rfdc.yearFraction(
                referenceDate, arguments_.exercise->date(0));

        // -----------------------------------------------------------------
        // Collect fixed cash dividends within (referenceDate, maturityDate].
        // Each entry is (time, dollar amount).  The BSM process carries any
        // continuous dividend yield; discrete dividends here are additive
        // cash amounts that reduce the stock price at the ex-date.
        // -----------------------------------------------------------------
        struct DivInfo { Time t; Real amount; };
        std::vector<DivInfo> cashDivs;

        for (const auto& div : dividends_) {
            if (div->date() > referenceDate && div->date() <= maturityDate) {
                Time td = rfdc.yearFraction(referenceDate, div->date());
                cashDivs.push_back({td, div->amount()});
            }
        }
        std::sort(cashDivs.begin(), cashDivs.end(),
                  [](const DivInfo& a, const DivInfo& b) {
                      return a.t < b.t;
                  });

        // -----------------------------------------------------------------
        // Leisen-Reimer tree geometry (odd number of steps required).
        // up/down/d2 use terminal flat rates — they set the recombining
        // grid and LR strike-centering.  Per-step probabilities and
        // discounts are derived below from the actual term structures
        // to eliminate flat-rate forward bias at intermediate nodes.
        // -----------------------------------------------------------------
        Size N = (timeSteps_ % 2 == 0) ? timeSteps_ + 1 : timeSteps_;
        Time dt = T / N;  // average dt — for LR geometry only

        Real variance     = sigma * sigma * T;
        Real sqrtVar      = std::sqrt(variance);
        // Terminal drift — for LR geometry only
        Real driftPerStep = (r - q - 0.5 * sigma * sigma) * dt;
        Real ermqdt       = std::exp(driftPerStep + 0.5 * variance / N);

        Real d2 = (std::log(S0 / K) + driftPerStep * N) / sqrtVar;
        Real pu_lr = PeizerPrattMethod2Inversion(d2, N);
        Real pdash = PeizerPrattMethod2Inversion(d2 + sqrtVar, N);
        Real up   = ermqdt * pdash / pu_lr;
        Real down = (ermqdt - pu_lr * up) / (1.0 - pu_lr);

        QL_ENSURE(up > down,
                  "LR tree: up (" << up << ") must exceed down (" << down << ")");

        // -----------------------------------------------------------------
        // Non-uniform time grid: insert exact ex-dividend dates as
        // mandatory nodes, then fill in uniform sub-steps within each
        // segment.  This eliminates snapping error that causes erratic
        // convergence with multiple dividends.  Total step count stays
        // at N; only the per-step dt varies.
        // -----------------------------------------------------------------
        std::vector<Time> stepTimes(N + 1);
        std::vector<Size> divStep;
        std::vector<Real> divAmt;

        if (cashDivs.empty()) {
            // Uniform grid — no dividends to align
            for (Size k = 0; k <= N; ++k)
                stepTimes[k] = k * dt;
        } else {
            // Build segment breakpoints: [0, t_d1, t_d2, ..., T]
            std::vector<Time> breaks;
            breaks.push_back(0.0);
            for (const auto& d : cashDivs)
                if (d.t > 1e-10 && d.t < T - 1e-10)
                    breaks.push_back(d.t);
            breaks.push_back(T);

            Size nSeg = breaks.size() - 1;

            // Allocate steps to each segment proportional to duration,
            // minimum 2 per segment to avoid degenerate single-step
            // segments that can't support interpolation.
            std::vector<Size> segSteps(nSeg);
            Size minPerSeg = 2;
            Size reserved = minPerSeg * nSeg;
            QL_REQUIRE(N >= reserved,
                       "VN tree: need at least " << reserved
                       << " steps for " << nSeg << " segments, got " << N);
            Size pool = N - reserved;
            Size allocated = 0;
            for (Size s = 0; s < nSeg; ++s) {
                Real frac = (breaks[s+1] - breaks[s]) / T;
                Size alloc = static_cast<Size>(std::round(pool * frac));
                segSteps[s] = minPerSeg + alloc;
                allocated += alloc;
            }
            // Reconcile rounding: give/take from longest segment
            Size longest = 0;
            for (Size s = 1; s < nSeg; ++s)
                if (breaks[s+1] - breaks[s] > breaks[longest+1] - breaks[longest])
                    longest = s;
            if (allocated < pool)
                segSteps[longest] += (pool - allocated);
            else if (allocated > pool && segSteps[longest] > minPerSeg + (allocated - pool))
                segSteps[longest] -= (allocated - pool);

            // Fill in step times: uniform within each segment
            Size idx = 0;
            stepTimes[0] = 0.0;
            for (Size s = 0; s < nSeg; ++s) {
                Time segStart = breaks[s];
                Time segLen = breaks[s+1] - breaks[s];
                Time segDt = segLen / segSteps[s];
                for (Size k = 1; k <= segSteps[s]; ++k)
                    stepTimes[++idx] = segStart + k * segDt;
            }

            // Map dividends to their exact step indices (segment boundaries)
            for (Size di = 0; di < cashDivs.size(); ++di) {
                Time td = cashDivs[di].t;
                // Find the step index at this dividend time
                auto it = std::lower_bound(
                    stepTimes.begin(), stepTimes.end(), td - 1e-10);
                Size s = static_cast<Size>(
                    std::distance(stepTimes.begin(), it));
                // Clamp to [1, N-1]
                s = std::max<Size>(1, std::min<Size>(s, N - 1));
                if (!divStep.empty() && divStep.back() == s)
                    divAmt.back() += cashDivs[di].amount;
                else {
                    divStep.push_back(s);
                    divAmt.push_back(cashDivs[di].amount);
                }
            }
        }
        Size nDiv = divStep.size();

        // -----------------------------------------------------------------
        // Per-step forward rates from actual term structures.
        // Precompute discount factors at each (non-uniform) step time,
        // then derive local growth/discount per step.  This matches
        // intermediate-node forwards to the true term structure.
        // -----------------------------------------------------------------
        std::vector<Real> dfRf(N + 1), dfQ(N + 1);
        dfRf[0] = dfQ[0] = 1.0;
        for (Size k = 1; k <= N; ++k) {
            dfRf[k] = process_->riskFreeRate()->discount(stepTimes[k]);
            dfQ[k]  = process_->dividendYield()->discount(stepTimes[k]);
        }

        std::vector<Real> discStep(N), puStep(N), pdStep(N);
        for (Size k = 0; k < N; ++k) {
            // One-step risk-free discount: df(t_{k+1}) / df(t_k)
            discStep[k] = dfRf[k + 1] / dfRf[k];
            // One-step forward growth of stock: F(t_{k+1}) / F(t_k)
            Real growth = (dfQ[k + 1] * dfRf[k]) / (dfQ[k] * dfRf[k + 1]);
            // Local risk-neutral probability matching the true forward.
            // Clamp to (0,1) — can be exceeded when the carry curve has
            // a proportional dividend jump that creates extreme growth in
            // a single step (e.g. short DTE near ex-date).  Clamping is
            // the tree's best approximation; error is localized to one
            // step and does not introduce arbitrage.
            puStep[k] = std::max(0.001, std::min(0.999,
                            (growth - down) / (up - down)));
            pdStep[k] = 1.0 - puStep[k];
        }

        // -----------------------------------------------------------------
        // Helper: stock price at node (step i, index j)
        //   S_{i,j} = S0 * up^j * down^(i-j)
        // Precompute power tables to avoid repeated std::pow in hot loops.
        // -----------------------------------------------------------------
        std::vector<Real> upPow(N + 1), dnPow(N + 1);
        upPow[0] = dnPow[0] = 1.0;
        for (Size k = 1; k <= N; ++k) {
            upPow[k] = upPow[k-1] * up;
            dnPow[k] = dnPow[k-1] * down;
        }

        auto stockPrice = [&](Size i, Size j) -> Real {
            return S0 * upPow[j] * dnPow[i - j];
        };

        // -----------------------------------------------------------------
        // Terminal payoff
        // -----------------------------------------------------------------
        Array V(N + 1);
        for (Size j = 0; j <= N; ++j)
            V[j] = (*payoff)(stockPrice(N, j));

        // -----------------------------------------------------------------
        // Backward induction with VN interpolation at ex-div dates
        // -----------------------------------------------------------------
        // Dividend pointer walks backwards through the sorted list.
        Size divIdx = nDiv;

        // Variables to capture values at steps 2, 1 for Greeks
        Real p2d = 0, p2m = 0, p2u = 0;
        Real s2d = 0, s2m = 0, s2u = 0;
        Real p1d = 0, p1u = 0;
        Real s1d = 0, s1u = 0;

        for (Integer i = N - 1; i >= 0; --i) {
            Size si = static_cast<Size>(i);

            // --- standard one-step backward induction ---
            // Uses per-step discount and probabilities from actual term
            // structures (not flat rates) to match local forwards.
            Array newV(si + 1);
            Real di = discStep[si], pui = puStep[si], pdi = pdStep[si];
            for (Size j = 0; j <= si; ++j)
                newV[j] = di * (pdi * V[j] + pui * V[j + 1]);

            // --- VN interpolation if this step is an ex-div date ---
            if (divIdx > 0 && divStep[divIdx - 1] == si) {
                --divIdx;
                Real D = divAmt[divIdx];

                // Stock prices at this step (monotonically increasing in j)
                std::vector<Real> prices(si + 1);
                for (Size j = 0; j <= si; ++j)
                    prices[j] = stockPrice(si, j);

                Array adjV(si + 1);
                for (Size j = 0; j <= si; ++j) {
                    Real Spost = prices[j] - D;

                    if (Spost <= 0.0) {
                        // dividend exceeds stock price
                        adjV[j] = (*payoff)(std::max(Spost, 0.0));
                    } else {
                        // linear interpolation / extrapolation
                        auto it = std::lower_bound(
                            prices.begin(), prices.end(), Spost);

                        Size j1, j0;
                        if (it == prices.begin()) {
                            j0 = 0; j1 = 1;      // extrapolate below
                        } else if (it == prices.end()) {
                            j1 = si; j0 = si - 1; // extrapolate above
                        } else {
                            j1 = static_cast<Size>(
                                std::distance(prices.begin(), it));
                            j0 = j1 - 1;
                        }
                        Real w = (Spost - prices[j0])
                               / (prices[j1] - prices[j0]);
                        adjV[j] = newV[j0] + w * (newV[j1] - newV[j0]);
                    }
                }
                newV = adjV;
            }

            // --- early exercise (American) ---
            if (isAmerican) {
                Time stepTime = stepTimes[si];
                if (stepTime >= earliestExercise) {
                    for (Size j = 0; j <= si; ++j)
                        newV[j] = std::max(newV[j],
                                           (*payoff)(stockPrice(si, j)));
                }
            }

            // --- capture nodes for Greeks ---
            if (si == 2) {
                p2d = newV[0]; p2m = newV[1]; p2u = newV[2];
                s2d = stockPrice(2, 0);
                s2m = stockPrice(2, 1);
                s2u = stockPrice(2, 2);
            }
            if (si == 1) {
                p1d = newV[0]; p1u = newV[1];
                s1d = stockPrice(1, 0);
                s1u = stockPrice(1, 1);
            }

            V = newV;
        }

        // -----------------------------------------------------------------
        // Results
        // -----------------------------------------------------------------
        results_.value = V[0];

        Real delta2u = (p2u - p2m) / (s2u - s2m);
        Real delta2d = (p2m - p2d) / (s2m - s2d);
        results_.gamma = (delta2u - delta2d) / ((s2u - s2d) / 2.0);

        results_.delta = (p1u - p1d) / (s1u - s1d);

        results_.theta = blackScholesTheta(process_,
                                           results_.value,
                                           results_.delta,
                                           results_.gamma);
    }

}
