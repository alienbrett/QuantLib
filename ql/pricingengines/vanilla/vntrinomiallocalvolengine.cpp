/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

#include <ql/pricingengines/vanilla/vntrinomiallocalvolengine.hpp>
#include <ql/exercise.hpp>
#include <ql/pricingengines/greeks.hpp>
#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    VNTrinomialLocalVolEngine::VNTrinomialLocalVolEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size timeSteps)
    : process_(std::move(process)),
      dividends_(std::move(dividends)),
      timeSteps_(timeSteps) {
        QL_REQUIRE(timeSteps >= 3,
                   "at least 3 time steps required, "
                   << timeSteps << " provided");
        registerWith(process_);
    }

    VNTrinomialLocalVolEngine::VNTrinomialLocalVolEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size timeSteps,
        ext::shared_ptr<LocalVolTermStructure> localVol)
    : process_(std::move(process)),
      dividends_(std::move(dividends)),
      timeSteps_(timeSteps),
      explicitLocalVol_(std::move(localVol)) {
        QL_REQUIRE(timeSteps >= 3,
                   "at least 3 time steps required, "
                   << timeSteps << " provided");
        registerWith(process_);
    }

    void VNTrinomialLocalVolEngine::calculate() const {

        DayCounter rfdc  = process_->riskFreeRate()->dayCounter();
        DayCounter divdc = process_->dividendYield()->dayCounter();

        Real S0 = process_->stateVariable()->value();
        QL_REQUIRE(S0 > 0.0, "negative or null underlying given");

        Date referenceDate = process_->riskFreeRate()->referenceDate();
        Date maturityDate  = arguments_.exercise->lastDate();
        Time T = rfdc.yearFraction(referenceDate, maturityDate);

        // Reference vol for grid spacing — ATM implied vol at maturity
        Volatility sigmaRef = process_->blackVolatility()->blackVol(
            maturityDate, S0);

        auto payoff =
            ext::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");
        Real K = payoff->strike();
        QL_REQUIRE(K > 0.0, "strike must be positive");

        bool isAmerican =
            (arguments_.exercise->type() == Exercise::American);

        Time earliestExercise = 0.0;
        if (arguments_.exercise->dates().size() > 1)
            earliestExercise = rfdc.yearFraction(
                referenceDate, arguments_.exercise->date(0));

        // Local vol surface: use explicit surface if provided (fast path
        // for analytic local vol, e.g. EssviLocalVolSurface), otherwise
        // fall back to process's generic Dupire extraction.
        auto lvSurface = explicitLocalVol_
            ? explicitLocalVol_
            : process_->localVolatility().currentLink();

        // -----------------------------------------------------------------
        // Collect discrete cash dividends
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
        // Trinomial grid geometry with adaptive branching.
        //
        // Log-price grid: x_j = ln(S0) + j * dx
        // dx = σ_ref * √(3dt) — fine ATM spacing for best accuracy
        //   near the strike.
        //
        // Standard trinomial: branch j → {j-1, j, j+1}.
        // Problem: when σ_loc > √3·σ_ref, the mid-probability pm goes
        // negative.  Clamping pm destroys variance matching → systematic
        // bias that doesn't converge away with more steps.
        //
        // Fix: adaptive branching.  At each node, compute stretch factor
        //   k = max(1, ceil(√(σ²dt / dx²)))
        // and branch j → {j-k, j, j+k}.  With effective spacing a=k·dx:
        //   pu = (σ²dt + μ² + μa) / (2a²)
        //   pd = (σ²dt + μ² - μa) / (2a²)
        //   pm = 1 - (σ²dt + μ²) / a²
        // which is always ≥ 0 by construction.
        //
        // Tree width at step i = i·K where K is the global max stretch.
        // K is determined by probing the local vol surface at startup.
        // The tree stays recombining — all children are on the same
        // uniform dx grid.
        // -----------------------------------------------------------------
        Size N = timeSteps_;
        Time dt = T / N;
        Real dx = sigmaRef * std::sqrt(3.0 * dt);

        QL_ENSURE(dx > 0.0, "grid spacing dx must be positive "
                  "(sigma=" << sigmaRef << ", dt=" << dt << ")");

        // Probe local vol to determine global max stretch K.
        Volatility sigMax = sigmaRef;
        {
            Real range = 5.0 * sigmaRef * std::sqrt(T);
            Real logProbes[] = {-range, -0.7*range, -0.4*range, 0.0,
                                 0.4*range, 0.7*range, range};
            Time timeProbes[] = {0.1*T, 0.3*T, 0.6*T, 0.9*T};
            for (Time tp : timeProbes) {
                if (tp < 1e-6) tp = 1e-6;
                for (Real lp : logProbes) {
                    Real Sp = S0 * std::exp(lp);
                    if (Sp <= 0.0) continue;
                    try {
                        Volatility v = lvSurface->localVol(tp, Sp, true);
                        sigMax = std::max(sigMax, v);
                    } catch (...) {}
                }
            }
            for (Time tp : timeProbes) {
                if (tp < 1e-6) tp = 1e-6;
                try {
                    Volatility v = lvSurface->localVol(tp, K, true);
                    sigMax = std::max(sigMax, v);
                } catch (...) {}
            }
        }

        Real lnS0 = std::log(S0);

        auto logPrice = [&](Integer j) -> Real {
            return lnS0 + j * dx;
        };
        auto stockPrice = [&](Integer j) -> Real {
            return std::exp(logPrice(j));
        };

        // -----------------------------------------------------------------
        // Non-uniform time grid with ex-div date alignment
        // (same logic as VNBinomialVanillaEngine)
        // -----------------------------------------------------------------
        std::vector<Time> stepTimes(N + 1);
        std::vector<Size> divStep;
        std::vector<Real> divAmt;

        if (cashDivs.empty()) {
            for (Size k = 0; k <= N; ++k)
                stepTimes[k] = k * dt;
        } else {
            std::vector<Time> breaks;
            breaks.push_back(0.0);
            for (const auto& d : cashDivs)
                if (d.t > 1e-10 && d.t < T - 1e-10)
                    breaks.push_back(d.t);
            breaks.push_back(T);

            Size nSeg = breaks.size() - 1;
            std::vector<Size> segSteps(nSeg);
            Size minPerSeg = 2;
            Size reserved = minPerSeg * nSeg;
            QL_REQUIRE(N >= reserved,
                       "VN trinomial: need at least " << reserved
                       << " steps for " << nSeg << " segments, got " << N);
            Size pool = N - reserved;
            Size allocated = 0;
            for (Size s = 0; s < nSeg; ++s) {
                Real frac = (breaks[s+1] - breaks[s]) / T;
                Size alloc = static_cast<Size>(std::round(pool * frac));
                segSteps[s] = minPerSeg + alloc;
                allocated += alloc;
            }
            Size longest = 0;
            for (Size s = 1; s < nSeg; ++s)
                if (breaks[s+1] - breaks[s] > breaks[longest+1] - breaks[longest])
                    longest = s;
            if (allocated < pool)
                segSteps[longest] += (pool - allocated);
            else if (allocated > pool && segSteps[longest] > minPerSeg + (allocated - pool))
                segSteps[longest] -= (allocated - pool);

            Size idx = 0;
            stepTimes[0] = 0.0;
            for (Size s = 0; s < nSeg; ++s) {
                Time segStart = breaks[s];
                Time segLen = breaks[s+1] - breaks[s];
                Time segDt = segLen / segSteps[s];
                for (Size k = 1; k <= segSteps[s]; ++k)
                    stepTimes[++idx] = segStart + k * segDt;
            }

            for (Size di = 0; di < cashDivs.size(); ++di) {
                Time td = cashDivs[di].t;
                auto it = std::lower_bound(
                    stepTimes.begin(), stepTimes.end(), td - 1e-10);
                Size s = static_cast<Size>(
                    std::distance(stepTimes.begin(), it));
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

        // Recompute dt per step from the (possibly non-uniform) grid
        std::vector<Time> dtStep(N);
        for (Size k = 0; k < N; ++k)
            dtStep[k] = stepTimes[k + 1] - stepTimes[k];

        // -----------------------------------------------------------------
        // Per-step discount factors from actual term structures
        // -----------------------------------------------------------------
        std::vector<Real> dfRf(N + 1), dfQ(N + 1);
        dfRf[0] = dfQ[0] = 1.0;
        for (Size k = 1; k <= N; ++k) {
            dfRf[k] = process_->riskFreeRate()->discount(stepTimes[k]);
            dfQ[k]  = process_->dividendYield()->discount(stepTimes[k]);
        }

        std::vector<Real> discStep(N), growthStep(N);
        for (Size k = 0; k < N; ++k) {
            discStep[k] = dfRf[k + 1] / dfRf[k];
            growthStep[k] = (dfQ[k + 1] * dfRf[k]) / (dfQ[k] * dfRf[k + 1]);
        }

        // -----------------------------------------------------------------
        // Adaptive branching: global max stretch factor KK.
        //
        // KK determines tree width at step i: nodes from -i*KK to +i*KK.
        // Computed from max local vol across the pricing-relevant range
        // and the worst-case dt (non-uniform grid may have longer steps).
        // +1 margin for probing gaps and drift contribution.
        // -----------------------------------------------------------------
        Real maxDt = *std::max_element(dtStep.begin(), dtStep.end());
        Integer KK = std::max(Integer(1),
            (Integer)std::ceil(
                std::sqrt(sigMax * sigMax * maxDt / (dx * dx))));
        KK += 1;  // safety margin

        // -----------------------------------------------------------------
        // Terminal payoff: step N has nodes j = -N*KK ... +N*KK
        // -----------------------------------------------------------------
        Integer termHalf = static_cast<Integer>(N) * KK;
        Size termNodes = 2 * termHalf + 1;
        Array V(termNodes);
        for (Integer j = -termHalf; j <= termHalf; ++j)
            V[j + termHalf] = (*payoff)(stockPrice(j));

        // -----------------------------------------------------------------
        // Backward induction with per-node adaptive branching
        // and VN dividend interpolation
        // -----------------------------------------------------------------
        Size divIdx = nDiv;

        // Greeks capture: values at steps 1 and 2
        Real p2d = 0, p2m = 0, p2u = 0;
        Real s2d = 0, s2m = 0, s2u = 0;
        Real p1d = 0, p1u = 0;
        Real s1d = 0, s1u = 0;

        for (Integer i = N - 1; i >= 0; --i) {
            Size si = static_cast<Size>(i);
            Integer halfWidth = static_cast<Integer>(si) * KK;
            Size nNodes = 2 * halfWidth + 1;
            Integer childHalf = static_cast<Integer>(si + 1) * KK;

            Real disc = discStep[si];
            Real growth = growthStep[si];
            Time dtk = dtStep[si];
            Time tMid = 0.5 * (stepTimes[si] + stepTimes[si + 1]);

            Real lnGrowthRaw = std::log(growth);

            Array newV(nNodes);

            for (Integer j = -halfWidth; j <= halfWidth; ++j) {
                Size idx = j + halfWidth;
                Real Sj = stockPrice(j);

                // Query local vol at this node
                Volatility sigLoc;
                try {
                    sigLoc = lvSurface->localVol(tMid, Sj, true);
                } catch (...) {
                    sigLoc = sigmaRef;
                }

                Real localVar = sigLoc * sigLoc * dtk;
                Real mu = lnGrowthRaw - 0.5 * sigLoc * sigLoc * dtk;
                Real mu2 = mu * mu;
                Real varTerm = localVar + mu2;

                // Adaptive stretch: choose kk so pm ≥ 0
                // pm = 1 - varTerm/(kk²dx²) ≥ 0  =>  kk ≥ √(varTerm/dx²)
                Integer kk = 1;
                Real dx2 = dx * dx;
                if (varTerm > dx2) {
                    kk = (Integer)std::ceil(std::sqrt(varTerm / dx2));
                }
                // Cap at global max (shouldn't bind if probe was accurate)
                kk = std::min(kk, KK);

                Real a = kk * dx;
                Real a2 = a * a;

                Real pu = (varTerm + mu * a) / (2.0 * a2);
                Real pd = (varTerm - mu * a) / (2.0 * a2);
                Real pm = 1.0 - pu - pd;

                // Safety: pm should be ≥ 0 by construction.  If floating
                // point puts it slightly negative, nudge.
                if (pm < 0.0) {
                    pm = 0.0;
                    Real psum = pu + pd;
                    if (psum > 0.0) { pu /= psum; pd /= psum; }
                }

                // Child node indices in step si+1 array
                Integer jUp  = j + kk;
                Integer jMid = j;
                Integer jDn  = j - kk;

                // Clamp to child bounds (only at very edge of tree)
                jUp  = std::min(jUp, childHalf);
                jDn  = std::max(jDn, -childHalf);

                Size idxUp  = jUp  + childHalf;
                Size idxMid = jMid + childHalf;
                Size idxDn  = jDn  + childHalf;

                newV[idx] = disc * (pu * V[idxUp]
                                  + pm * V[idxMid]
                                  + pd * V[idxDn]);
            }

            // --- VN interpolation at ex-dividend dates ---
            if (divIdx > 0 && divStep[divIdx - 1] == si) {
                --divIdx;
                Real D = divAmt[divIdx];

                std::vector<Real> prices(nNodes);
                for (Integer j = -halfWidth; j <= halfWidth; ++j)
                    prices[j + halfWidth] = stockPrice(j);

                Array adjV(nNodes);
                for (Size idx = 0; idx < nNodes; ++idx) {
                    Real Spost = prices[idx] - D;
                    if (Spost <= 0.0) {
                        adjV[idx] = (*payoff)(std::max(Spost, 0.0));
                    } else {
                        auto it = std::lower_bound(
                            prices.begin(), prices.end(), Spost);
                        Size j1, j0;
                        if (it == prices.begin()) {
                            j0 = 0; j1 = 1;
                        } else if (it == prices.end()) {
                            j1 = nNodes - 1; j0 = nNodes - 2;
                        } else {
                            j1 = static_cast<Size>(
                                std::distance(prices.begin(), it));
                            j0 = j1 - 1;
                        }
                        Real w = (Spost - prices[j0])
                               / (prices[j1] - prices[j0]);
                        adjV[idx] = newV[j0] + w * (newV[j1] - newV[j0]);
                    }
                }
                newV = adjV;
            }

            // --- early exercise ---
            if (isAmerican) {
                Time stepTime = stepTimes[si];
                if (stepTime >= earliestExercise) {
                    for (Integer j = -halfWidth; j <= halfWidth; ++j) {
                        Size idx = j + halfWidth;
                        newV[idx] = std::max(newV[idx],
                                             (*payoff)(stockPrice(j)));
                    }
                }
            }

            // --- capture nodes for Greeks ---
            // j = -2, 0, +2 always exist at step 2 (halfWidth = 2*KK ≥ 2)
            if (si == 2) {
                p2d = newV[(-2) + halfWidth];
                p2m = newV[0 + halfWidth];
                p2u = newV[2 + halfWidth];
                s2d = stockPrice(-2);
                s2m = stockPrice(0);
                s2u = stockPrice(2);
            }
            // j = -1, +1 always exist at step 1 (halfWidth = KK ≥ 1)
            if (si == 1) {
                p1d = newV[(-1) + halfWidth];
                p1u = newV[1 + halfWidth];
                s1d = stockPrice(-1);
                s1u = stockPrice(1);
            }

            V = newV;
        }

        // -----------------------------------------------------------------
        // Results — V has 1 element (root: halfWidth = 0*KK = 0)
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
