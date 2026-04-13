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

        // Local vol surface from the process
        auto lvSurface = process_->localVolatility().currentLink();

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
        // Trinomial grid geometry.
        // Log-price grid: x_j = ln(S0) + j * dx, j = -N ... +N
        // dx = sigmaRef * sqrt(3 * dt) — standard choice ensuring
        // probabilities stay in (0,1) when local vol ≈ sigmaRef.
        // At step k there are 2k+1 nodes (j = -k ... +k).
        // -----------------------------------------------------------------
        Size N = timeSteps_;
        Time dt = T / N;
        Real dx = sigmaRef * std::sqrt(3.0 * dt);

        QL_ENSURE(dx > 0.0, "grid spacing dx must be positive "
                  "(sigma=" << sigmaRef << ", dt=" << dt << ")");

        Real lnS0 = std::log(S0);

        // Stock price at node (step k, index j) where j in [-k, k]
        // Stored with offset: array index = j + k
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
        // Terminal payoff: step N has 2N+1 nodes, j = -N ... +N
        // Array index = j + N
        // -----------------------------------------------------------------
        Size termNodes = 2 * N + 1;
        Array V(termNodes);
        for (Integer j = -(Integer)N; j <= (Integer)N; ++j)
            V[j + N] = (*payoff)(stockPrice(j));

        // -----------------------------------------------------------------
        // Backward induction with per-node local vol probabilities
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
            Size nNodes = 2 * si + 1;  // nodes at step si: j = -si ... +si

            Real disc = discStep[si];
            Real growth = growthStep[si];
            Time dtk = dtStep[si];
            Time tMid = 0.5 * (stepTimes[si] + stepTimes[si + 1]);

            // ln(growth) = (r-q)*dt from term structures.
            // Log-space drift at each node needs Ito correction:
            //   mu(S)*dt = ln(growth) - 0.5*σ_local(t,S)²*dt
            // Applied per-node below.
            Real lnGrowthRaw = std::log(growth);

            Array newV(nNodes);

            for (Integer j = -(Integer)si; j <= (Integer)si; ++j) {
                Size idx = j + si;  // array index in newV
                Real Sj = stockPrice(j);

                // Query local vol at this node
                Volatility sigLoc;
                try {
                    sigLoc = lvSurface->localVol(tMid, Sj, true);
                } catch (...) {
                    sigLoc = sigmaRef;
                }

                // Trinomial probabilities matching forward + variance in log-space:
                //   x transitions to x+dx (up), x (mid), x-dx (down)
                //   mu = ln(growth) - 0.5*sigLoc^2*dt  (Ito-corrected drift)
                //   E[Δx] = mu*dt     => pu*dx - pd*dx = mu*dt
                //   Var[Δx] = sig^2*dt => pu*dx^2 + pd*dx^2 = sig^2*dt + (mu*dt)^2
                //   pu + pm + pd = 1
                //
                // Solving:
                //   pu = (sig^2*dt + mu² + mu*dx) / (2*dx^2)
                //   pd = (sig^2*dt + mu² - mu*dx) / (2*dx^2)
                //   pm = 1 - pu - pd
                Real localVar = sigLoc * sigLoc * dtk;
                Real mu = lnGrowthRaw - 0.5 * sigLoc * sigLoc * dtk;
                Real dx2 = dx * dx;
                Real mu2 = mu * mu;

                Real pu = (localVar + mu2 + mu * dx) / (2.0 * dx2);
                Real pd = (localVar + mu2 - mu * dx) / (2.0 * dx2);
                Real pm = 1.0 - pu - pd;

                // Clamp probabilities to avoid negative values at extreme nodes
                pu = std::max(0.0, std::min(1.0, pu));
                pd = std::max(0.0, std::min(1.0, pd));
                pm = std::max(0.0, 1.0 - pu - pd);
                // Renormalize after clamping
                Real psum = pu + pm + pd;
                if (psum > 0.0) {
                    pu /= psum; pm /= psum; pd /= psum;
                } else {
                    pu = pd = 0.0; pm = 1.0;
                }

                // V at step (si+1) has 2*(si+1)+1 nodes, index = j + (si+1)
                // Transitions: j+1 (up), j (mid), j-1 (down) in the (si+1) grid
                Size idxUp  = (j + 1) + (si + 1);
                Size idxMid = j + (si + 1);
                Size idxDn  = (j - 1) + (si + 1);

                newV[idx] = disc * (pu * V[idxUp] + pm * V[idxMid] + pd * V[idxDn]);
            }

            // --- VN interpolation at ex-dividend dates ---
            if (divIdx > 0 && divStep[divIdx - 1] == si) {
                --divIdx;
                Real D = divAmt[divIdx];

                // Stock prices at this step (monotonically increasing)
                std::vector<Real> prices(nNodes);
                for (Integer j = -(Integer)si; j <= (Integer)si; ++j)
                    prices[j + si] = stockPrice(j);

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
                    for (Integer j = -(Integer)si; j <= (Integer)si; ++j) {
                        Size idx = j + si;
                        newV[idx] = std::max(newV[idx],
                                             (*payoff)(stockPrice(j)));
                    }
                }
            }

            // --- capture nodes for Greeks ---
            if (si == 2) {
                // j = -2, 0, +2 (every other node for gamma)
                p2d = newV[0]; p2m = newV[2]; p2u = newV[4];
                s2d = stockPrice(-2);
                s2m = stockPrice(0);
                s2u = stockPrice(2);
            }
            if (si == 1) {
                p1d = newV[0]; p1u = newV[2];
                s1d = stockPrice(-1);
                s1u = stockPrice(1);
            }

            V = newV;
        }

        // -----------------------------------------------------------------
        // Results — V has 1 element (the root node)
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
