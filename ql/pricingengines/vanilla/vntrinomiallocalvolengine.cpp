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
        Size timeSteps,
        std::vector<Time> breakTimes)
    : process_(std::move(process)),
      dividends_(std::move(dividends)),
      timeSteps_(timeSteps),
      breakTimes_(std::move(breakTimes)) {
        QL_REQUIRE(timeSteps >= 3,
                   "at least 3 time steps required, "
                   << timeSteps << " provided");
        registerWith(process_);
    }

    VNTrinomialLocalVolEngine::VNTrinomialLocalVolEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size timeSteps,
        ext::shared_ptr<LocalVolTermStructure> localVol,
        std::vector<Time> breakTimes)
    : process_(std::move(process)),
      dividends_(std::move(dividends)),
      timeSteps_(timeSteps),
      explicitLocalVol_(std::move(localVol)),
      breakTimes_(std::move(breakTimes)) {
        QL_REQUIRE(timeSteps >= 3,
                   "at least 3 time steps required, "
                   << timeSteps << " provided");
        registerWith(process_);
    }

    VNTrinomialLocalVolEngine::VNTrinomialLocalVolEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        DividendSchedule dividends,
        Size timeSteps,
        ext::shared_ptr<LocalVolTermStructure> localVol,
        Size lvGridStride,
        std::vector<Time> breakTimes)
    : process_(std::move(process)),
      dividends_(std::move(dividends)),
      timeSteps_(timeSteps),
      explicitLocalVol_(std::move(localVol)),
      lvGridStride_(lvGridStride),
      breakTimes_(std::move(breakTimes)) {
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
        // Non-uniform time grid with break-time alignment.
        //
        // Break times come from two sources:
        //   1. Ex-dividend dates (cash dividends → VN interpolation)
        //   2. User-supplied break times (e.g. vol surface slice dates)
        // Aligning step boundaries with these times prevents local vol
        // discontinuities from causing non-monotonic convergence.
        // -----------------------------------------------------------------
        std::vector<Time> stepTimes(N + 1);
        std::vector<Size> divStep;
        std::vector<Real> divAmt;

        // Collect all break times (dividends + explicit)
        std::vector<Time> allBreaks;
        for (const auto& d : cashDivs)
            if (d.t > 1e-10 && d.t < T - 1e-10)
                allBreaks.push_back(d.t);
        for (Time bt : breakTimes_)
            if (bt > 1e-10 && bt < T - 1e-10)
                allBreaks.push_back(bt);
        std::sort(allBreaks.begin(), allBreaks.end());
        allBreaks.erase(std::unique(allBreaks.begin(), allBreaks.end(),
            [](Time a, Time b) { return std::abs(a - b) < 1e-10; }),
            allBreaks.end());

        if (allBreaks.empty()) {
            for (Size k = 0; k <= N; ++k)
                stepTimes[k] = k * dt;
        } else {
            std::vector<Time> breaks;
            breaks.push_back(0.0);
            for (Time bt : allBreaks)
                breaks.push_back(bt);
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

        std::vector<Time> dtStep(N);
        for (Size k = 0; k < N; ++k)
            dtStep[k] = stepTimes[k + 1] - stepTimes[k];

        // -----------------------------------------------------------------
        // Per-step discount factors
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
        // -----------------------------------------------------------------
        Real maxDt = *std::max_element(dtStep.begin(), dtStep.end());
        Integer KK = std::max(Integer(1),
            (Integer)std::ceil(
                std::sqrt(sigMax * sigMax * maxDt / (dx * dx))));
        KK += 1;  // safety margin

        // -----------------------------------------------------------------
        // Coarse-grid local vol precomputation (when lvGridStride_ > 0).
        //
        // For each time step, evaluate local vol at every stride-th
        // tree node and store.  During backward induction, linearly
        // interpolate between the two nearest grid points.
        //
        // Cost: O(KN²/S) evaluations vs O(KN²) exact, where S=stride.
        // Interpolation error is negligible for smooth surfaces (eSSVI).
        // -----------------------------------------------------------------
        Size stride = lvGridStride_;
        bool useGrid = (stride > 0);

        // Per-step grid: lvGrid_[si] has local vols at j = jMin, jMin+stride, ...
        // lvGridJMin_[si] = first j in the grid for step si
        // lvGridSize_[si] = number of grid points for step si
        std::vector<std::vector<Volatility>> lvGrid_;
        std::vector<Integer> lvGridJMin_;

        if (useGrid) {
            lvGrid_.resize(N);
            lvGridJMin_.resize(N);
            Integer iStride = static_cast<Integer>(stride);
            for (Size si = 0; si < N; ++si) {
                Integer halfW = static_cast<Integer>(si) * KK;
                Time tMid = 0.5 * (stepTimes[si] + stepTimes[si + 1]);

                // Align grid to stride boundaries so spacing is truly
                // uniform — required for correct Lagrange interpolation.
                // Round halfW up to next multiple of stride.
                Integer halfAligned = ((halfW + iStride - 1) / iStride) * iStride;
                Integer jMin = -halfAligned;
                Integer jMax =  halfAligned;
                lvGridJMin_[si] = jMin;

                Size nGrid = static_cast<Size>((jMax - jMin) / iStride) + 1;
                lvGrid_[si].resize(nGrid);

                for (Size g = 0; g < nGrid; ++g) {
                    Integer j = jMin + static_cast<Integer>(g) * iStride;
                    Real Sj = stockPrice(j);
                    try {
                        lvGrid_[si][g] = lvSurface->localVol(tMid, Sj, true);
                    } catch (...) {
                        lvGrid_[si][g] = sigmaRef;
                    }
                }
            }
        }

        // Local vol lookup: cubic Lagrange interpolation or exact evaluation.
        // Cubic (4-point) gives O(h⁴) error vs O(h²) for linear, critical
        // for eSSVI surfaces with non-trivial curvature.  Falls back to
        // linear when the per-step grid has < 4 points (early time steps).
        auto getLocalVol = [&](Size si, Integer j, Time tMid) -> Volatility {
            if (useGrid) {
                Integer jMin = lvGridJMin_[si];
                const auto& grid = lvGrid_[si];
                Size n = grid.size();
                Real fIdx = static_cast<Real>(j - jMin)
                          / static_cast<Real>(stride);
                if (fIdx < 0.0) fIdx = 0.0;
                if (fIdx > static_cast<Real>(n - 1))
                    fIdx = static_cast<Real>(n - 1);

                if (n >= 4) {
                    // Cubic Lagrange interpolation over 4 nearest points.
                    // Choose i0 so points i0..i0+3 straddle fIdx.
                    Integer i0 = static_cast<Integer>(fIdx) - 1;
                    if (i0 < 0) i0 = 0;
                    if (i0 > static_cast<Integer>(n) - 4)
                        i0 = static_cast<Integer>(n) - 4;
                    Real x = fIdx - static_cast<Real>(i0);
                    // x ∈ [0, 3], interpolation points at x=0,1,2,3
                    Real x0 = x, x1 = x - 1.0, x2 = x - 2.0, x3 = x - 3.0;
                    Real L0 = (-x1 * x2 * x3) / 6.0;
                    Real L1 = ( x0 * x2 * x3) / 2.0;
                    Real L2 = (-x0 * x1 * x3) / 2.0;
                    Real L3 = ( x0 * x1 * x2) / 6.0;
                    return L0 * grid[i0] + L1 * grid[i0+1]
                         + L2 * grid[i0+2] + L3 * grid[i0+3];
                } else {
                    // Linear fallback for small grids
                    Size lo = static_cast<Size>(fIdx);
                    if (lo >= n - 1) lo = (n > 1) ? n - 2 : 0;
                    Real w = fIdx - static_cast<Real>(lo);
                    return grid[lo] + w * (grid[lo + 1] - grid[lo]);
                }
            } else {
                Real Sj = stockPrice(j);
                try {
                    return lvSurface->localVol(tMid, Sj, true);
                } catch (...) {
                    return sigmaRef;
                }
            }
        };

        // -----------------------------------------------------------------
        // Terminal payoff
        // -----------------------------------------------------------------
        Integer termHalf = static_cast<Integer>(N) * KK;
        Size termNodes = 2 * termHalf + 1;
        Array V(termNodes);
        for (Integer j = -termHalf; j <= termHalf; ++j)
            V[j + termHalf] = (*payoff)(stockPrice(j));

        // -----------------------------------------------------------------
        // Backward induction
        // -----------------------------------------------------------------
        Size divIdx = nDiv;

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

                Volatility sigLoc = getLocalVol(si, j, tMid);

                Real localVar = sigLoc * sigLoc * dtk;
                Real mu = lnGrowthRaw - 0.5 * sigLoc * sigLoc * dtk;
                Real mu2 = mu * mu;
                Real varTerm = localVar + mu2;

                Integer kk = 1;
                Real dx2 = dx * dx;
                if (varTerm > dx2) {
                    kk = (Integer)std::ceil(std::sqrt(varTerm / dx2));
                }
                kk = std::min(kk, KK);

                Real a = kk * dx;
                Real a2 = a * a;

                Real pu = (varTerm + mu * a) / (2.0 * a2);
                Real pd = (varTerm - mu * a) / (2.0 * a2);
                Real pm = 1.0 - pu - pd;

                if (pm < 0.0) {
                    pm = 0.0;
                    Real psum = pu + pd;
                    if (psum > 0.0) { pu /= psum; pd /= psum; }
                }

                Integer jUp  = std::min(j + kk, childHalf);
                Integer jDn  = std::max(j - kk, -childHalf);

                Size idxUp  = jUp  + childHalf;
                Size idxMid = j    + childHalf;
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
            if (si == 2) {
                p2d = newV[(-2) + halfWidth];
                p2m = newV[0 + halfWidth];
                p2u = newV[2 + halfWidth];
                s2d = stockPrice(-2);
                s2m = stockPrice(0);
                s2u = stockPrice(2);
            }
            if (si == 1) {
                p1d = newV[(-1) + halfWidth];
                p1u = newV[1 + halfWidth];
                s1d = stockPrice(-1);
                s1u = stockPrice(1);
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
