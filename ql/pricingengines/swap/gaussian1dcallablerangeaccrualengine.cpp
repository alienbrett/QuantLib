/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 chloride

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

#include <ql/pricingengines/swap/gaussian1dcallablerangeaccrualengine.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    Gaussian1dCallableRangeAccrualEngine::Gaussian1dCallableRangeAccrualEngine(
        const ext::shared_ptr<Gaussian1dModel>& model,
        const int integrationPoints,
        const Real stddevs,
        const bool extrapolatePayoff,
        const bool flatPayoffExtrapolation,
        const Handle<YieldTermStructure>& discountCurve,
        const int innerPoints,
        const Real innerStddevs)
    : GenericModelEngine<Gaussian1dModel,
                         CallableRangeAccrualSwap::arguments,
                         CallableRangeAccrualSwap::results>(model),
      integrationPoints_(integrationPoints),
      stddevs_(stddevs),
      extrapolatePayoff_(extrapolatePayoff),
      flatPayoffExtrapolation_(flatPayoffExtrapolation),
      discountCurve_(discountCurve),
      innerPoints_(innerPoints),
      innerStddevs_(innerStddevs) {

        if (!discountCurve_.empty())
            registerWith(discountCurve_);
    }

    Gaussian1dCallableRangeAccrualEngine::Gaussian1dCallableRangeAccrualEngine(
        const Handle<Gaussian1dModel>& model,
        const int integrationPoints,
        const Real stddevs,
        const bool extrapolatePayoff,
        const bool flatPayoffExtrapolation,
        const Handle<YieldTermStructure>& discountCurve,
        const int innerPoints,
        const Real innerStddevs)
    : GenericModelEngine<Gaussian1dModel,
                         CallableRangeAccrualSwap::arguments,
                         CallableRangeAccrualSwap::results>(model),
      integrationPoints_(integrationPoints),
      stddevs_(stddevs),
      extrapolatePayoff_(extrapolatePayoff),
      flatPayoffExtrapolation_(flatPayoffExtrapolation),
      discountCurve_(discountCurve),
      innerPoints_(innerPoints),
      innerStddevs_(innerStddevs) {

        if (!discountCurve_.empty())
            registerWith(discountCurve_);
    }

    std::vector<Gaussian1dCallableRangeAccrualEngine::ObsRoot>
    Gaussian1dCallableRangeAccrualEngine::precomputeRoots(
        const Size couponIndex,
        const Date& settlement) const {

        const auto& obsDates = arguments_.observationDates[couponIndex];
        const Period& tenor = arguments_.observationIndex->tenor();
        Rate lower = arguments_.lowerTriggers[couponIndex];
        Rate upper = arguments_.upperTriggers[couponIndex];

        const auto& yts = model_->termStructure();

        // Wide y-range for root finding (unconditional distribution)
        Real ySearchLo = -innerStddevs_;
        Real ySearchHi =  innerStddevs_;

        Brent solver;
        solver.setMaxEvaluations(50);

        std::vector<ObsRoot> roots(obsDates.size());

        for (Size j = 0; j < obsDates.size(); ++j) {
            auto& r = roots[j];
            r.obsTime = yts->timeFromReference(obsDates[j]);
            Time settlementTime = yts->timeFromReference(settlement);

            // Past observation — deterministic
            if (r.obsTime <= settlementTime + QL_EPSILON) {
                r.fixedObs = true;
                r.fixedRate = model_->swapRate(
                    obsDates[j], tenor, settlement, 0.0,
                    arguments_.observationIndex);
                continue;
            }

            // Evaluate swap rate at wide bounds of standardized y space
            Rate rateLo = model_->swapRate(
                obsDates[j], tenor, obsDates[j], ySearchLo,
                arguments_.observationIndex);
            Rate rateHi = model_->swapRate(
                obsDates[j], tenor, obsDates[j], ySearchHi,
                arguments_.observationIndex);

            if (rateLo > upper || rateHi < lower) {
                r.alwaysOut = true;
                continue;
            }
            if (rateLo >= lower && rateHi <= upper) {
                r.alwaysIn = true;
                continue;
            }

            // Find y_L: swapRate(y_L) = lower
            if (rateLo >= lower) {
                r.y_L = ySearchLo;
            } else {
                auto fLower = [&](Real yy) -> Real {
                    return model_->swapRate(
                               obsDates[j], tenor, obsDates[j], yy,
                               arguments_.observationIndex) -
                           lower;
                };
                r.y_L = solver.solve(fLower, 1e-7,
                                     0.5 * (ySearchLo + ySearchHi),
                                     ySearchLo, ySearchHi);
            }

            // Find y_U: swapRate(y_U) = upper
            if (rateHi <= upper) {
                r.y_U = ySearchHi;
            } else {
                auto fUpper = [&](Real yy) -> Real {
                    return model_->swapRate(
                               obsDates[j], tenor, obsDates[j], yy,
                               arguments_.observationIndex) -
                           upper;
                };
                r.y_U = solver.solve(fUpper, 1e-7,
                                     0.5 * (ySearchLo + ySearchHi),
                                     ySearchLo, ySearchHi);
            }
        }

        return roots;
    }

    Real Gaussian1dCallableRangeAccrualEngine::rangeProbabilityFromRoot(
        const ObsRoot& root,
        const Time refTime,
        const Real y) const {

        // Past observation
        if (root.fixedObs)
            return (root.fixedRate >= arguments_.lowerTriggers[0] &&
                    root.fixedRate <= arguments_.upperTriggers[0])
                       ? 1.0
                       : 0.0;

        if (root.alwaysOut)
            return 0.0;
        if (root.alwaysIn)
            return 1.0;

        // Conditional distribution of y_obs given y at refTime
        const auto& proc = model_->stateProcess();

        Real stdDev_0_ref = proc->stdDeviation(0.0, 0.0, refTime);
        Real e_0_ref = proc->expectation(0.0, 0.0, refTime);
        Real x_ref = y * stdDev_0_ref + e_0_ref;

        Real dt = root.obsTime - refTime;
        Real stdDev_ref_obs = proc->stdDeviation(refTime, 0.0, dt);
        Real e_ref_obs = proc->expectation(refTime, x_ref, dt);

        Real stdDev_0_obs = proc->stdDeviation(0.0, 0.0, root.obsTime);
        Real e_0_obs = proc->expectation(0.0, 0.0, root.obsTime);

        Real mu_y = (e_ref_obs - e_0_obs) / stdDev_0_obs;
        Real sigma_y = stdDev_ref_obs / stdDev_0_obs;

        CumulativeNormalDistribution cdf;
        Real z_L = (root.y_L - mu_y) / sigma_y;
        Real z_U = (root.y_U - mu_y) / sigma_y;

        return std::max(0.0, cdf(z_U) - cdf(z_L));
    }

    Real Gaussian1dCallableRangeAccrualEngine::expectedAccrualFraction(
        const std::vector<ObsRoot>& roots,
        const Time refTime,
        const Real y) const {

        Real accrual = 0.0;
        for (const auto& root : roots) {
            accrual += rangeProbabilityFromRoot(root, refTime, y);
        }
        return accrual / static_cast<Real>(roots.size());
    }

    Real Gaussian1dCallableRangeAccrualEngine::underlyingNpv(
        const Date& expiry,
        const Real y,
        const std::vector<std::vector<ObsRoot>>& allRoots) const {

        Time expiryTime =
            model_->termStructure()->timeFromReference(expiry);

        // Fixed leg NPV
        Size fixedIdx =
            std::upper_bound(arguments_.fixedResetDates.begin(),
                             arguments_.fixedResetDates.end(), expiry - 1) -
            arguments_.fixedResetDates.begin();

        Real fixedNpv = 0.0;
        for (Size i = fixedIdx; i < arguments_.fixedResetDates.size(); ++i) {
            fixedNpv += arguments_.fixedCoupons[i] *
                        model_->zerobond(arguments_.fixedPayDates[i],
                                         expiry, y, discountCurve_);
        }

        // Range accrual leg NPV
        Size raIdx =
            std::upper_bound(arguments_.raResetDates.begin(),
                             arguments_.raResetDates.end(), expiry - 1) -
            arguments_.raResetDates.begin();

        Real raNpv = 0.0;
        for (Size i = raIdx; i < arguments_.raResetDates.size(); ++i) {
            Real accrualFrac =
                expectedAccrualFraction(allRoots[i], expiryTime, y);

            // Standard RA coupon:
            //   nominal * dcf * (gearing * accrualFrac + spread)
            Real amount = arguments_.raNominal[i] *
                          arguments_.raAccrualTimes[i] *
                          (arguments_.raGearings[i] * accrualFrac +
                           arguments_.raSpreads[i]);

            raNpv += amount *
                     model_->zerobond(arguments_.raPayDates[i],
                                      expiry, y, discountCurve_);
        }

        Real type = static_cast<Real>(arguments_.type);
        return type * (raNpv - fixedNpv);
    }

    void Gaussian1dCallableRangeAccrualEngine::calculate() const {

        QL_REQUIRE(arguments_.settlementType == Settlement::Physical,
                   "cash settlement not supported");

        Date settlement = model_->termStructure()->referenceDate();

        if (arguments_.exercise->dates().back() <= settlement) {
            results_.value = 0.0;
            return;
        }

        // --- Pre-compute observation roots (O(obs_dates) swap rate evals) ---
        Size nRA = arguments_.raPayDates.size();
        std::vector<std::vector<ObsRoot>> allRoots(nRA);
        for (Size i = 0; i < nRA; ++i) {
            allRoots[i] = precomputeRoots(i, settlement);
        }

        // --- Backward induction ---
        int idx = arguments_.exercise->dates().size() - 1;
        int minIdxAlive = static_cast<int>(
            std::upper_bound(arguments_.exercise->dates().begin(),
                             arguments_.exercise->dates().end(), settlement) -
            arguments_.exercise->dates().begin());

        Array npv0(2 * integrationPoints_ + 1, 0.0),
            npv1(2 * integrationPoints_ + 1, 0.0);
        Array z = model_->yGrid(stddevs_, integrationPoints_);
        Array p(z.size(), 0.0);

        Date expiry1 = Date(), expiry0;
        Time expiry1Time = Null<Real>(), expiry0Time;

        do {

            if (idx == minIdxAlive - 1)
                expiry0 = settlement;
            else
                expiry0 = arguments_.exercise->dates()[idx];

            expiry0Time = std::max(
                model_->termStructure()->timeFromReference(expiry0), 0.0);

            for (Size k = 0; k < (expiry0 > settlement ? npv0.size() : 1);
                 k++) {

                Real price = 0.0;
                if (expiry1Time != Null<Real>()) {
                    Array yg = model_->yGrid(stddevs_, integrationPoints_,
                                              expiry1Time, expiry0Time,
                                              expiry0 > settlement ? z[k]
                                                                   : 0.0);
                    CubicInterpolation payoff0(
                        z.begin(), z.end(), npv1.begin(),
                        CubicInterpolation::Spline, true,
                        CubicInterpolation::Lagrange, 0.0,
                        CubicInterpolation::Lagrange, 0.0);
                    for (Size i = 0; i < yg.size(); i++) {
                        p[i] = payoff0(yg[i], true);
                    }
                    CubicInterpolation payoff1(
                        z.begin(), z.end(), p.begin(),
                        CubicInterpolation::Spline, true,
                        CubicInterpolation::Lagrange, 0.0,
                        CubicInterpolation::Lagrange, 0.0);
                    for (Size i = 0; i < z.size() - 1; i++) {
                        price +=
                            Gaussian1dModel::gaussianShiftedPolynomialIntegral(
                                0.0, payoff1.cCoefficients()[i],
                                payoff1.bCoefficients()[i],
                                payoff1.aCoefficients()[i], p[i], z[i], z[i],
                                z[i + 1]);
                    }
                    if (extrapolatePayoff_) {
                        if (flatPayoffExtrapolation_) {
                            price +=
                                Gaussian1dModel::gaussianShiftedPolynomialIntegral(
                                    0.0, 0.0, 0.0, 0.0, p[z.size() - 2],
                                    z[z.size() - 2], z[z.size() - 1], 100.0);
                            price +=
                                Gaussian1dModel::gaussianShiftedPolynomialIntegral(
                                    0.0, 0.0, 0.0, 0.0, p[0], z[0], -100.0,
                                    z[0]);
                        } else {
                            price +=
                                Gaussian1dModel::gaussianShiftedPolynomialIntegral(
                                    0.0,
                                    payoff1.cCoefficients()[z.size() - 2],
                                    payoff1.bCoefficients()[z.size() - 2],
                                    payoff1.aCoefficients()[z.size() - 2],
                                    p[z.size() - 2], z[z.size() - 2],
                                    z[z.size() - 1], 100.0);
                            price +=
                                Gaussian1dModel::gaussianShiftedPolynomialIntegral(
                                    0.0, payoff1.cCoefficients()[0],
                                    payoff1.bCoefficients()[0],
                                    payoff1.aCoefficients()[0], p[0], z[0],
                                    -100.0, z[0]);
                        }
                    }
                }

                npv0[k] = price;

                // Exercise decision
                if (expiry0 > settlement) {
                    Real exerciseValue =
                        underlyingNpv(expiry0, z[k], allRoots) /
                        model_->numeraire(expiry0Time, z[k], discountCurve_);
                    npv0[k] = std::max(npv0[k], exerciseValue);
                }
            }

            npv1.swap(npv0);

            expiry1 = expiry0;
            expiry1Time = expiry0Time;

        } while (--idx >= minIdxAlive - 1);

        results_.value =
            npv1[0] * model_->numeraire(0.0, 0.0, discountCurve_);
    }
}
