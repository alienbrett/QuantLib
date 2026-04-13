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
#include <algorithm>

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

    Real Gaussian1dCallableRangeAccrualEngine::rangeProbability(
        const Date& observationDate,
        const Period& tenor,
        const Rate lower,
        const Rate upper,
        const Date& referenceDate,
        const Real y) const {

        Time obsTime = model_->termStructure()->timeFromReference(observationDate);
        Time refTime = model_->termStructure()->timeFromReference(referenceDate);

        // If observation is in the past or at reference date, use fixing
        if (obsTime <= refTime + QL_EPSILON) {
            Real rate = model_->swapRate(observationDate, tenor,
                                         referenceDate, y,
                                         arguments_.observationIndex);
            return (rate >= lower && rate <= upper) ? 1.0 : 0.0;
        }

        // Inner quadrature: integrate over conditional distribution at obsTime
        // given state y at refTime.
        Array innerZ = model_->yGrid(innerStddevs_, innerPoints_);
        Array innerY = model_->yGrid(innerStddevs_, innerPoints_,
                                      obsTime, refTime, y);

        // Evaluate indicator at each inner grid point
        Array indicator(innerZ.size());
        for (Size i = 0; i < innerZ.size(); ++i) {
            Real rate = model_->swapRate(observationDate, tenor,
                                         observationDate, innerY[i],
                                         arguments_.observationIndex);
            indicator[i] = (rate >= lower && rate <= upper) ? 1.0 : 0.0;
        }

        // Integrate indicator against standard normal density using
        // cubic interpolation + Gaussian polynomial integration
        // (same technique as the main backward induction loop)
        CubicInterpolation interp(
            innerZ.begin(), innerZ.end(), indicator.begin(),
            CubicInterpolation::Spline, true,
            CubicInterpolation::Lagrange, 0.0,
            CubicInterpolation::Lagrange, 0.0);

        Real prob = 0.0;
        for (Size i = 0; i < innerZ.size() - 1; ++i) {
            prob += Gaussian1dModel::gaussianShiftedPolynomialIntegral(
                0.0, interp.cCoefficients()[i],
                interp.bCoefficients()[i],
                interp.aCoefficients()[i],
                indicator[i], innerZ[i],
                innerZ[i], innerZ[i + 1]);
        }
        // Flat extrapolation of indicator in tails
        prob += Gaussian1dModel::gaussianShiftedPolynomialIntegral(
            0.0, 0.0, 0.0, 0.0,
            indicator[innerZ.size() - 1], innerZ[innerZ.size() - 1],
            innerZ[innerZ.size() - 1], 100.0);
        prob += Gaussian1dModel::gaussianShiftedPolynomialIntegral(
            0.0, 0.0, 0.0, 0.0,
            indicator[0], innerZ[0],
            -100.0, innerZ[0]);

        return std::max(0.0, std::min(1.0, prob));
    }

    Real Gaussian1dCallableRangeAccrualEngine::expectedAccrualFraction(
        const Size couponIndex,
        const Date& referenceDate,
        const Real y) const {

        const auto& obsDates = arguments_.observationDates[couponIndex];
        const Period& tenor = arguments_.observationIndex->tenor();
        Rate lower = arguments_.lowerTriggers[couponIndex];
        Rate upper = arguments_.upperTriggers[couponIndex];

        Real accrual = 0.0;
        for (const auto& obsDate : obsDates) {
            accrual += rangeProbability(obsDate, tenor, lower, upper,
                                        referenceDate, y);
        }
        return accrual / static_cast<Real>(obsDates.size());
    }

    Real Gaussian1dCallableRangeAccrualEngine::underlyingNpv(
        const Date& expiry,
        const Real y) const {

        Real type = (Real)arguments_.type;

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
            Real accrualFrac = expectedAccrualFraction(i, expiry, y);
            Real couponRate = arguments_.raGearings[i] * accrualFrac +
                              arguments_.raSpreads[i];
            // Note: accrualFrac modulates the gearing (the "in-range" fraction),
            // and the spread accrues unconditionally.  Coupon amount is:
            // nominal * accrualTime * (gearing * swapRate_avg * accrualFrac + spread)
            // But the standard RA payoff is:
            // nominal * accrualTime * (gearing * refRate + spread) * accrualFrac
            // We use the latter: everything scales by the in-range fraction.
            Real amount = arguments_.raNominal[i] *
                          arguments_.raAccrualTimes[i] *
                          couponRate * accrualFrac;

            // Simpler formulation: the RA coupon pays
            // nominal * dcf * (gearing * swap_rate_at_fixing + spread) * accrualFraction
            // The swap rate at fixing is a single rate, not averaged over obs dates.
            // But for range accruals the coupon pays a fixed spread (no index rate
            // in the coupon itself — the index only determines whether you accrue).
            // Standard RA coupon: nominal * dcf * (gearing * fixedCouponRate + spread) * accrualFrac
            // where accrualFrac = (1/N) * sum(1{L <= index(t_i) <= U})
            //
            // Recompute: the coupon amount is simply:
            amount = arguments_.raNominal[i] *
                     arguments_.raAccrualTimes[i] *
                     (arguments_.raGearings[i] * accrualFrac +
                      arguments_.raSpreads[i]);

            raNpv += amount *
                     model_->zerobond(arguments_.raPayDates[i],
                                      expiry, y, discountCurve_);
        }

        // type: Receiver receives RA, pays fixed
        //       => value = raNpv - fixedNpv  (for Receiver)
        //       Payer pays RA, receives fixed
        //       => value = fixedNpv - raNpv  (for Payer)
        // Convention: type = +1 for Payer, -1 for Receiver
        // The callable option holder can cancel: they exercise when
        // continuation value < 0 (they're losing money on the swap).
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
                            // Extrapolate using last/first cubic segment
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
                        underlyingNpv(expiry0, z[k]) /
                        model_->numeraire(expiry0Time, z[k], discountCurve_);

                    // The callable holder exercises when it is beneficial to
                    // cancel the swap.  If the swap NPV (from the holder's
                    // perspective) is negative, they want to cancel => exercise
                    // value is the avoided future loss.  The continuation value
                    // (npv0[k]) represents the value of holding the option.
                    //
                    // For a callable RA swap, the issuer (short the option)
                    // can call when the swap is in their favor.  The option
                    // value to the holder is:
                    //   max(continuation_value, termination_value)
                    // where termination_value = 0 (swap ceases, no further
                    // payments).  So:
                    //   npv0[k] = max(npv0[k], 0)
                    // But we also need the swap value from the non-exercise
                    // perspective.  The standard approach (matching the
                    // NonstandardSwaption engine):
                    //   exerciseValue = swap_npv / numeraire
                    //   npv0[k] = max(continuation, exerciseValue)
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
