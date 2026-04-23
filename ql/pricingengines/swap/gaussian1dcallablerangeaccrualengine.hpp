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

/*! \file gaussian1dcallablerangeaccrualengine.hpp
    \brief Gaussian 1-d model engine for callable range accrual swaps.

    Designed for use with the MarkovFunctional model, which recovers
    correct CMS forward rates (including convexity adjustment) at each
    grid point via its calibration to the swaption smile.

    Backward induction over Bermudan exercise dates.  Range accrual
    coupon values are computed by evaluating the conditional probability
    that the observation index falls within [lower, upper] at each
    observation date, given the state variable at the exercise date.

    Exploits the monotonicity of swapRate(t, y) in y to replace
    inner quadrature with pre-computed root finding + CDF evaluation.
    Observation roots y_L, y_U are computed once per obs date, then
    conditional probabilities for each grid point are O(1) CDF lookups.
*/

#ifndef quantlib_gaussian1d_callable_range_accrual_engine_hpp
#define quantlib_gaussian1d_callable_range_accrual_engine_hpp

#include <ql/instruments/callablerangeaccrualswap.hpp>
#include <ql/models/shortrate/onefactormodels/gaussian1dmodel.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>

namespace QuantLib {

    class Gaussian1dCallableRangeAccrualEngine
        : public GenericModelEngine<Gaussian1dModel,
                                    CallableRangeAccrualSwap::arguments,
                                    CallableRangeAccrualSwap::results> {
      public:
        Gaussian1dCallableRangeAccrualEngine(
            const ext::shared_ptr<Gaussian1dModel>& model,
            int integrationPoints = 64,
            Real stddevs = 7.0,
            bool extrapolatePayoff = true,
            bool flatPayoffExtrapolation = false,
            const Handle<YieldTermStructure>& discountCurve =
                Handle<YieldTermStructure>(),
            int innerPoints = 32,
            Real innerStddevs = 7.0);

        Gaussian1dCallableRangeAccrualEngine(
            const Handle<Gaussian1dModel>& model,
            int integrationPoints = 64,
            Real stddevs = 7.0,
            bool extrapolatePayoff = true,
            bool flatPayoffExtrapolation = false,
            const Handle<YieldTermStructure>& discountCurve =
                Handle<YieldTermStructure>(),
            int innerPoints = 32,
            Real innerStddevs = 7.0);

        void calculate() const override;

      private:

        //! Pre-computed root info for one observation date.
        //! y_L and y_U are in standardized state space (y-coordinate).
        //! alwaysIn/alwaysOut flags bypass root finding when the entire
        //! reachable rate range is inside or outside the trigger band.
        struct ObsRoot {
            Time obsTime = 0.0;
            Real y_L = 0.0;    //!< swapRate(obsDate, y_L) = lower
            Real y_U = 0.0;    //!< swapRate(obsDate, y_U) = upper
            bool alwaysIn = false;   //!< rate always in [lower, upper]
            bool alwaysOut = false;  //!< rate always outside
            bool fixedObs = false;   //!< obs in the past, deterministic
            Real fixedRate = 0.0;    //!< rate at fixing (if fixedObs)
        };

        //! Pre-compute roots for all observation dates in a coupon period.
        std::vector<ObsRoot> precomputeRoots(
            Size couponIndex,
            const Date& settlement) const;

        //! Compute P(lower <= swapRate <= upper | y at referenceDate)
        //! using pre-computed roots (CDF-only, no swap rate evaluation).
        Real rangeProbabilityFromRoot(
            const ObsRoot& root,
            Time refTime,
            Real y) const;

        //! Compute the expected range accrual fraction for RA coupon i,
        //! conditional on state y at referenceDate, using pre-computed roots.
        Real expectedAccrualFraction(
            const std::vector<ObsRoot>& roots,
            Time refTime,
            Real y) const;

        //! Compute the NPV of remaining swap cash flows given state y
        //! at the exercise date.
        Real underlyingNpv(
            const Date& expiry,
            Real y,
            const std::vector<std::vector<ObsRoot>>& allRoots) const;

        const int integrationPoints_;
        const Real stddevs_;
        const bool extrapolatePayoff_;
        const bool flatPayoffExtrapolation_;
        const Handle<YieldTermStructure> discountCurve_;
        const int innerPoints_;
        const Real innerStddevs_;
    };
}

#endif
