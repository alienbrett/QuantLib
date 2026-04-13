/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Chloride Project

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

/*! \file dualwingessvi.hpp
    \brief Dual-wing eSSVI volatility term structure

    Extension of eSSVI with separate curvature parameters (psi) for
    put and call wings.  C^0 at ATMF (w(0)=theta), kink when psi_lo != psi_hi.
    Each wing is independently arb-free via the Mingone global parametrization.
*/

#ifndef quantlib_dual_wing_essvi_vol_term_structure_hpp
#define quantlib_dual_wing_essvi_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/essvihelpers.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

namespace QuantLib {

    //! Dual-wing eSSVI volatility term structure
    class DualWingEssviVolatilityTermStructure : public BlackVolatilityTermStructure {
      public:
        //! Construct from native per-slice dual-wing parameters
        DualWingEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos,
            const std::vector<Real>& psis_lo,
            const std::vector<Real>& psis_hi,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc = Actual365Fixed());

        //! Construct from global arb-free dual-wing parameters
        DualWingEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos,
            Real theta1,
            const std::vector<Real>& as,
            const std::vector<Real>& cs_lo,
            const std::vector<Real>& cs_hi,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            EssviButterflyCondition::Type bflyType
                = EssviButterflyCondition::GatheralJacquier,
            const DayCounter& dc = Actual365Fixed());

        //! \name TermStructure interface
        //@{
        Date maxDate() const override;
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override;
        Real maxStrike() const override;
        //@}

        //! \name Inspectors
        //@{
        Size numSlices() const;
        const std::vector<DualWingEssviSliceParams>& slices() const;
        const DualWingEssviSurface& surface() const { return surface_; }
        //@}

        //! \name Gradient
        //@{
        DualWingEssviSliceGradient impliedVolGradient(Size sliceIdx, Real strike) const;

        std::vector<Real> impliedVolGlobalGradient(
            Size sliceIdx, Real strike,
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;
        //@}

        //! \name Batch evaluation
        //@{
        std::vector<Real> batchBlackVol(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes) const;

        std::vector<Real> batchImpliedVolGlobalGradient(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes,
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;
        //@}

        //! \name Chain Jacobian
        //@{
        std::vector<Real> chainJacobian(
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const {
            return surface_.chainJacobian(gp, bflyCond);
        }
        //@}

        void accept(AcyclicVisitor&) override;

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Real forward(Time t) const;

        DualWingEssviSurface surface_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
    };

} // namespace QuantLib

#endif
