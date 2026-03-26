/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Chloride Project
*/

#ifndef quantlib_split_rho_essvi_vol_term_structure_hpp
#define quantlib_split_rho_essvi_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/essvihelpers.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

namespace QuantLib {

    //! Split-rho eSSVI volatility term structure (5 params per slice)
    class SplitRhoEssviVolatilityTermStructure : public BlackVolatilityTermStructure {
      public:
        //! Construct from native per-slice parameters
        SplitRhoEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos_lo,
            const std::vector<Real>& psis_lo,
            const std::vector<Real>& rhos_hi,
            const std::vector<Real>& psis_hi,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc = Actual365Fixed());

        //! Construct from global arb-free parameters
        SplitRhoEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos_lo,
            const std::vector<Real>& rhos_hi,
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

        Date maxDate() const override;
        Real minStrike() const override;
        Real maxStrike() const override;

        Size numSlices() const;
        const std::vector<SplitRhoEssviSliceParams>& slices() const;
        const SplitRhoEssviSurface& surface() const { return surface_; }

        SplitRhoEssviSliceGradient impliedVolGradient(Size sliceIdx, Real strike) const;
        std::vector<Real> impliedVolGlobalGradient(
            Size sliceIdx, Real strike,
            const SplitRhoEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;

        void accept(AcyclicVisitor&) override;

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Real forward(Time t) const;

        SplitRhoEssviSurface surface_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
    };

} // namespace QuantLib

#endif
