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

/*! \file essvivoltermstructure.hpp
    \brief eSSVI volatility term structure

    Global arbitrage-free eSSVI volatility surface as a QuantLib
    BlackVolatilityTermStructure. Based on:
      Mingone (2022), "No arbitrage global parametrization for the eSSVI
      volatility surface", arXiv:2204.00312v1.

    Supports discrete dividends: when a DividendSchedule is provided,
    the forward is computed as:
      F(T) = (S - PV(fixed divs before T)) * exp(-q*T) / exp(-r*T)
    where PV discounts each dividend at the risk-free rate.
*/

#ifndef quantlib_essvi_vol_term_structure_hpp
#define quantlib_essvi_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/essvihelpers.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/instruments/dividendschedule.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

namespace QuantLib {

    //! eSSVI volatility term structure
    /*! Implements the extended Surface SVI (eSSVI) volatility surface
        with the global arbitrage-free parametrization from Mingone (2022).

        The surface can be constructed from either:
        - Native per-slice parameters (theta, rho, psi)
        - Global arb-free parameters (rhos, theta1, as, cs)

        The surface requires a spot price and yield curves to convert
        between strikes and log-forward moneyness internally.

        When a DividendSchedule is provided, discrete cash dividends
        are subtracted from the spot (PV'd at the risk-free rate)
        before computing the forward. This ensures log-moneyness is
        computed against the correct forward price.
    */
    class EssviVolatilityTermStructure : public BlackVolatilityTermStructure {
      public:
        //! Construct from native per-slice SSVI parameters
        EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos,
            const std::vector<Real>& psis,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc = Actual365Fixed());

        //! Construct from native params with discrete dividends
        EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos,
            const std::vector<Real>& psis,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
            const DayCounter& dc = Actual365Fixed());

        //! Construct from global arb-free parameters (Proposition 3.1)
        EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos,
            Real theta1,
            const std::vector<Real>& as,
            const std::vector<Real>& cs,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            EssviButterflyCondition::Type bflyType
                = EssviButterflyCondition::GatheralJacquier,
            const DayCounter& dc = Actual365Fixed());

        //! Construct from global arb-free parameters with discrete dividends
        EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos,
            Real theta1,
            const std::vector<Real>& as,
            const std::vector<Real>& cs,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
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
        const std::vector<EssviSliceParams>& slices() const;
        const EssviSurface& surface() const { return surface_; }
        //@}

        //! \name Gradient
        //@{
        //! Gradient of blackVol w.r.t. native (theta_i, rho_i, psi_i) for slice i
        EssviSliceGradient impliedVolGradient(Size sliceIdx, Real strike) const;

        //! Gradient of blackVol w.r.t. global params [rhos, theta1, as, cs]
        std::vector<Real> impliedVolGlobalGradient(
            Size sliceIdx, Real strike,
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;
        //@}

        //! \name Batch evaluation (vectorized, avoids per-call SWIG overhead)
        //@{
        /*! Evaluate blackVol for multiple (sliceIdx, strike) pairs.
            Returns vector of length N = sliceIndices.size().
        */
        std::vector<Real> batchBlackVol(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes) const;

        /*! Evaluate impliedVolGlobalGradient for multiple observations.
            Returns flat row-major matrix: N_obs rows × N_params cols.
        */
        std::vector<Real> batchImpliedVolGlobalGradient(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes,
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Real forward(Time t) const;

        EssviSurface surface_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        DividendSchedule dividends_;
    };

} // namespace QuantLib

#endif
