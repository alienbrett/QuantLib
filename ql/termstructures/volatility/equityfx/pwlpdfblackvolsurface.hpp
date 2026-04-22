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

/*! \file pwlpdfblackvolsurface.hpp
    \brief Forward-aware black vol surface with per-slice log-moneyness grids
*/

#ifndef quantlib_pwl_pdf_black_vol_surface_hpp
#define quantlib_pwl_pdf_black_vol_surface_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/math/interpolation.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/instruments/dividendschedule.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

namespace QuantLib {

    //! Forward-aware black vol surface for tabulated per-slice smiles
    /*! Stores per-expiry implied vol smiles on log-moneyness grids
        (k = ln(K/F)).  All interpolation happens in moneyness space,
        avoiding the negative-variance artifacts that BlackVarianceSurface
        produces when bilinearly interpolating in absolute strike space.

        Forward-aware: takes spot + rate handles, computes F(T) = S·dq/dr
        internally — same pattern as EssviVolatilityTermStructure.

        Between expiries, total variance w(k,T) = σ²T is linearly
        interpolated in T at matching log-moneyness.

        Per-slice interpolation is linear by default.  A template
        setInterpolation<Interpolator>() method is declared for future
        cubic/monotone variants but only linear is exposed for now.

        \ingroup termstructures
    */
    class PwlPdfBlackVolSurface : public BlackVolatilityTermStructure {
      public:
        //! Continuous dividend yield only
        PwlPdfBlackVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& kGrids,
            const std::vector<std::vector<Real>>& impliedVols,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc = Actual365Fixed(),
            const std::vector<Real>& calibrationForwards
                = std::vector<Real>());

        //! With discrete dividends (same signature as eSSVI)
        PwlPdfBlackVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& kGrids,
            const std::vector<std::vector<Real>>& impliedVols,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
            const DayCounter& dc = Actual365Fixed(),
            const std::vector<Real>& calibrationForwards
                = std::vector<Real>());

        //! \name TermStructure interface
        //@{
        Date maxDate() const override;
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override;
        Real maxStrike() const override;
        //@}

        //! \name Modifiers
        //@{
        /*! Set per-slice interpolation type.  Only Linear is currently
            exposed via SWIG; the template is here for future extensions
            (e.g. Cubic, MonotonicCubic). */
        template <class Interpolator>
        void setInterpolation(const Interpolator& i = Interpolator());
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

        //! \name Internal types (public for static init helper)
        //@{
        struct Slice {
            Time T;
            Real calForward;              // calibration forward (0 = use handles)
            std::vector<Real> k;          // log-moneyness grid
            std::vector<Real> totalVar;   // σ²T at each grid node
            mutable Interpolation interp; // lazy-built interpolation
        };
        //@}

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Real forward(Time t) const;
        Real forwardAtSlice(Size idx) const;
        Real interpolateSlice(Size idx, Real k) const;

        std::vector<Slice> slices_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        DividendSchedule dividends_;
    };


    // =====================================================================
    // inline / template definitions
    // =====================================================================

    template <class Interpolator>
    void PwlPdfBlackVolSurface::setInterpolation(const Interpolator& i) {
        for (auto& sl : slices_) {
            sl.interp = i.interpolate(sl.k.begin(), sl.k.end(),
                                      sl.totalVar.begin());
            sl.interp.enableExtrapolation();
        }
        notifyObservers();
    }

    inline void PwlPdfBlackVolSurface::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<PwlPdfBlackVolSurface>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

} // namespace QuantLib

#endif
