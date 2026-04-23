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

/*! \file pwlpdfsurface.hpp
    \brief Black vol surface from piecewise-linear PDF in moneyness space

    Takes per-expiry moneyness grids (x = K/F), moneyness-space PDF values
    q_x (the QP fitter's native variable), and calibration forwards.
    Computes call forwards by exact polynomial integration, then inverts
    to Black vol via LiRS solver.

    Paired with PwlPdfLocalVolSurface for analytical Dupire local vol.
*/

#ifndef quantlib_pwl_pdf_surface_hpp
#define quantlib_pwl_pdf_surface_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/instruments/dividendschedule.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

namespace QuantLib {

    //! Black vol surface from piecewise-linear moneyness-space PDF
    /*! Stores per-expiry PWL PDFs on moneyness grids (x = K/F).
        Call forward prices CF/F are computed by exact polynomial
        integration of the PWL density — no transcendental functions.

        Between expiries, total variance (σ²·T) at fixed strike is
        linearly interpolated in time — calendar-arb-free by construction.

        \ingroup termstructures
    */
    class PwlPdfVolSurface : public BlackVolatilityTermStructure {
      public:
        //! Per-expiry PDF slice
        struct PdfSlice {
            Time T;
            Real calForward;
            std::vector<Real> x;    // moneyness grid x = K/F (sorted)
            std::vector<Real> q;    // moneyness-space PDF q_x values
            std::vector<Real> h;    // x[i+1] - x[i] (spacings)
            // Wing extrapolation: linear total variance w(k) = σ²T.
            // Boundary vol + slope computed from the outermost grid
            // points where IV inversion succeeds.  Linear w(k) gives
            // C¹-smooth vol surface → well-defined Dupire everywhere.
            Real leftWingVol = 0.20;
            Real rightWingVol = 0.20;
            Real leftBoundaryK = -1.0;    // log-moneyness at left boundary
            Real rightBoundaryK = 1.0;    // log-moneyness at right boundary
            Real leftWingSlope = 0.0;     // dw/dk at left boundary (≤0)
            Real rightWingSlope = 0.0;    // dw/dk at right boundary (≥0)
        };

        //! Continuous dividend yield only
        PwlPdfVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& xGrids,
            const std::vector<std::vector<Real>>& qValues,
            const std::vector<Real>& calibrationForwards,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc = Actual365Fixed());

        //! With discrete dividends
        PwlPdfVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& xGrids,
            const std::vector<std::vector<Real>>& qValues,
            const std::vector<Real>& calibrationForwards,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
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
        /*! Normalized call forward CF/F at moneyness x and time t. */
        Real callForward(Real x, Time t) const;

        /*! Moneyness-space PDF q_x at moneyness x and time t. */
        Real pdfValue(Real x, Time t) const;

        /*! Survival function S(x) = integral from x to x_max of q(u) du.
            dC/dK = -S(K/F). */
        Real survivalFunction(Real x, Time t) const;

        const std::vector<PdfSlice>& pdfSlices() const { return slices_; }
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;
        //! Override default vol²×t to add epsilon-safe calendar floor.
        //! blackForwardVariance (non-virtual in base) asserts var(t2)>=var(t1).
        //! Our blackVolImpl enforces monotonicity but floating-point noise
        //! can produce sub-epsilon violations → crash in FD engine.
        //! This override adds a tiny upward nudge proportional to t so
        //! var(t+dt) > var(t) always holds numerically.
        Real blackVarianceImpl(Time t, Real strike) const override;

      private:
        Real forward(Time t) const;
        Real forwardAtSlice(Size idx) const;
        Real callForwardSlice(Size idx, Real xEval) const;
        Real interpolatePdf(Size idx, Real xEval) const;
        Real survivalSlice(Size idx, Real xEval) const;

        std::vector<PdfSlice> slices_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        DividendSchedule dividends_;
    };

    inline void PwlPdfVolSurface::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<PwlPdfVolSurface>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

} // namespace QuantLib

#endif
