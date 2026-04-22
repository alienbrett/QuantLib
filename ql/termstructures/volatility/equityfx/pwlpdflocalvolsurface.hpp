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

/*! \file pwlpdflocalvolsurface.hpp
    \brief Analytical Dupire local vol from piecewise-linear PDF surface

    Uses Dupire's formula with the PDF available analytically:
      sigma_loc^2 = (dC/dT) / (0.5 * K^2 * p(K,T))
    where p(K,T) is the strike-space density = q(k) / (F * K),
    and dC/dT is a finite difference of analytically-computed call prices
    between adjacent expiry slices.

    Eliminates the double finite-difference noise of NoExceptLocalVolSurface
    (which bumps the black vol surface to estimate d^2C/dK^2 and dC/dT).
    Here the denominator (butterfly density) is exact from the QP.
*/

#ifndef quantlib_pwl_pdf_local_vol_surface_hpp
#define quantlib_pwl_pdf_local_vol_surface_hpp

#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/pwlpdfsurface.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>

namespace QuantLib {

    //! Analytical Dupire local vol from a PwlPdfVolSurface
    class PwlPdfLocalVolSurface : public LocalVolTermStructure {
      public:
        PwlPdfLocalVolSurface(
            ext::shared_ptr<PwlPdfVolSurface> pdfSurface,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Handle<Quote> spot);

        //! \name TermStructure interface
        //@{
        const Date& referenceDate() const override;
        DayCounter dayCounter() const override;
        Date maxDate() const override;
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override;
        Real maxStrike() const override;
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility localVolImpl(Time t, Real underlyingLevel) const override;

      private:
        ext::shared_ptr<PwlPdfVolSurface> pdfSurface_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        Handle<Quote> spot_;
    };

} // namespace QuantLib

#endif
