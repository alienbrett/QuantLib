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

/*! \file essvilocalvolsurface.hpp
    \brief Analytic Dupire local volatility from eSSVI surface
*/

#ifndef quantlib_essvi_local_vol_surface_hpp
#define quantlib_essvi_local_vol_surface_hpp

#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/essvivoltermstructure.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>

namespace QuantLib {

    //! Analytic Dupire local vol from an eSSVI implied vol surface
    /*! Computes local volatility via Dupire's formula using analytic
        derivatives of the eSSVI total variance:
          - dw/dk and d²w/dk² from closed-form eSSVI strike derivatives
          - dw/dT from differentiation of the linear param interpolation

        This avoids the finite-difference bumping used by LocalVolSurface
        and is both faster and more accurate.
    */
    class EssviLocalVolSurface : public LocalVolTermStructure {
      public:
        EssviLocalVolSurface(
            ext::shared_ptr<EssviVolatilityTermStructure> essviSurface,
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
        ext::shared_ptr<EssviVolatilityTermStructure> essviSurface_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        Handle<Quote> spot_;
    };

} // namespace QuantLib

#endif
