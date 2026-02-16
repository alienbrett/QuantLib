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

/*! \file mynewvolsurface.hpp
    \brief Custom pedagogical volatility surface for testing build pipeline
*/

#ifndef quantlib_mynewvolsurface_hpp
#define quantlib_mynewvolsurface_hpp

#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <utility>

namespace QuantLib {

    //! My New Custom Volatility Surface
    /*! This is a pedagogical example to demonstrate adding custom code to QuantLib.
        It implements a simple constant volatility surface for testing purposes.

        This class demonstrates the full workflow:
        1. Adding C++ code to QuantLib
        2. Building it into the library
        3. Exposing it via SWIG
        4. Using it from Python
    */
    class MyNewVolSurface : public BlackVolatilityTermStructure {
      public:
        //! Constructor with explicit volatility value
        MyNewVolSurface(const Date& referenceDate,
                        const Calendar& calendar,
                        Volatility volatility,
                        const DayCounter& dayCounter);

        //! Constructor with Quote handle (allows dynamic updates)
        MyNewVolSurface(const Date& referenceDate,
                        const Calendar& calendar,
                        Handle<Quote> volatility,
                        const DayCounter& dayCounter);

        //! Constructor with settlement days and explicit volatility
        MyNewVolSurface(Natural settlementDays,
                        const Calendar& calendar,
                        Volatility volatility,
                        const DayCounter& dayCounter);

        //! Constructor with settlement days and Quote handle
        MyNewVolSurface(Natural settlementDays,
                        const Calendar& calendar,
                        Handle<Quote> volatility,
                        const DayCounter& dayCounter);

        //! \name TermStructure interface
        //@{
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

        //! Custom method: get the base volatility value
        /*! This demonstrates adding custom functionality */
        Volatility baseVolatility() const;

      protected:
        //! Implementation of BlackVolatilityTermStructure interface
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Handle<Quote> volatility_;
    };


    // Inline implementations

    inline MyNewVolSurface::MyNewVolSurface(const Date& referenceDate,
                                            const Calendar& cal,
                                            Volatility volatility,
                                            const DayCounter& dc)
    : BlackVolatilityTermStructure(referenceDate, cal, Following, dc),
      volatility_(ext::shared_ptr<Quote>(new SimpleQuote(volatility))) {}

    inline MyNewVolSurface::MyNewVolSurface(const Date& referenceDate,
                                            const Calendar& cal,
                                            Handle<Quote> volatility,
                                            const DayCounter& dc)
    : BlackVolatilityTermStructure(referenceDate, cal, Following, dc),
      volatility_(std::move(volatility)) {
        registerWith(volatility_);
    }

    inline MyNewVolSurface::MyNewVolSurface(Natural settlementDays,
                                            const Calendar& cal,
                                            Volatility volatility,
                                            const DayCounter& dc)
    : BlackVolatilityTermStructure(settlementDays, cal, Following, dc),
      volatility_(ext::shared_ptr<Quote>(new SimpleQuote(volatility))) {}

    inline MyNewVolSurface::MyNewVolSurface(Natural settlementDays,
                                            const Calendar& cal,
                                            Handle<Quote> volatility,
                                            const DayCounter& dc)
    : BlackVolatilityTermStructure(settlementDays, cal, Following, dc),
      volatility_(std::move(volatility)) {
        registerWith(volatility_);
    }

    inline Date MyNewVolSurface::maxDate() const {
        return Date::maxDate();
    }

    inline Real MyNewVolSurface::minStrike() const {
        return QL_MIN_REAL;
    }

    inline Real MyNewVolSurface::maxStrike() const {
        return QL_MAX_REAL;
    }

    inline void MyNewVolSurface::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<MyNewVolSurface>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

    inline Volatility MyNewVolSurface::blackVolImpl(Time, Real) const {
        // Returns constant volatility regardless of time and strike
        return volatility_->value();
    }

    inline Volatility MyNewVolSurface::baseVolatility() const {
        return volatility_->value();
    }

}

#endif