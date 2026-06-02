/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 chloride contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.
*/

/*! \file scaledvoltermstructure.hpp
    \brief Multiplicative scaling wrappers around Black/local vol term structures.

    Three adapters that scale a base volatility surface:

    1. ScaledBlackVolTermStructure / ScaledLocalVolTermStructure:
       Flat multiplicative scaling by a Quote-driven multiplier alpha:
         sigma_scaled(t, K)   = alpha * sigma_base(t, K)
         sigmaLV_scaled(t, K) = alpha * sigmaLV_base(t, K)

    2. SqrtTimeScaledBlackVolTermStructure:
       Time-decaying multiplicative scaling:
         sigma_scaled(t, K) = m^sqrt(T_ref / t) * sigma_base(t, K)
       where m is a Quote-driven multiplier and T_ref is a reference tenor.
       At t = T_ref the multiplier is m (full shift); short-dated vol is
       amplified, long-dated vol is dampened. This captures the empirical
       observation that instantaneous vol shocks propagate to the term
       structure with 1/sqrt(T) decay.

    For Black vol, scaling implied vol by alpha is equivalent to scaling
    total variance by alpha^2. Under Dupire, local variance also scales
    by alpha^2 (numerator scales; denominator dw/dk terms are invariant).

    Intended use: scenario state generation where many vol-perturbed
    surfaces share a single calibrated base, and the perturbation factor
    is a small scalar driven by a state-vector coordinate.

    All wrappers are observable: changes to the base term structure or
    the Quote(s) propagate to dependent pricers automatically.
*/

#ifndef quantlib_scaled_vol_term_structure_hpp
#define quantlib_scaled_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <cmath>
#include <utility>

namespace QuantLib {

    //! Multiplicative scaling of a Black volatility term structure.
    /*! Wraps a base BlackVolTermStructure handle and scales every Black
        vol query by the value of a Quote. Scaling implied vol by alpha
        is mathematically equivalent to scaling total variance by alpha^2.

        Reference date, day counter, calendar, day-count, and strike range
        are inherited from the base term structure.
    */
    class ScaledBlackVolTermStructure : public BlackVolatilityTermStructure {
      public:
        ScaledBlackVolTermStructure(Handle<BlackVolTermStructure> base,
                                    Handle<Quote> multiplier);

        //! \name TermStructure interface
        //@{
        DayCounter dayCounter() const override { return base_->dayCounter(); }
        Date maxDate() const override { return base_->maxDate(); }
        const Date& referenceDate() const override {
            return base_->referenceDate();
        }
        Calendar calendar() const override { return base_->calendar(); }
        Natural settlementDays() const override {
            return base_->settlementDays();
        }
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override { return base_->minStrike(); }
        Real maxStrike() const override { return base_->maxStrike(); }
        //@}

        //! \name Scaled-surface accessors
        //@{
        const Handle<BlackVolTermStructure>& base() const { return base_; }
        const Handle<Quote>& multiplier() const { return multiplier_; }
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Handle<BlackVolTermStructure> base_;
        Handle<Quote> multiplier_;
    };


    //! Multiplicative scaling of a local volatility term structure.
    /*! Wraps a base LocalVolTermStructure handle and scales every local vol
        query by the value of a Quote. Use the same multiplier as the
        accompanying ScaledBlackVolTermStructure to keep implied vol and
        local vol scaled consistently.
    */
    class ScaledLocalVolTermStructure : public LocalVolTermStructure {
      public:
        ScaledLocalVolTermStructure(Handle<LocalVolTermStructure> base,
                                    Handle<Quote> multiplier);

        //! \name TermStructure interface
        //@{
        DayCounter dayCounter() const override { return base_->dayCounter(); }
        Date maxDate() const override { return base_->maxDate(); }
        const Date& referenceDate() const override {
            return base_->referenceDate();
        }
        Calendar calendar() const override { return base_->calendar(); }
        Natural settlementDays() const override {
            return base_->settlementDays();
        }
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override { return base_->minStrike(); }
        Real maxStrike() const override { return base_->maxStrike(); }
        //@}

        //! \name Scaled-surface accessors
        //@{
        const Handle<LocalVolTermStructure>& base() const { return base_; }
        const Handle<Quote>& multiplier() const { return multiplier_; }
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility localVolImpl(Time t, Real underlyingLevel) const override;

      private:
        Handle<LocalVolTermStructure> base_;
        Handle<Quote> multiplier_;
    };


    //! Time-decaying multiplicative scaling of a Black vol term structure.
    /*! Wraps a base BlackVolTermStructure handle and applies a
        time-dependent multiplier:

            sigma_scaled(t, K) = m^sqrt(T_ref / t) * sigma_base(t, K)

        where m is a Quote-driven multiplier (e.g. 1.10 for +10% vol at
        the reference tenor) and T_ref is a fixed reference tenor in year
        fractions. The exponent sqrt(T_ref / t) makes the scaling
        strongest at short tenors and decay toward long tenors, consistent
        with how instantaneous variance shocks propagate through the term
        structure.

        At t = T_ref: multiplier = m (full shift).
        At t = 4*T_ref: multiplier = sqrt(m) (half the log-shift).
        m = 1 is a no-op.

        For t < tFloor (default 1/365), t is clamped to tFloor to prevent
        blow-up at very short tenors.

        Reference date, day counter, calendar, and strike range are
        inherited from the base term structure.
    */
    class SqrtTimeScaledBlackVolTermStructure
        : public BlackVolatilityTermStructure {
      public:
        SqrtTimeScaledBlackVolTermStructure(
            Handle<BlackVolTermStructure> base,
            Handle<Quote> multiplier,
            Time tRef,
            Time tFloor = 1.0 / 365.0);

        //! \name TermStructure interface
        //@{
        DayCounter dayCounter() const override { return base_->dayCounter(); }
        Date maxDate() const override { return base_->maxDate(); }
        const Date& referenceDate() const override {
            return base_->referenceDate();
        }
        Calendar calendar() const override { return base_->calendar(); }
        Natural settlementDays() const override {
            return base_->settlementDays();
        }
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override { return base_->minStrike(); }
        Real maxStrike() const override { return base_->maxStrike(); }
        //@}

        //! \name Accessors
        //@{
        const Handle<BlackVolTermStructure>& base() const { return base_; }
        const Handle<Quote>& multiplier() const { return multiplier_; }
        Time tRef() const { return tRef_; }
        Time tFloor() const { return tFloor_; }
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Handle<BlackVolTermStructure> base_;
        Handle<Quote> multiplier_;
        Time tRef_;
        Time tFloor_;
    };


    // ── inline definitions ────────────────────────────────────────────────

    inline ScaledBlackVolTermStructure::ScaledBlackVolTermStructure(
        Handle<BlackVolTermStructure> base,
        Handle<Quote> multiplier)
    : BlackVolatilityTermStructure(base->businessDayConvention(),
                                   base->dayCounter()),
      base_(std::move(base)),
      multiplier_(std::move(multiplier)) {
        registerWith(base_);
        registerWith(multiplier_);
    }

    inline Volatility
    ScaledBlackVolTermStructure::blackVolImpl(Time t, Real strike) const {
        Real m = multiplier_->value();
        // Extrapolation is safe: the wrapper has already passed range checks.
        Volatility v = base_->blackVol(t, strike, true);
        return m * v;
    }

    inline void ScaledBlackVolTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<ScaledBlackVolTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }


    inline ScaledLocalVolTermStructure::ScaledLocalVolTermStructure(
        Handle<LocalVolTermStructure> base,
        Handle<Quote> multiplier)
    : LocalVolTermStructure(base->businessDayConvention(),
                            base->dayCounter()),
      base_(std::move(base)),
      multiplier_(std::move(multiplier)) {
        registerWith(base_);
        registerWith(multiplier_);
    }

    inline Volatility
    ScaledLocalVolTermStructure::localVolImpl(Time t,
                                              Real underlyingLevel) const {
        Real m = multiplier_->value();
        Volatility v = base_->localVol(t, underlyingLevel, true);
        return m * v;
    }

    inline void ScaledLocalVolTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<ScaledLocalVolTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            LocalVolTermStructure::accept(v);
    }


    inline SqrtTimeScaledBlackVolTermStructure::
        SqrtTimeScaledBlackVolTermStructure(
            Handle<BlackVolTermStructure> base,
            Handle<Quote> multiplier,
            Time tRef,
            Time tFloor)
    : BlackVolatilityTermStructure(base->businessDayConvention(),
                                   base->dayCounter()),
      base_(std::move(base)),
      multiplier_(std::move(multiplier)),
      tRef_(tRef),
      tFloor_(tFloor) {
        QL_REQUIRE(tRef_ > 0.0,
                   "SqrtTimeScaledBlackVolTermStructure: tRef must be > 0, "
                   "got " << tRef_);
        QL_REQUIRE(tFloor_ > 0.0,
                   "SqrtTimeScaledBlackVolTermStructure: tFloor must be > 0, "
                   "got " << tFloor_);
        registerWith(base_);
        registerWith(multiplier_);
    }

    inline Volatility
    SqrtTimeScaledBlackVolTermStructure::blackVolImpl(Time t,
                                                      Real strike) const {
        Real m = multiplier_->value();
        Volatility v = base_->blackVol(t, strike, true);
        if (m == 1.0)
            return v;
        Time tEff = std::max(t, tFloor_);
        // m^sqrt(T_ref / t) = exp(ln(m) * sqrt(T_ref / t))
        Real alpha = std::exp(std::log(m) * std::sqrt(tRef_ / tEff));
        return alpha * v;
    }

    inline void
    SqrtTimeScaledBlackVolTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 =
            dynamic_cast<Visitor<SqrtTimeScaledBlackVolTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

}

#endif
