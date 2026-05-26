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
*/

/*! \file eventvollocalvoltermstructure.hpp
    \brief Event-aware local-vol surface (atomic variance jumps at known dates)

    Wraps a base ``LocalVolTermStructure`` and injects discrete
    variance jumps ``Δ_k`` at known event times ``τ_k``.  Used when
    the ATM term structure has been decomposed (e.g., via chloride's
    ``event_vol`` QP) into a smooth base + per-event jumps, and the
    FD pricer needs to know about the jumps explicitly.

    Implementation: LV-spike smearing.  In a narrow window
    ``[τ_k - δ, τ_k + δ]`` of width ``2δ`` around each event date,
    the local variance is augmented:

        σ_LV²(t, K) = σ_base²(t, K) + Δ_k / (2δ)

    Integrating over the impulse window:
    ``∫_{τ_k − δ}^{τ_k + δ} (σ_base² + Δ_k / 2δ) dt
       = σ_base²(τ_k, K) · 2δ + Δ_k``

    so the event contributes exactly ``Δ_k`` of total variance plus
    the continuous base variance over the window.

    Caveats:

    - This is a *smeared* impulse — the jump happens over ``2δ``
      rather than instantly at ``τ_k``.  Accurate to within FD
      time-step resolution for vanilla payoffs.  For path-dependent
      products (barriers, Asians) where the exact instant matters,
      a true ``FdmStepCondition``-based convolution is required.
      See docs/tasks/abstract_mixed_vol_surface.md §C for the true
      impulse path.
    - Default ``δ = 0.5 / 365`` (half a day).  The FD time grid
      must resolve at least one step inside ``[τ_k − δ, τ_k + δ]``
      or the spike is silently averaged out; use QL's
      ``ConcentratingCurve`` time grid feature when configuring the
      FD engine to concentrate steps around event dates.
*/

#ifndef quantlib_event_vol_local_vol_term_structure_hpp
#define quantlib_event_vol_local_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/handle.hpp>
#include <vector>

namespace QuantLib {

    //! Event-aware local-vol surface with discrete variance jumps
    /*! See file comment for the math and caveats. */
    class EventVolLocalVolTermStructure : public LocalVolTermStructure {
      public:
        EventVolLocalVolTermStructure(
            Handle<LocalVolTermStructure> baseLV,
            std::vector<Time> eventTimes,
            std::vector<Real> eventVariances,
            Real impulseHalfWidth = 0.5 / 365.0);

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
        const std::vector<Time>& eventTimes() const { return eventTimes_; }
        const std::vector<Real>& eventVariances() const { return eventVariances_; }
        Real impulseHalfWidth() const { return impulseHalfWidth_; }
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility localVolImpl(Time t, Real strike) const override;

      private:
        void validate_() const;

        Handle<LocalVolTermStructure> baseLV_;
        std::vector<Time> eventTimes_;
        std::vector<Real> eventVariances_;
        Real impulseHalfWidth_;
    };

}

#endif
