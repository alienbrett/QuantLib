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

#include <ql/termstructures/volatility/equityfx/eventvollocalvoltermstructure.hpp>
#include <ql/errors.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    namespace {
        constexpr Real kTimeFloor = 1e-14;
        constexpr Real kVolFloor  = 1e-6;
        constexpr Real kVolCap    = 100.0;   // 10000% LV ceiling (safety on spike)
    }

    EventVolLocalVolTermStructure::EventVolLocalVolTermStructure(
            Handle<LocalVolTermStructure> baseLV,
            std::vector<Time> eventTimes,
            std::vector<Real> eventVariances,
            Real impulseHalfWidth)
        : LocalVolTermStructure(
              baseLV->referenceDate(),
              baseLV->calendar(),
              baseLV->businessDayConvention(),
              baseLV->dayCounter()),
          baseLV_(std::move(baseLV)),
          eventTimes_(std::move(eventTimes)),
          eventVariances_(std::move(eventVariances)),
          impulseHalfWidth_(impulseHalfWidth) {
        validate_();
        registerWith(baseLV_);
    }

    void EventVolLocalVolTermStructure::validate_() const {
        QL_REQUIRE(!baseLV_.empty(), "baseLV handle must not be empty");
        QL_REQUIRE(eventTimes_.size() == eventVariances_.size(),
                   "eventTimes (" << eventTimes_.size()
                   << ") and eventVariances (" << eventVariances_.size()
                   << ") must have equal length");
        QL_REQUIRE(impulseHalfWidth_ > 0.0,
                   "impulseHalfWidth must be > 0, got " << impulseHalfWidth_);
        for (Size i = 0; i < eventTimes_.size(); ++i) {
            QL_REQUIRE(eventTimes_[i] > 0.0,
                       "eventTimes[" << i << "] = " << eventTimes_[i]
                       << " must be > 0");
            QL_REQUIRE(eventVariances_[i] >= 0.0,
                       "eventVariances[" << i << "] = " << eventVariances_[i]
                       << " must be >= 0");
        }
        // Sort events by time (helper for window lookup) — caller responsibility
        // but we don't enforce it here since the search is linear in event count.
    }

    Date EventVolLocalVolTermStructure::maxDate() const {
        return baseLV_->maxDate();
    }

    Real EventVolLocalVolTermStructure::minStrike() const {
        return baseLV_->minStrike();
    }

    Real EventVolLocalVolTermStructure::maxStrike() const {
        return baseLV_->maxStrike();
    }

    void EventVolLocalVolTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 =
            dynamic_cast<Visitor<EventVolLocalVolTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            LocalVolTermStructure::accept(v);
    }

    Volatility EventVolLocalVolTermStructure::localVolImpl(
            Time t, Real strike) const {
        if (t < kTimeFloor) t = kTimeFloor;

        Real base_lv = baseLV_->localVol(t, strike, true);
        Real base_var = base_lv * base_lv;

        // Find any event whose impulse window covers t
        Real spike_var = 0.0;
        for (Size i = 0; i < eventTimes_.size(); ++i) {
            Real dt = std::abs(t - eventTimes_[i]);
            if (dt < impulseHalfWidth_) {
                // Per-event variance rate over the 2δ window:
                //   ν_event = Δ_k / (2δ)
                spike_var += eventVariances_[i] / (2.0 * impulseHalfWidth_);
            }
        }

        Real total_var = base_var + spike_var;
        Real lv = std::sqrt(std::max(total_var, kVolFloor * kVolFloor));
        return std::min(lv, kVolCap);
    }

} // namespace QuantLib
