/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 chloride

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

/*! \file rangeaccrualswap.hpp
    \brief swap with a range accrual leg and a fixed funding leg,
           supporting bespoke per-period parameters
*/

#ifndef quantlib_range_accrual_swap_hpp
#define quantlib_range_accrual_swap_hpp

#include <ql/instruments/swap.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/schedule.hpp>

namespace QuantLib {

    class SwapIndex;

    //! Range accrual swap
    /*! Fixed leg vs range-accrual leg.  Each RA coupon accrues at
        (gearing * index_rate + spread) * accrual_fraction, where
        accrual_fraction is the proportion of observation dates for
        which the swap rate (or other index) falls within [lower, upper].

        All vectors are per-period and must have size = schedule.size()-1.

        \ingroup instruments
    */

    class RangeAccrualSwap : public Swap {
      public:
        class arguments;
        class results;
        class engine;

        //! Receiver receives the RA leg, Payer pays the RA leg
        RangeAccrualSwap(
            Swap::Type type,
            // fixed (funding) leg
            std::vector<Real> fixedNominal,
            Schedule fixedSchedule,
            std::vector<Real> fixedRate,
            DayCounter fixedDayCount,
            // range accrual leg
            std::vector<Real> raNominal,
            Schedule raSchedule,
            ext::shared_ptr<SwapIndex> observationIndex,
            Period observationTenor,
            std::vector<Real> raGearings,
            std::vector<Spread> raSpreads,
            std::vector<Rate> lowerTriggers,
            std::vector<Rate> upperTriggers,
            DayCounter raDayCount,
            BusinessDayConvention paymentConvention = Following,
            BusinessDayConvention observationConvention = ModifiedFollowing);

        //! \name Inspectors
        //@{
        Swap::Type type() const;
        const std::vector<Real>& fixedNominal() const;
        const Schedule& fixedSchedule() const;
        const std::vector<Real>& fixedRate() const;
        const DayCounter& fixedDayCount() const;

        const std::vector<Real>& raNominal() const;
        const Schedule& raSchedule() const;
        const ext::shared_ptr<SwapIndex>& observationIndex() const;
        const Period& observationTenor() const;
        const std::vector<Real>& raGearings() const;
        const std::vector<Spread>& raSpreads() const;
        const std::vector<Rate>& lowerTriggers() const;
        const std::vector<Rate>& upperTriggers() const;
        const DayCounter& raDayCount() const;
        BusinessDayConvention paymentConvention() const;
        BusinessDayConvention observationConvention() const;

        const Leg& fixedLeg() const;
        const Leg& raLeg() const;

        //! Observation dates for the i-th RA coupon period
        const std::vector<Date>& observationDates(Size i) const;
        //@}

        void setupArguments(PricingEngine::arguments* args) const override;
        void fetchResults(const PricingEngine::results*) const override;

      private:
        void init();
        void setupExpired() const override;

        Swap::Type type_;
        // fixed leg
        std::vector<Real> fixedNominal_;
        Schedule fixedSchedule_;
        std::vector<Real> fixedRate_;
        DayCounter fixedDayCount_;
        // RA leg
        std::vector<Real> raNominal_;
        Schedule raSchedule_;
        ext::shared_ptr<SwapIndex> observationIndex_;
        Period observationTenor_;
        std::vector<Real> raGearings_;
        std::vector<Spread> raSpreads_;
        std::vector<Rate> lowerTriggers_;
        std::vector<Rate> upperTriggers_;
        DayCounter raDayCount_;
        BusinessDayConvention paymentConvention_;
        BusinessDayConvention observationConvention_;

        // cached per-period observation dates
        std::vector<std::vector<Date>> observationDates_;
    };

    //! %Arguments for range accrual swap calculation
    class RangeAccrualSwap::arguments : public Swap::arguments {
      public:
        arguments() = default;
        Swap::Type type = Swap::Receiver;

        // fixed leg
        std::vector<Real> fixedNominal;
        std::vector<Date> fixedResetDates;
        std::vector<Date> fixedPayDates;
        std::vector<Real> fixedCoupons;
        std::vector<Real> fixedRate;

        // range accrual leg
        std::vector<Real> raNominal;
        std::vector<Date> raResetDates;
        std::vector<Date> raPayDates;
        std::vector<Time> raAccrualTimes;
        std::vector<Real> raGearings;
        std::vector<Spread> raSpreads;
        std::vector<Rate> lowerTriggers;
        std::vector<Rate> upperTriggers;

        // per-period observation dates (outer: coupon index, inner: obs dates)
        std::vector<std::vector<Date>> observationDates;

        ext::shared_ptr<SwapIndex> observationIndex;
        Period observationTenor;

        void validate() const override;
    };

    //! %Results from range accrual swap calculation
    class RangeAccrualSwap::results : public Swap::results {
      public:
        void reset() override;
    };

    class RangeAccrualSwap::engine
        : public GenericEngine<RangeAccrualSwap::arguments,
                               RangeAccrualSwap::results> {};

    // inline definitions

    inline Swap::Type RangeAccrualSwap::type() const { return type_; }

    inline const std::vector<Real>& RangeAccrualSwap::fixedNominal() const {
        return fixedNominal_;
    }

    inline const Schedule& RangeAccrualSwap::fixedSchedule() const {
        return fixedSchedule_;
    }

    inline const std::vector<Real>& RangeAccrualSwap::fixedRate() const {
        return fixedRate_;
    }

    inline const DayCounter& RangeAccrualSwap::fixedDayCount() const {
        return fixedDayCount_;
    }

    inline const std::vector<Real>& RangeAccrualSwap::raNominal() const {
        return raNominal_;
    }

    inline const Schedule& RangeAccrualSwap::raSchedule() const {
        return raSchedule_;
    }

    inline const ext::shared_ptr<SwapIndex>&
    RangeAccrualSwap::observationIndex() const {
        return observationIndex_;
    }

    inline const Period& RangeAccrualSwap::observationTenor() const {
        return observationTenor_;
    }

    inline const std::vector<Real>& RangeAccrualSwap::raGearings() const {
        return raGearings_;
    }

    inline const std::vector<Spread>& RangeAccrualSwap::raSpreads() const {
        return raSpreads_;
    }

    inline const std::vector<Rate>& RangeAccrualSwap::lowerTriggers() const {
        return lowerTriggers_;
    }

    inline const std::vector<Rate>& RangeAccrualSwap::upperTriggers() const {
        return upperTriggers_;
    }

    inline const DayCounter& RangeAccrualSwap::raDayCount() const {
        return raDayCount_;
    }

    inline BusinessDayConvention RangeAccrualSwap::paymentConvention() const {
        return paymentConvention_;
    }

    inline BusinessDayConvention
    RangeAccrualSwap::observationConvention() const {
        return observationConvention_;
    }

    inline const Leg& RangeAccrualSwap::fixedLeg() const { return legs_[0]; }

    inline const Leg& RangeAccrualSwap::raLeg() const { return legs_[1]; }

    inline const std::vector<Date>&
    RangeAccrualSwap::observationDates(Size i) const {
        return observationDates_[i];
    }
}

#endif
