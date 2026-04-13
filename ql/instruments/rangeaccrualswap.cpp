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

#include <ql/instruments/rangeaccrualswap.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/time/schedule.hpp>
#include <utility>

namespace QuantLib {

    RangeAccrualSwap::RangeAccrualSwap(
        Swap::Type type,
        std::vector<Real> fixedNominal,
        Schedule fixedSchedule,
        std::vector<Real> fixedRate,
        DayCounter fixedDayCount,
        std::vector<Real> raNominal,
        Schedule raSchedule,
        ext::shared_ptr<SwapIndex> observationIndex,
        Period observationTenor,
        std::vector<Real> raGearings,
        std::vector<Spread> raSpreads,
        std::vector<Rate> lowerTriggers,
        std::vector<Rate> upperTriggers,
        DayCounter raDayCount,
        BusinessDayConvention paymentConvention,
        BusinessDayConvention observationConvention)
    : Swap(2), type_(type),
      fixedNominal_(std::move(fixedNominal)),
      fixedSchedule_(std::move(fixedSchedule)),
      fixedRate_(std::move(fixedRate)),
      fixedDayCount_(std::move(fixedDayCount)),
      raNominal_(std::move(raNominal)),
      raSchedule_(std::move(raSchedule)),
      observationIndex_(std::move(observationIndex)),
      observationTenor_(std::move(observationTenor)),
      raGearings_(std::move(raGearings)),
      raSpreads_(std::move(raSpreads)),
      lowerTriggers_(std::move(lowerTriggers)),
      upperTriggers_(std::move(upperTriggers)),
      raDayCount_(std::move(raDayCount)),
      paymentConvention_(paymentConvention),
      observationConvention_(observationConvention) {

        init();
    }

    void RangeAccrualSwap::init() {

        Size nFixed = fixedSchedule_.size() - 1;
        Size nRA = raSchedule_.size() - 1;

        QL_REQUIRE(fixedNominal_.size() == nFixed,
                   "Fixed nominal size (" << fixedNominal_.size()
                   << ") does not match schedule size ("
                   << fixedSchedule_.size() << ") - 1");
        QL_REQUIRE(fixedRate_.size() == nFixed,
                   "Fixed rate size (" << fixedRate_.size()
                   << ") does not match schedule size ("
                   << fixedSchedule_.size() << ") - 1");
        QL_REQUIRE(raNominal_.size() == nRA,
                   "RA nominal size (" << raNominal_.size()
                   << ") does not match RA schedule size ("
                   << raSchedule_.size() << ") - 1");
        QL_REQUIRE(raGearings_.size() == nRA,
                   "RA gearings size (" << raGearings_.size()
                   << ") does not match RA schedule size - 1");
        QL_REQUIRE(raSpreads_.size() == nRA,
                   "RA spreads size (" << raSpreads_.size()
                   << ") does not match RA schedule size - 1");
        QL_REQUIRE(lowerTriggers_.size() == nRA,
                   "Lower triggers size (" << lowerTriggers_.size()
                   << ") does not match RA schedule size - 1");
        QL_REQUIRE(upperTriggers_.size() == nRA,
                   "Upper triggers size (" << upperTriggers_.size()
                   << ") does not match RA schedule size - 1");

        // Build fixed leg
        legs_[0] = FixedRateLeg(fixedSchedule_)
                       .withNotionals(fixedNominal_)
                       .withCouponRates(fixedRate_, fixedDayCount_)
                       .withPaymentAdjustment(paymentConvention_);

        // Build RA leg as simple fixed-rate placeholders.
        // The actual RA coupon values are computed by the pricing engine
        // using the observation index and range triggers.  We store
        // zero-rate coupons here so that the leg has the right structure
        // (dates, nominals) for setupArguments.
        legs_[1] = FixedRateLeg(raSchedule_)
                       .withNotionals(raNominal_)
                       .withCouponRates(std::vector<Real>(nRA, 0.0), raDayCount_)
                       .withPaymentAdjustment(paymentConvention_);

        // Build per-period observation date schedules
        Calendar obsCalendar = raSchedule_.calendar();
        observationDates_.resize(nRA);
        for (Size i = 0; i < nRA; ++i) {
            Date start = raSchedule_.date(i);
            Date end = raSchedule_.date(i + 1);
            Schedule obsSched = MakeSchedule()
                                    .from(start)
                                    .to(end)
                                    .withTenor(observationTenor_)
                                    .withCalendar(obsCalendar)
                                    .withConvention(observationConvention_)
                                    .forwards();
            // Observation dates are the interior dates (exclude the end date
            // of the accrual period itself to avoid double-counting).
            observationDates_[i].clear();
            for (Size j = 0; j < obsSched.size() - 1; ++j) {
                observationDates_[i].push_back(obsSched.date(j));
            }
            QL_REQUIRE(!observationDates_[i].empty(),
                       "No observation dates for RA coupon period " << i);
        }

        switch (type_) {
        case Swap::Payer:
            payer_[0] = +1.0;   // pays fixed
            payer_[1] = -1.0;   // receives RA
            break;
        case Swap::Receiver:
            payer_[0] = -1.0;   // receives fixed
            payer_[1] = +1.0;   // pays RA
            break;
        default:
            QL_FAIL("Unknown swap type");
        }
    }

    void RangeAccrualSwap::setupArguments(PricingEngine::arguments* args) const {

        Swap::setupArguments(args);

        auto* arguments = dynamic_cast<RangeAccrualSwap::arguments*>(args);
        if (arguments == nullptr)
            return;

        arguments->type = type_;

        // Fixed leg
        const Leg& fixedCoupons = fixedLeg();
        Size nFixed = fixedCoupons.size();
        arguments->fixedNominal = fixedNominal_;
        arguments->fixedRate = fixedRate_;
        arguments->fixedResetDates.resize(nFixed);
        arguments->fixedPayDates.resize(nFixed);
        arguments->fixedCoupons.resize(nFixed);

        for (Size i = 0; i < nFixed; ++i) {
            auto coupon =
                ext::dynamic_pointer_cast<FixedRateCoupon>(fixedCoupons[i]);
            QL_REQUIRE(coupon, "expected FixedRateCoupon on fixed leg");
            arguments->fixedResetDates[i] = coupon->accrualStartDate();
            arguments->fixedPayDates[i] = coupon->date();
            arguments->fixedCoupons[i] = coupon->amount();
        }

        // RA leg
        const Leg& raCoupons = raLeg();
        Size nRA = raCoupons.size();
        arguments->raNominal = raNominal_;
        arguments->raGearings = raGearings_;
        arguments->raSpreads = raSpreads_;
        arguments->lowerTriggers = lowerTriggers_;
        arguments->upperTriggers = upperTriggers_;
        arguments->raResetDates.resize(nRA);
        arguments->raPayDates.resize(nRA);
        arguments->raAccrualTimes.resize(nRA);

        for (Size i = 0; i < nRA; ++i) {
            auto coupon =
                ext::dynamic_pointer_cast<FixedRateCoupon>(raCoupons[i]);
            QL_REQUIRE(coupon, "expected coupon on RA leg");
            arguments->raResetDates[i] = coupon->accrualStartDate();
            arguments->raPayDates[i] = coupon->date();
            arguments->raAccrualTimes[i] = coupon->accrualPeriod();
        }

        arguments->observationDates = observationDates_;
        arguments->observationIndex = observationIndex_;
        arguments->observationTenor = observationTenor_;
    }

    void RangeAccrualSwap::setupExpired() const { Swap::setupExpired(); }

    void RangeAccrualSwap::fetchResults(
        const PricingEngine::results* r) const {
        Swap::fetchResults(r);
    }

    void RangeAccrualSwap::arguments::validate() const {
        Swap::arguments::validate();
        QL_REQUIRE(fixedNominal.size() == fixedPayDates.size(),
                   "fixed nominal size != fixed pay dates size");
        QL_REQUIRE(fixedRate.size() == fixedPayDates.size(),
                   "fixed rate size != fixed pay dates size");
        QL_REQUIRE(raNominal.size() == raPayDates.size(),
                   "RA nominal size != RA pay dates size");
        QL_REQUIRE(raGearings.size() == raPayDates.size(),
                   "RA gearings size != RA pay dates size");
        QL_REQUIRE(raSpreads.size() == raPayDates.size(),
                   "RA spreads size != RA pay dates size");
        QL_REQUIRE(lowerTriggers.size() == raPayDates.size(),
                   "lower triggers size != RA pay dates size");
        QL_REQUIRE(upperTriggers.size() == raPayDates.size(),
                   "upper triggers size != RA pay dates size");
        QL_REQUIRE(observationDates.size() == raPayDates.size(),
                   "observation dates size != RA pay dates size");
        QL_REQUIRE(observationIndex != nullptr,
                   "observation index not set");
    }

    void RangeAccrualSwap::results::reset() { Swap::results::reset(); }
}
