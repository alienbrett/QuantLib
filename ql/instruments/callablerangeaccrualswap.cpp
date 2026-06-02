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

#include <ql/instruments/callablerangeaccrualswap.hpp>
#include <ql/exercise.hpp>
#include <utility>

namespace QuantLib {

    CallableRangeAccrualSwap::CallableRangeAccrualSwap(
        ext::shared_ptr<RangeAccrualSwap> swap,
        const ext::shared_ptr<Exercise>& exercise,
        Settlement::Type delivery,
        Settlement::Method settlementMethod,
        Real callPrice)
    : Option(ext::shared_ptr<Payoff>(), exercise),
      swap_(std::move(swap)),
      settlementType_(delivery),
      settlementMethod_(settlementMethod),
      callPrice_(callPrice) {

        registerWith(swap_);
        swap_->alwaysForwardNotifications();
    }

    bool CallableRangeAccrualSwap::isExpired() const {
        return detail::simple_event(exercise_->dates().back()).hasOccurred();
    }

    void CallableRangeAccrualSwap::setupArguments(
        PricingEngine::arguments* args) const {

        swap_->setupArguments(args);

        auto* arguments =
            dynamic_cast<CallableRangeAccrualSwap::arguments*>(args);
        QL_REQUIRE(arguments != nullptr, "argument types do not match");

        arguments->swap = swap_;
        arguments->exercise = exercise_;
        arguments->settlementType = settlementType_;
        arguments->settlementMethod = settlementMethod_;
        arguments->callPrice = callPrice_;
    }

    void CallableRangeAccrualSwap::arguments::validate() const {
        RangeAccrualSwap::arguments::validate();
        QL_REQUIRE(swap, "underlying range accrual swap not set");
        QL_REQUIRE(exercise, "exercise not set");
        Settlement::checkTypeAndMethodConsistency(settlementType,
                                                  settlementMethod);
    }
}
