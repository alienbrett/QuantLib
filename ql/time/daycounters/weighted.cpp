/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Chloride

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

#include <ql/time/daycounters/weighted.hpp>
#include <ql/errors.hpp>
#include <utility>

namespace QuantLib {

    WeightedDayCounter::Impl::Impl(Date anchor,
                                   std::vector<Real> cumWeights,
                                   Real annualWeight,
                                   std::string name)
    : anchor_(anchor),
      cumWeights_(std::move(cumWeights)),
      annualWeight_(annualWeight),
      name_(std::move(name)) {
        QL_REQUIRE(cumWeights_.size() >= 2,
                   "WeightedDayCounter: cumWeights must have at least 2 entries");
        QL_REQUIRE(cumWeights_.front() == 0.0,
                   "WeightedDayCounter: cumWeights[0] must be 0.0");
        QL_REQUIRE(annualWeight_ > 0.0,
                   "WeightedDayCounter: annualWeight must be > 0");
    }

    Date::serial_type WeightedDayCounter::Impl::dayCount(const Date& d1,
                                                         const Date& d2) const {
        return d2 - d1;
    }

    Time WeightedDayCounter::Impl::yearFraction(const Date& d1,
                                                const Date& d2,
                                                const Date&,
                                                const Date&) const {
        if (d1 == d2) return 0.0;
        if (d1 > d2) return -yearFraction(d2, d1, Date(), Date());

        const Date::serial_type k1 = d1 - anchor_;
        const Date::serial_type k2 = d2 - anchor_;
        const Date::serial_type n = static_cast<Date::serial_type>(cumWeights_.size()) - 1;

        QL_REQUIRE(k1 >= 0 && k2 <= n,
                   "WeightedDayCounter: dates ("
                   << d1 << ", " << d2 << ") out of range [anchor="
                   << anchor_ << ", anchor+" << n << "]");

        return (cumWeights_[k2] - cumWeights_[k1]) / annualWeight_;
    }

}
