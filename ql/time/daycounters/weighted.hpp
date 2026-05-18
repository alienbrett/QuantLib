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

/*! \file weighted.hpp
    \brief variance-time day counter driven by a precomputed cumulative
           weight table
*/

#ifndef quantlib_weighted_day_counter_hpp
#define quantlib_weighted_day_counter_hpp

#include <ql/time/daycounter.hpp>
#include <string>
#include <utility>
#include <vector>

namespace QuantLib {

    //! Variance-time day counter
    /*! Maps calendar dates to "variance time" via a precomputed
        cumulative-weight array. ``yearFraction(d1, d2)`` is
        ``(cum[d2 - anchor] - cum[d1 - anchor]) / annual_weight``.

        Typical use: the caller (Python) fits per-day variance weights
        from realized returns (e.g. weekend << NYSE session, FOMC days
        elevated), then materialises a cumulative table over a horizon
        that comfortably covers all option expiries. QuantLib pricing
        engines then see this day counter and integrate volatility in
        variance time rather than calendar time.

        \ingroup daycounters
    */
    class WeightedDayCounter : public DayCounter {
      private:
        class Impl final : public DayCounter::Impl {
          private:
            Date anchor_;
            std::vector<Real> cumWeights_;
            Real annualWeight_;
            std::string name_;
          public:
            Impl(Date anchor,
                 std::vector<Real> cumWeights,
                 Real annualWeight,
                 std::string name);
            std::string name() const override { return name_; }
            Date::serial_type dayCount(const Date& d1, const Date& d2) const override;
            Time yearFraction(const Date& d1, const Date& d2,
                              const Date&, const Date&) const override;
        };
      public:
        WeightedDayCounter(const Date& anchor,
                           const std::vector<Real>& cumWeights,
                           Real annualWeight,
                           const std::string& name = "WeightedDayCounter")
        : DayCounter(ext::shared_ptr<DayCounter::Impl>(
              new WeightedDayCounter::Impl(anchor, cumWeights, annualWeight, name))) {}
    };

}

#endif
