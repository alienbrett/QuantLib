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

/*! \file callablerangeaccrualswap.hpp
    \brief callable (Bermudan) option on a range accrual swap
*/

#ifndef quantlib_callable_range_accrual_swap_hpp
#define quantlib_callable_range_accrual_swap_hpp

#include <ql/option.hpp>
#include <ql/instruments/rangeaccrualswap.hpp>
#include <ql/instruments/swaption.hpp>

namespace QuantLib {

    //! Callable range accrual swap
    /*! Option (typically Bermudan) on a RangeAccrualSwap.

        The holder has the right to terminate (call) the swap at each
        exercise date.  Exercise is physical settlement: on exercise
        the remaining swap cash flows cease.

        Follows the same pattern as NonstandardSwaption wrapping
        NonstandardSwap.

        \ingroup instruments
    */

    class CallableRangeAccrualSwap : public Option {
      public:
        class arguments;
        class engine;

        CallableRangeAccrualSwap(
            ext::shared_ptr<RangeAccrualSwap> swap,
            const ext::shared_ptr<Exercise>& exercise,
            Settlement::Type delivery = Settlement::Physical,
            Settlement::Method settlementMethod = Settlement::PhysicalOTC,
            Real callPrice = 0.0);

        //! \name Instrument interface
        //@{
        bool isExpired() const override;
        void setupArguments(PricingEngine::arguments*) const override;
        //@}

        //! \name Inspectors
        //@{
        Settlement::Type settlementType() const { return settlementType_; }
        Settlement::Method settlementMethod() const {
            return settlementMethod_;
        }
        Swap::Type type() const { return swap_->type(); }

        const ext::shared_ptr<RangeAccrualSwap>& underlyingSwap() const {
            return swap_;
        }
        //@}

        Real callPrice() const { return callPrice_; }

      private:
        ext::shared_ptr<RangeAccrualSwap> swap_;
        Settlement::Type settlementType_;
        Settlement::Method settlementMethod_;
        Real callPrice_;
    };

    //! %Arguments for callable range accrual swap calculation
    class CallableRangeAccrualSwap::arguments
        : public RangeAccrualSwap::arguments,
          public Option::arguments {
      public:
        arguments() = default;
        ext::shared_ptr<RangeAccrualSwap> swap;
        Settlement::Type settlementType = Settlement::Physical;
        Settlement::Method settlementMethod = Settlement::PhysicalOTC;
        Real callPrice = 0.0;
        void validate() const override;
    };

    //! base class for callable range accrual swap engines
    class CallableRangeAccrualSwap::engine
        : public GenericEngine<CallableRangeAccrualSwap::arguments,
                               CallableRangeAccrualSwap::results> {};
}

#endif
