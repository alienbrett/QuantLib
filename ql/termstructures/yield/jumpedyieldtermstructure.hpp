/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Jumped yield term structure — multiplicative DF adjustments at
 discrete dates.

 Used for proportional equity dividends: at each ex-date the forward
 drops by (1 − q), modeled as a step in the discount factor.

 DF(t) = base_DF(t) × prod(factor_i  for  jumpTime_i ≤ t)

 Forwards between jumps are preserved exactly because QL computes them
 from the discount function.
*/

#ifndef quantlib_jumped_yield_term_structure_hpp
#define quantlib_jumped_yield_term_structure_hpp

#include <ql/termstructures/yieldtermstructure.hpp>
#include <algorithm>
#include <utility>
#include <vector>

namespace QuantLib {

    //! Yield term structure with multiplicative discount factor jumps
    /*! At specified times, the discount factor is multiplied by a jump
        factor (typically 1 − q for a proportional dividend q).

        The base curve's forwards are preserved exactly between jumps.

        \note This term structure remains linked to the base — any
              changes in the base are reflected here.

        \ingroup yieldtermstructures
    */
    class JumpedYieldTermStructure : public YieldTermStructure {
      public:
        /*! \param baseHandle      Base yield term structure.
            \param jumpTimes       Times (year fractions from base ref date)
                                   at which jumps occur. Must be sorted.
            \param jumpFactors     Multiplicative factor at each jump time
                                   (e.g. 0.98 for a 2% proportional div).
        */
        JumpedYieldTermStructure(
            Handle<YieldTermStructure> baseHandle,
            std::vector<Time> jumpTimes,
            std::vector<DiscountFactor> jumpFactors);

        //! \name YieldTermStructure interface
        //@{
        DayCounter dayCounter() const override;
        Calendar calendar() const override;
        Natural settlementDays() const override;
        Date maxDate() const override;
        //@}

      protected:
        DiscountFactor discountImpl(Time) const override;

      private:
        Handle<YieldTermStructure> base_;
        std::vector<Time> jumpTimes_;
        // Cumulative product of jump factors: cumFactors_[i] =
        // prod(jumpFactors_[0..i]).  Used for O(log n) lookup.
        std::vector<DiscountFactor> cumFactors_;
    };


    // inline definitions

    inline JumpedYieldTermStructure::JumpedYieldTermStructure(
        Handle<YieldTermStructure> baseHandle,
        std::vector<Time> jumpTimes,
        std::vector<DiscountFactor> jumpFactors)
    : YieldTermStructure(baseHandle->referenceDate()),
      base_(std::move(baseHandle)),
      jumpTimes_(std::move(jumpTimes))
    {
        QL_REQUIRE(jumpTimes_.size() == jumpFactors.size(),
                   "jumpTimes size (" << jumpTimes_.size()
                   << ") must match jumpFactors size ("
                   << jumpFactors.size() << ")");

        // Build cumulative factor array
        cumFactors_.resize(jumpTimes_.size());
        DiscountFactor cum = 1.0;
        for (Size i = 0; i < jumpFactors.size(); ++i) {
            cum *= jumpFactors[i];
            cumFactors_[i] = cum;
        }

        registerWith(base_);
    }

    inline DayCounter JumpedYieldTermStructure::dayCounter() const {
        return base_->dayCounter();
    }

    inline Calendar JumpedYieldTermStructure::calendar() const {
        return base_->calendar();
    }

    inline Natural JumpedYieldTermStructure::settlementDays() const {
        return base_->settlementDays();
    }

    inline Date JumpedYieldTermStructure::maxDate() const {
        return base_->maxDate();
    }

    inline DiscountFactor JumpedYieldTermStructure::discountImpl(Time t) const {
        // Find the cumulative jump factor for all jumps at or before t.
        // jumpTimes_ is sorted; upper_bound gives first element > t.
        auto it = std::upper_bound(jumpTimes_.begin(), jumpTimes_.end(), t);
        DiscountFactor jumpFactor = 1.0;
        if (it != jumpTimes_.begin()) {
            jumpFactor = cumFactors_[std::distance(jumpTimes_.begin(), it) - 1];
        }
        return base_->discount(t, true) * jumpFactor;
    }

}

#endif
