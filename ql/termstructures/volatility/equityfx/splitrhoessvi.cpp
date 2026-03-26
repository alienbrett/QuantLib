/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Chloride Project
*/

#include <ql/termstructures/volatility/equityfx/splitrhoessvi.hpp>

namespace QuantLib {

    namespace {
        std::vector<Real> datesToTimes(const Date& ref,
                                        const std::vector<Date>& dates,
                                        const DayCounter& dc) {
            std::vector<Real> t(dates.size());
            for (Size i = 0; i < dates.size(); ++i) {
                t[i] = dc.yearFraction(ref, dates[i]);
                QL_REQUIRE(t[i] > 0.0, "date[" << i << "] must be after reference date");
            }
            return t;
        }

        std::vector<SplitRhoEssviSliceParams> buildSlices(
                const std::vector<Real>& thetas,
                const std::vector<Real>& rhos_lo,
                const std::vector<Real>& psis_lo,
                const std::vector<Real>& rhos_hi,
                const std::vector<Real>& psis_hi) {
            Size N = thetas.size();
            QL_REQUIRE(rhos_lo.size() == N && psis_lo.size() == N
                       && rhos_hi.size() == N && psis_hi.size() == N,
                       "parameter vector sizes must match");
            std::vector<SplitRhoEssviSliceParams> s(N);
            for (Size i = 0; i < N; ++i) {
                s[i].theta  = thetas[i];
                s[i].rho_lo = rhos_lo[i];
                s[i].psi_lo = psis_lo[i];
                s[i].rho_hi = rhos_hi[i];
                s[i].psi_hi = psis_hi[i];
            }
            return s;
        }
    }

    // Native params constructor
    SplitRhoEssviVolatilityTermStructure::SplitRhoEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos_lo,
            const std::vector<Real>& psis_lo,
            const std::vector<Real>& rhos_hi,
            const std::vector<Real>& psis_hi,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          surface_(datesToTimes(referenceDate, dates, dc),
                   buildSlices(thetas, rhos_lo, psis_lo, rhos_hi, psis_hi)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield))
    {
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    // Global params constructor
    SplitRhoEssviVolatilityTermStructure::SplitRhoEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos_lo,
            const std::vector<Real>& rhos_hi,
            Real theta1,
            const std::vector<Real>& as,
            const std::vector<Real>& cs_lo,
            const std::vector<Real>& cs_hi,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            EssviButterflyCondition::Type bflyType,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          surface_(datesToTimes(referenceDate, dates, dc),
                   [&]() {
                       SplitRhoEssviGlobalParams gp;
                       gp.rhos_lo = rhos_lo;
                       gp.rhos_hi = rhos_hi;
                       gp.theta1  = theta1;
                       gp.as      = as;
                       gp.cs_lo   = cs_lo;
                       gp.cs_hi   = cs_hi;
                       return gp;
                   }(),
                   bflyType),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield))
    {
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    Date SplitRhoEssviVolatilityTermStructure::maxDate() const { return Date::maxDate(); }
    Real SplitRhoEssviVolatilityTermStructure::minStrike() const { return QL_MIN_REAL; }
    Real SplitRhoEssviVolatilityTermStructure::maxStrike() const { return QL_MAX_REAL; }
    Size SplitRhoEssviVolatilityTermStructure::numSlices() const { return surface_.numSlices(); }

    const std::vector<SplitRhoEssviSliceParams>&
    SplitRhoEssviVolatilityTermStructure::slices() const { return surface_.slices(); }

    void SplitRhoEssviVolatilityTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<SplitRhoEssviVolatilityTermStructure>*>(&v);
        if (v1) v1->visit(*this);
        else BlackVolatilityTermStructure::accept(v);
    }

    SplitRhoEssviSliceGradient
    SplitRhoEssviVolatilityTermStructure::impliedVolGradient(Size sliceIdx, Real strike) const {
        Time t = surface_.maturities().at(sliceIdx);
        Real k = std::log(strike / forward(t));
        return surface_.impliedVolGradient(sliceIdx, k);
    }

    std::vector<Real>
    SplitRhoEssviVolatilityTermStructure::impliedVolGlobalGradient(
            Size sliceIdx, Real strike,
            const SplitRhoEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {
        Time t = surface_.maturities().at(sliceIdx);
        Real k = std::log(strike / forward(t));
        return surface_.impliedVolGlobalGradient(sliceIdx, k, gp, bflyCond);
    }

    Real SplitRhoEssviVolatilityTermStructure::forward(Time t) const {
        return spot_->value() * dividendYield_->discount(t) / riskFreeRate_->discount(t);
    }

    Volatility SplitRhoEssviVolatilityTermStructure::blackVolImpl(Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;
        Real k = std::log(strike / forward(t));
        return surface_.impliedVol(k, t);
    }

} // namespace QuantLib
