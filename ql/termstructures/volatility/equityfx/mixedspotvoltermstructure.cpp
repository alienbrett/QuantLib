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

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/termstructures/volatility/equityfx/mixedspotvoltermstructure.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/instruments/payoffs.hpp>
#include <ql/errors.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    namespace {
        constexpr Real kSigmaFFloor   = 1e-12;
        constexpr Real kTimeFloor     = 1e-14;
        constexpr Real kMinVolReturn  = 1e-6;
        constexpr Real kMaxVolReturn  = 10.0;
        constexpr Real kIvAccuracy    = 1e-10;
        constexpr Natural kIvMaxIter  = 200;
    }

    MixedSpotVolTermStructure::MixedSpotVolTermStructure(
            Handle<BlackVolTermStructure> base,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau)
        : BlackVolatilityTermStructure(
              base->referenceDate(),
              base->calendar(),
              base->businessDayConvention(),
              base->dayCounter()),
          base_(std::move(base)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          sigmaF_(sigmaF),
          fwdDecayTau_(fwdDecayTau),
          nodes_({-1.0, 1.0}),
          weights_({0.5, 0.5}) {
        validate_();
        registerWith(base_);
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    MixedSpotVolTermStructure::MixedSpotVolTermStructure(
            Handle<BlackVolTermStructure> base,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau,
            std::vector<Real> quadNodes,
            std::vector<Real> quadWeights)
        : BlackVolatilityTermStructure(
              base->referenceDate(),
              base->calendar(),
              base->businessDayConvention(),
              base->dayCounter()),
          base_(std::move(base)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          sigmaF_(sigmaF),
          fwdDecayTau_(fwdDecayTau),
          nodes_(std::move(quadNodes)),
          weights_(std::move(quadWeights)) {
        validate_();
        registerWith(base_);
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    void MixedSpotVolTermStructure::validate_() const {
        QL_REQUIRE(!base_.empty(), "base vol handle must not be empty");
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(),
                   "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(),
                   "dividend yield handle must not be empty");
        QL_REQUIRE(fwdDecayTau_ >= 0.0,
                   "fwdDecayTau must be >= 0 (0 disables decay), got "
                       << fwdDecayTau_);
        QL_REQUIRE(nodes_.size() >= 2,
                   "quadrature must have >= 2 nodes, got " << nodes_.size());
        QL_REQUIRE(nodes_.size() == weights_.size(),
                   "quadrature nodes (" << nodes_.size()
                   << ") and weights (" << weights_.size()
                   << ") must have equal length");
        Real ws = 0.0;
        for (Real w : weights_) {
            QL_REQUIRE(w > 0.0, "quadrature weights must be > 0");
            ws += w;
        }
        QL_REQUIRE(std::abs(ws - 1.0) < 1e-12,
                   "quadrature weights must sum to 1, got " << ws);
    }

    Date MixedSpotVolTermStructure::maxDate() const {
        return base_->maxDate();
    }

    Real MixedSpotVolTermStructure::minStrike() const {
        return base_->minStrike();
    }

    Real MixedSpotVolTermStructure::maxStrike() const {
        return base_->maxStrike();
    }

    void MixedSpotVolTermStructure::setSigmaF(Real sigmaF) {
        sigmaF_ = sigmaF;
        notifyObservers();
    }

    void MixedSpotVolTermStructure::setFwdDecayTau(Real tau) {
        QL_REQUIRE(tau >= 0.0, "fwdDecayTau must be >= 0, got " << tau);
        fwdDecayTau_ = tau;
        notifyObservers();
    }

    void MixedSpotVolTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<MixedSpotVolTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

    Real MixedSpotVolTermStructure::forward_(Time t) const {
        Real S  = spot_->value();
        Real Dr = riskFreeRate_->discount(t);
        Real Dq = dividendYield_->discount(t);
        QL_REQUIRE(Dr > 0.0, "risk-free discount must be > 0");
        return S * Dq / Dr;
    }

    Real MixedSpotVolTermStructure::sigmaFAt_(Time t) const {
        if (fwdDecayTau_ > kSigmaFFloor) {
            return sigmaF_ * std::exp(-t / fwdDecayTau_);
        }
        return sigmaF_;
    }

    Volatility MixedSpotVolTermStructure::blackVolImpl(Time t, Real strike) const {
        // Degenerate-T floor (mirrors EssviVolatilityTermStructure)
        if (t < kTimeFloor) t = kTimeFloor;

        Real sf = sigmaFAt_(t);

        // sigma_F == 0 ⇒ identity overlay; pass through base directly
        if (std::abs(sf) < kSigmaFFloor) {
            return base_->blackVol(t, strike, true);
        }

        Real F  = forward_(t);
        Real Dr = riskFreeRate_->discount(t);
        Real norm = std::cosh(sf);

        // Build mixture price by 2-point quadrature on the spot factor
        Real priceMix = 0.0;
        for (Size i = 0; i < nodes_.size(); ++i) {
            Real F_i = F * std::exp(sf * nodes_[i]) / norm;
            Real K_q = strike * F / F_i;   // base sees log(K_q/F) = log(K/F_i)
            Real vol_i = base_->blackVol(t, K_q, true);
            Real stddev = vol_i * std::sqrt(t);
            priceMix += weights_[i] *
                        blackFormula(Option::Call, strike, F_i, stddev, Dr);
        }

        // Invert blended call price to BS std-dev at the TRUE forward F
        Real intrinsic = std::max(F - strike, 0.0) * Dr;
        if (priceMix <= intrinsic + 1e-14) {
            // Numerical degeneracy (e.g. very deep ITM call with both
            // F_i below strike): fall back to the base vol
            return base_->blackVol(t, strike, true);
        }

        Real stddevMix;
        try {
            stddevMix = blackFormulaImpliedStdDev(
                Option::Call,
                strike,
                F,
                priceMix,
                Dr,
                0.0,             // displacement
                Null<Real>(),    // guess
                kIvAccuracy,
                kIvMaxIter);
        } catch (const std::exception&) {
            // Fall back to base if BS inversion fails (e.g. wing pathology)
            return base_->blackVol(t, strike, true);
        }

        Real vol = stddevMix / std::sqrt(t);
        // Clamp to a sensible range so downstream FD engines don't blow up
        return std::max(kMinVolReturn, std::min(kMaxVolReturn, vol));
    }

} // namespace QuantLib
