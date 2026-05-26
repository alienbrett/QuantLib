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

#include <ql/termstructures/volatility/equityfx/mixedspotlocalvoltermstructure.hpp>
#include <ql/errors.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    namespace {
        constexpr Real kSigmaFFloor   = 1e-12;
        constexpr Real kTimeFloor     = 1e-14;
        constexpr Real kVolFloor      = 1e-6;
        constexpr Real kDensityFloor  = 1e-300;

        // Black-Scholes log-normal density at S=K under fwd Fᵢ, vol σᵢ, T.
        // p(K) = exp(-d2²/2) / (K · σ · √(2π·T))
        Real lognormalDensity(Real K, Real F_i, Real sigma_i, Time T) {
            if (K <= 0.0 || F_i <= 0.0 || sigma_i < kVolFloor || T < kTimeFloor)
                return 0.0;
            Real stddev = sigma_i * std::sqrt(T);
            Real d2 = (std::log(F_i / K) - 0.5 * stddev * stddev) / stddev;
            Real pdf = std::exp(-0.5 * d2 * d2) /
                       (K * stddev * std::sqrt(2.0 * M_PI));
            return pdf;
        }
    }

    MixedSpotLocalVolTermStructure::MixedSpotLocalVolTermStructure(
            Handle<LocalVolTermStructure> baseLV,
            Handle<BlackVolTermStructure> baseBlack,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau)
        : LocalVolTermStructure(
              baseLV->referenceDate(),
              baseLV->calendar(),
              baseLV->businessDayConvention(),
              baseLV->dayCounter()),
          baseLV_(std::move(baseLV)),
          baseBlack_(std::move(baseBlack)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          sigmaF_(sigmaF),
          fwdDecayTau_(fwdDecayTau),
          nodes_({-1.0, 1.0}),
          weights_({0.5, 0.5}) {
        validate_();
        registerWith(baseLV_);
        registerWith(baseBlack_);
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    MixedSpotLocalVolTermStructure::MixedSpotLocalVolTermStructure(
            Handle<LocalVolTermStructure> baseLV,
            Handle<BlackVolTermStructure> baseBlack,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau,
            std::vector<Real> quadNodes,
            std::vector<Real> quadWeights)
        : LocalVolTermStructure(
              baseLV->referenceDate(),
              baseLV->calendar(),
              baseLV->businessDayConvention(),
              baseLV->dayCounter()),
          baseLV_(std::move(baseLV)),
          baseBlack_(std::move(baseBlack)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          sigmaF_(sigmaF),
          fwdDecayTau_(fwdDecayTau),
          nodes_(std::move(quadNodes)),
          weights_(std::move(quadWeights)) {
        validate_();
        registerWith(baseLV_);
        registerWith(baseBlack_);
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    void MixedSpotLocalVolTermStructure::validate_() const {
        QL_REQUIRE(!baseLV_.empty(), "base LV handle must not be empty");
        QL_REQUIRE(!baseBlack_.empty(), "base Black handle must not be empty");
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(),
                   "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(),
                   "dividend yield handle must not be empty");
        QL_REQUIRE(fwdDecayTau_ >= 0.0,
                   "fwdDecayTau must be >= 0, got " << fwdDecayTau_);
        QL_REQUIRE(nodes_.size() >= 2,
                   "quadrature must have >= 2 nodes");
        QL_REQUIRE(nodes_.size() == weights_.size(),
                   "quadrature nodes/weights size mismatch");
        Real ws = 0.0;
        for (Real w : weights_) {
            QL_REQUIRE(w > 0.0, "quadrature weights must be > 0");
            ws += w;
        }
        QL_REQUIRE(std::abs(ws - 1.0) < 1e-12,
                   "quadrature weights must sum to 1");
    }

    Date MixedSpotLocalVolTermStructure::maxDate() const {
        return baseLV_->maxDate();
    }

    Real MixedSpotLocalVolTermStructure::minStrike() const {
        return baseLV_->minStrike();
    }

    Real MixedSpotLocalVolTermStructure::maxStrike() const {
        return baseLV_->maxStrike();
    }

    void MixedSpotLocalVolTermStructure::setSigmaF(Real sigmaF) {
        sigmaF_ = sigmaF;
        notifyObservers();
    }

    void MixedSpotLocalVolTermStructure::setFwdDecayTau(Real tau) {
        QL_REQUIRE(tau >= 0.0, "fwdDecayTau must be >= 0");
        fwdDecayTau_ = tau;
        notifyObservers();
    }

    void MixedSpotLocalVolTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 =
            dynamic_cast<Visitor<MixedSpotLocalVolTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            LocalVolTermStructure::accept(v);
    }

    Real MixedSpotLocalVolTermStructure::forward_(Time t) const {
        Real S  = spot_->value();
        Real Dr = riskFreeRate_->discount(t);
        Real Dq = dividendYield_->discount(t);
        QL_REQUIRE(Dr > 0.0, "risk-free discount must be > 0");
        return S * Dq / Dr;
    }

    Real MixedSpotLocalVolTermStructure::sigmaFAt_(Time t) const {
        if (fwdDecayTau_ > kSigmaFFloor) {
            return sigmaF_ * std::exp(-t / fwdDecayTau_);
        }
        return sigmaF_;
    }

    Volatility MixedSpotLocalVolTermStructure::localVolImpl(
            Time t, Real strike) const {
        if (t < kTimeFloor) t = kTimeFloor;

        Real sf = sigmaFAt_(t);

        // sigma_F == 0 ⇒ identity overlay; pass through base directly
        if (std::abs(sf) < kSigmaFFloor) {
            return baseLV_->localVol(t, strike, true);
        }

        Real F = forward_(t);
        Real norm = std::cosh(sf);

        // Density-weighted local vol:
        //   σ_LV,mix²(K,T) = (Σᵢ wᵢ · ρᵢ(K) · σ_LV,base²(K·F/Fᵢ, T))
        //                 / (Σᵢ wᵢ · ρᵢ(K))
        Real sumNum = 0.0, sumDen = 0.0;
        for (Size i = 0; i < nodes_.size(); ++i) {
            Real F_i = F * std::exp(sf * nodes_[i]) / norm;
            Real K_q = strike * F / F_i;

            // Branch i's BS vol at the SYNTHETIC strike — same as overlay's
            // internal query.
            Real vol_i = baseBlack_->blackVol(t, K_q, true);
            if (vol_i < kVolFloor) continue;

            Real lv_i = baseLV_->localVol(t, K_q, true);
            if (lv_i < kVolFloor) continue;

            Real density_i = lognormalDensity(strike, F_i, vol_i, t);
            if (density_i < kDensityFloor) continue;

            Real w_rho = weights_[i] * density_i;
            sumNum += w_rho * lv_i * lv_i;
            sumDen += w_rho;
        }

        if (sumDen < kDensityFloor) {
            // Degenerate: fall back to base LV at the true strike
            return baseLV_->localVol(t, strike, true);
        }

        return std::sqrt(sumNum / sumDen);
    }

} // namespace QuantLib
