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

#include <ql/termstructures/volatility/equityfx/dualwingessvi.hpp>

namespace QuantLib {

    namespace {
        std::vector<Real> datesToTimes(const Date& referenceDate,
                                        const std::vector<Date>& dates,
                                        const DayCounter& dc) {
            std::vector<Real> times(dates.size());
            for (Size i = 0; i < dates.size(); ++i) {
                times[i] = dc.yearFraction(referenceDate, dates[i]);
                QL_REQUIRE(times[i] > 0.0,
                           "date[" << i << "] must be after reference date");
            }
            return times;
        }

        std::vector<DualWingEssviSliceParams> buildDualWingSlices(
                const std::vector<Real>& thetas,
                const std::vector<Real>& rhos,
                const std::vector<Real>& psis_lo,
                const std::vector<Real>& psis_hi) {
            Size N = thetas.size();
            QL_REQUIRE(rhos.size() == N && psis_lo.size() == N && psis_hi.size() == N,
                       "parameter vector sizes must match");
            std::vector<DualWingEssviSliceParams> slices(N);
            for (Size i = 0; i < N; ++i) {
                slices[i].theta  = thetas[i];
                slices[i].rho    = rhos[i];
                slices[i].psi_lo = psis_lo[i];
                slices[i].psi_hi = psis_hi[i];
            }
            return slices;
        }
    }

    // Native params constructor
    DualWingEssviVolatilityTermStructure::DualWingEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos,
            const std::vector<Real>& psis_lo,
            const std::vector<Real>& psis_hi,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          surface_(datesToTimes(referenceDate, dates, dc),
                   buildDualWingSlices(thetas, rhos, psis_lo, psis_hi)),
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
    DualWingEssviVolatilityTermStructure::DualWingEssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos,
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
                       DualWingEssviGlobalParams gp;
                       gp.rhos   = rhos;
                       gp.theta1 = theta1;
                       gp.as     = as;
                       gp.cs_lo  = cs_lo;
                       gp.cs_hi  = cs_hi;
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

    Date DualWingEssviVolatilityTermStructure::maxDate() const {
        return Date::maxDate();
    }

    Real DualWingEssviVolatilityTermStructure::minStrike() const {
        return QL_MIN_REAL;
    }

    Real DualWingEssviVolatilityTermStructure::maxStrike() const {
        return QL_MAX_REAL;
    }

    Size DualWingEssviVolatilityTermStructure::numSlices() const {
        return surface_.numSlices();
    }

    const std::vector<DualWingEssviSliceParams>&
    DualWingEssviVolatilityTermStructure::slices() const {
        return surface_.slices();
    }

    void DualWingEssviVolatilityTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<DualWingEssviVolatilityTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

    DualWingEssviSliceGradient
    DualWingEssviVolatilityTermStructure::impliedVolGradient(
            Size sliceIdx, Real strike) const {
        Time t = surface_.maturities().at(sliceIdx);
        Real fwd = forward(t);
        Real k = std::log(strike / fwd);
        return surface_.impliedVolGradient(sliceIdx, k);
    }

    std::vector<Real>
    DualWingEssviVolatilityTermStructure::impliedVolGlobalGradient(
            Size sliceIdx, Real strike,
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {
        Time t = surface_.maturities().at(sliceIdx);
        Real fwd = forward(t);
        Real k = std::log(strike / fwd);
        return surface_.impliedVolGlobalGradient(sliceIdx, k, gp, bflyCond);
    }

    Real DualWingEssviVolatilityTermStructure::forward(Time t) const {
        Real S = spot_->value();
        Real df = riskFreeRate_->discount(t);
        Real dq = dividendYield_->discount(t);
        return S * dq / df;
    }

    Volatility DualWingEssviVolatilityTermStructure::blackVolImpl(
            Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;
        Real fwd = forward(t);
        Real k = std::log(strike / fwd);
        return surface_.impliedVol(k, t);
    }

    std::vector<Real>
    DualWingEssviVolatilityTermStructure::batchBlackVol(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes) const {
        const Size nObs = sliceIndices.size();
        QL_REQUIRE(strikes.size() == nObs, "strikes size mismatch");

        const auto& T = surface_.maturities();
        const Size nSlices = T.size();
        std::vector<Real> fwds(nSlices);
        for (Size s = 0; s < nSlices; ++s)
            fwds[s] = forward(T[s]);

        std::vector<Real> vols(nObs);
        for (Size i = 0; i < nObs; ++i) {
            Size si = sliceIndices[i];
            Real k = std::log(strikes[i] / fwds[si]);
            Real w = surface_.totalVariance(si, k);
            vols[i] = std::sqrt(w / T[si]);
        }
        return vols;
    }

    std::vector<Real>
    DualWingEssviVolatilityTermStructure::batchImpliedVolGlobalGradient(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes,
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {
        const Size nObs = sliceIndices.size();
        QL_REQUIRE(strikes.size() == nObs, "strikes size mismatch");
        const Size N = gp.numSlices();
        const Size nParams = 4 * N;

        const auto& T = surface_.maturities();
        std::vector<Real> fwds(N);
        for (Size s = 0; s < N; ++s)
            fwds[s] = forward(T[s]);

        // Factored: compute chain Jacobian once
        // Layout: [dTheta(N×P), dPsiLo(N×P), dPsiHi(N×P)]
        auto chainJac = surface_.chainJacobian(gp, bflyCond);
        const Real* dTheta = chainJac.data();
        const Real* dPsiLo = chainJac.data() + N * nParams;
        const Real* dPsiHi = chainJac.data() + 2 * N * nParams;
        Size rhoOff = 0;

        std::vector<Real> J(nObs * nParams, 0.0);
        for (Size i = 0; i < nObs; ++i) {
            Size si = sliceIndices[i];
            Real k = std::log(strikes[i] / fwds[si]);
            DualWingEssviSliceGradient sg =
                dualWingEssviImpliedVolGradient(surface_.slice(si), k, T[si]);

            Real* row = &J[i * nParams];
            const Real* dTh_row = &dTheta[si * nParams];
            // Active wing's psi chain
            Real dSigma_dPsi;
            const Real* dPs_row;
            if (k < 0.0) {
                dSigma_dPsi = sg.dSigma_dPsiLo;
                dPs_row = &dPsiLo[si * nParams];
            } else {
                dSigma_dPsi = sg.dSigma_dPsiHi;
                dPs_row = &dPsiHi[si * nParams];
            }

            for (Size j = 0; j < nParams; ++j) {
                row[j] = sg.dSigma_dTheta * dTh_row[j]
                       + dSigma_dPsi * dPs_row[j];
            }
            row[rhoOff + si] += sg.dSigma_dRho;
        }
        return J;
    }

} // namespace QuantLib
