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

#include <ql/termstructures/volatility/equityfx/essvivoltermstructure.hpp>
#include <ql/cashflows/dividend.hpp>

namespace QuantLib {

    namespace {
        // Convert dates to year fractions relative to referenceDate
        std::vector<Real> datesToTimes(const Date& referenceDate,
                                        const std::vector<Date>& dates,
                                        const DayCounter& dc) {
            std::vector<Real> times(dates.size());
            for (Size i = 0; i < dates.size(); ++i) {
                times[i] = dc.yearFraction(referenceDate, dates[i]);
                QL_REQUIRE(times[i] > 0.0,
                           "date[" << i << "] (" << dates[i]
                           << ") must be after reference date ("
                           << referenceDate << ")");
            }
            return times;
        }

        // Build EssviSliceParams from parallel vectors
        std::vector<EssviSliceParams> buildSlices(
                const std::vector<Real>& thetas,
                const std::vector<Real>& rhos,
                const std::vector<Real>& psis) {
            Size N = thetas.size();
            QL_REQUIRE(rhos.size() == N,
                       "rhos size (" << rhos.size()
                       << ") must match thetas size (" << N << ")");
            QL_REQUIRE(psis.size() == N,
                       "psis size (" << psis.size()
                       << ") must match thetas size (" << N << ")");
            std::vector<EssviSliceParams> slices(N);
            for (Size i = 0; i < N; ++i) {
                slices[i].theta = thetas[i];
                slices[i].rho   = rhos[i];
                slices[i].psi   = psis[i];
            }
            return slices;
        }
    }

    // =================================================================
    // Native params constructor (continuous dividends only)
    // =================================================================

    EssviVolatilityTermStructure::EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos,
            const std::vector<Real>& psis,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          surface_(datesToTimes(referenceDate, dates, dc),
                   buildSlices(thetas, rhos, psis)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_()
    {
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    // =================================================================
    // Native params constructor with discrete dividends
    // =================================================================

    EssviVolatilityTermStructure::EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& thetas,
            const std::vector<Real>& rhos,
            const std::vector<Real>& psis,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          surface_(datesToTimes(referenceDate, dates, dc),
                   buildSlices(thetas, rhos, psis)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_(std::move(dividends))
    {
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    // =================================================================
    // Global params constructor
    // =================================================================

    EssviVolatilityTermStructure::EssviVolatilityTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<Real>& rhos,
            Real theta1,
            const std::vector<Real>& as,
            const std::vector<Real>& cs,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            EssviButterflyCondition::Type bflyType,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          surface_(datesToTimes(referenceDate, dates, dc),
                   [&]() {
                       EssviGlobalParams gp;
                       gp.rhos   = rhos;
                       gp.theta1 = theta1;
                       gp.as     = as;
                       gp.cs     = cs;
                       return gp;
                   }(),
                   bflyType),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_()
    {
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    // =================================================================
    // TermStructure interface
    // =================================================================

    Date EssviVolatilityTermStructure::maxDate() const {
        return Date::maxDate();
    }

    Real EssviVolatilityTermStructure::minStrike() const {
        return QL_MIN_REAL;
    }

    Real EssviVolatilityTermStructure::maxStrike() const {
        return QL_MAX_REAL;
    }

    // =================================================================
    // Inspectors
    // =================================================================

    Size EssviVolatilityTermStructure::numSlices() const {
        return surface_.numSlices();
    }

    const std::vector<EssviSliceParams>&
    EssviVolatilityTermStructure::slices() const {
        return surface_.slices();
    }

    // =================================================================
    // Visitability
    // =================================================================

    void EssviVolatilityTermStructure::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<EssviVolatilityTermStructure>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            BlackVolatilityTermStructure::accept(v);
    }

    // =================================================================
    // Gradient
    // =================================================================

    EssviSliceGradient
    EssviVolatilityTermStructure::impliedVolGradient(Size sliceIdx,
                                                      Real strike) const {
        Time t = surface_.maturities().at(sliceIdx);
        Real fwd = forward(t);
        Real k = std::log(strike / fwd);
        return surface_.impliedVolGradient(sliceIdx, k);
    }

    std::vector<Real>
    EssviVolatilityTermStructure::impliedVolGlobalGradient(
            Size sliceIdx, Real strike,
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {
        Time t = surface_.maturities().at(sliceIdx);
        Real fwd = forward(t);
        Real k = std::log(strike / fwd);
        return surface_.impliedVolGlobalGradient(sliceIdx, k, gp, bflyCond);
    }

    // =================================================================
    // Core implementation
    // =================================================================

    Real EssviVolatilityTermStructure::forward(Time t) const {
        Real S = spot_->value();
        Real df = riskFreeRate_->discount(t);
        Real dq = dividendYield_->discount(t);

        // Subtract PV of discrete dividends that ex before time t
        Real pvDivs = 0.0;
        if (!dividends_.empty()) {
            Date refDate = referenceDate();
            DayCounter dc = dayCounter();
            for (const auto& div : dividends_) {
                Date exDate = div->date();
                if (exDate <= refDate)
                    continue;
                Time tDiv = dc.yearFraction(refDate, exDate);
                if (tDiv < t) {
                    // PV of fixed dividend at risk-free rate
                    pvDivs += div->amount() * riskFreeRate_->discount(tDiv);
                }
            }
        }

        return (S - pvDivs) * dq / df;
    }

    Volatility EssviVolatilityTermStructure::blackVolImpl(Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;  // avoid division by zero
        Real fwd = forward(t);
        Real k = std::log(strike / fwd);
        return surface_.impliedVol(k, t);
    }

    // =================================================================
    // Batch evaluation
    // =================================================================

    std::vector<Real>
    EssviVolatilityTermStructure::batchBlackVol(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes) const {
        const Size nObs = sliceIndices.size();
        QL_REQUIRE(strikes.size() == nObs, "strikes size mismatch");

        // Precompute forwards per slice
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
    EssviVolatilityTermStructure::batchImpliedVolGlobalGradient(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes,
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {
        const Size nObs = sliceIndices.size();
        QL_REQUIRE(strikes.size() == nObs, "strikes size mismatch");
        const Size N = gp.numSlices();
        const Size nParams = 3 * N;

        // Precompute forwards
        const auto& T = surface_.maturities();
        std::vector<Real> fwds(N);
        for (Size s = 0; s < N; ++s)
            fwds[s] = forward(T[s]);

        // ── Factored Jacobian: compute chain Jacobian ONCE ──
        // chainJac layout: [dTheta(N×P), dPsi(N×P)] where P = nParams
        auto chainJac = surface_.chainJacobian(gp, bflyCond);
        const Real* dTheta = chainJac.data();
        const Real* dPsi   = chainJac.data() + N * nParams;

        // For each observation, combine per-slice gradient with chain Jacobian.
        // grad[j] = dSigma/dTheta * dTheta[si,j] + dSigma/dRho * delta(j==rhoOff+si)
        //         + dSigma/dPsi * dPsi[si,j]
        Size rhoOff = 0;

        std::vector<Real> J(nObs * nParams, 0.0);
        for (Size i = 0; i < nObs; ++i) {
            Size si = sliceIndices[i];
            Real k = std::log(strikes[i] / fwds[si]);
            EssviSliceGradient sg = essviImpliedVolGradient(
                surface_.slice(si), k, T[si]);

            Real* row = &J[i * nParams];
            const Real* dTh_row = &dTheta[si * nParams];
            const Real* dPs_row = &dPsi[si * nParams];

            for (Size j = 0; j < nParams; ++j) {
                row[j] = sg.dSigma_dTheta * dTh_row[j]
                       + sg.dSigma_dPsi   * dPs_row[j];
            }
            // Direct rho contribution (rho_si is a direct global param)
            row[rhoOff + si] += sg.dSigma_dRho;
        }
        return J;
    }

} // namespace QuantLib
