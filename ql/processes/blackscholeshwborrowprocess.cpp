/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2025 Brett Graves

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

#include <ql/processes/blackscholeshwborrowprocess.hpp>
#include <utility>

namespace QuantLib {

    ext::shared_ptr<StochasticProcess1D>
    BlackScholesHWBorrowProcess::extractProcess(
        const ext::shared_ptr<OneFactorAffineModel>& model)
    {
        auto ofm = ext::dynamic_pointer_cast<OneFactorModel>(model);
        QL_REQUIRE(ofm, "model must be a OneFactorModel");
        return ofm->dynamics()->process();
    }

    BlackScholesHWBorrowProcess::BlackScholesHWBorrowProcess(
        Handle<Quote> spot,
        Handle<BlackVolTermStructure> volSurface,
        ext::shared_ptr<OneFactorAffineModel> rateModel,
        ext::shared_ptr<OneFactorAffineModel> borrowModel,
        Handle<Quote> corrSR,
        Handle<Quote> corrSB,
        Handle<Quote> corrRB)
        : spot_(std::move(spot)),
          volSurface_(std::move(volSurface)),
          rateModel_(std::move(rateModel)),
          borrowModel_(std::move(borrowModel)),
          corrSR_(std::move(corrSR)),
          corrSB_(std::move(corrSB)),
          corrRB_(std::move(corrRB)),
          choleskyValid_(false)
    {
        rateStateProcess_ = extractProcess(rateModel_);
        borrowStateProcess_ = extractProcess(borrowModel_);
        rateX0_ = rateStateProcess_->x0();
        borrowX0_ = borrowStateProcess_->x0();

        // Get initial short rates via the model dynamics mapping
        auto rateOFM = ext::dynamic_pointer_cast<OneFactorModel>(rateModel_);
        auto borrowOFM = ext::dynamic_pointer_cast<OneFactorModel>(borrowModel_);
        rateR0_ = rateOFM->dynamics()->shortRate(0.0, rateX0_);
        borrowR0_ = borrowOFM->dynamics()->shortRate(0.0, borrowX0_);

        registerWith(spot_);
        registerWith(volSurface_);
        registerWith(corrSR_);
        registerWith(corrSB_);
        registerWith(corrRB_);
    }

    void BlackScholesHWBorrowProcess::update() {
        choleskyValid_ = false;
        StochasticProcess::update();
    }

    Array BlackScholesHWBorrowProcess::initialValues() const {
        Real S0 = spot_->value();
        Real r0 = rateAlpha(0.0);
        Real b0 = borrowAlpha(0.0);
        Array v(3);
        v[0] = S0;
        v[1] = r0;
        v[2] = b0;
        return v;
    }

    Array BlackScholesHWBorrowProcess::drift(Time t, const Array& x) const {
        Real S = x[0];
        Real r = x[1];
        Real b = x[2];
        Real sigma = equityVol(t, S);

        Array d(3);
        d[0] = (r - b) * S;
        // For the rate/borrow drift, we'd need theta(t) - a*x.
        // This is approximate; the engine uses exact OU transitions instead.
        d[1] = 0.0;
        d[2] = 0.0;
        return d;
    }

    Matrix BlackScholesHWBorrowProcess::diffusion(Time t,
                                                   const Array& x) const {
        Real S = x[0];
        Real sigma = equityVol(t, S);

        auto rp = rateOUParams(t, eps_);
        auto bp = borrowOUParams(t, eps_);
        Real sigR = rp.volstep / std::sqrt(eps_);
        Real sigB = bp.volstep / std::sqrt(eps_);

        const Matrix& L = choleskyDecomposition();

        Matrix diff(3, 3, 0.0);
        diff[0][0] = sigma * S * L[0][0];
        diff[0][1] = sigma * S * L[0][1];
        diff[0][2] = sigma * S * L[0][2];
        diff[1][0] = sigR * L[1][0];
        diff[1][1] = sigR * L[1][1];
        diff[1][2] = sigR * L[1][2];
        diff[2][0] = sigB * L[2][0];
        diff[2][1] = sigB * L[2][1];
        diff[2][2] = sigB * L[2][2];
        return diff;
    }

    Array BlackScholesHWBorrowProcess::evolve(Time t0, const Array& x0,
                                               Time dt,
                                               const Array& dw) const {
        const Matrix& L = choleskyDecomposition();

        // Correlate the independent normals
        Real w_S = L[0][0]*dw[0];
        Real w_r = L[1][0]*dw[0] + L[1][1]*dw[1];
        Real w_b = L[2][0]*dw[0] + L[2][1]*dw[1] + L[2][2]*dw[2];

        Real S = x0[0];
        Real r_t = x0[1];
        Real b_t = x0[2];

        // OU state deviation from initial
        Real dx_r = r_t - rateAlpha(t0);
        Real dx_b = b_t - borrowAlpha(t0);

        // Exact OU transition (under Q)
        auto rp = rateOUParams(t0, dt);
        auto bp = borrowOUParams(t0, dt);
        Real new_dx_r = rp.decay * dx_r + rp.volstep * w_r;
        Real new_dx_b = bp.decay * dx_b + bp.volstep * w_b;

        Real sigma = equityVol(t0, S);
        Real sqrt_dt = std::sqrt(dt);
        Real logS = std::log(S)
                   + (r_t - b_t - 0.5*sigma*sigma) * dt
                   + sigma * sqrt_dt * w_S;

        Array x1(3);
        x1[0] = std::exp(logS);
        x1[1] = rateAlpha(t0 + dt) + new_dx_r;
        x1[2] = borrowAlpha(t0 + dt) + new_dx_b;
        return x1;
    }

    const Matrix&
    BlackScholesHWBorrowProcess::choleskyDecomposition() const {
        if (!choleskyValid_)
            rebuildCholesky();
        return choleskyDecomp_;
    }

    void BlackScholesHWBorrowProcess::rebuildCholesky() const {
        Real rhoSR = corrSR_->value();
        Real rhoSB = corrSB_->value();
        Real rhoRB = corrRB_->value();

        Matrix corrMatrix(3, 3, 0.0);
        corrMatrix[0][0] = corrMatrix[1][1] = corrMatrix[2][2] = 1.0;
        corrMatrix[0][1] = corrMatrix[1][0] = rhoSR;
        corrMatrix[0][2] = corrMatrix[2][0] = rhoSB;
        corrMatrix[1][2] = corrMatrix[2][1] = rhoRB;
        choleskyDecomp_ = CholeskyDecomposition(corrMatrix);
        choleskyValid_ = true;
    }

    Real BlackScholesHWBorrowProcess::rateAlpha(Time t) const {
        // Forward rate from the initial curve at time t, extracted
        // model-agnostically via discountBond(0, T, r0).
        Time t_safe = std::max(t, eps_);
        return -(std::log(rateModel_->discountBond(0, t_safe + eps_, rateR0_))
               - std::log(rateModel_->discountBond(0, t_safe, rateR0_)))
               / eps_;
    }

    Real BlackScholesHWBorrowProcess::borrowAlpha(Time t) const {
        Time t_safe = std::max(t, eps_);
        return -(std::log(borrowModel_->discountBond(0, t_safe + eps_, borrowR0_))
               - std::log(borrowModel_->discountBond(0, t_safe, borrowR0_)))
               / eps_;
    }

    BlackScholesHWBorrowProcess::OUParams
    BlackScholesHWBorrowProcess::rateOUParams(Time t, Time dt) const {
        OUParams p;
        p.decay = rateStateProcess_->expectation(t, rateX0_ + 1.0, dt)
                - rateStateProcess_->expectation(t, rateX0_, dt);
        p.volstep = rateStateProcess_->stdDeviation(t, rateX0_, dt);
        return p;
    }

    BlackScholesHWBorrowProcess::OUParams
    BlackScholesHWBorrowProcess::borrowOUParams(Time t, Time dt) const {
        OUParams p;
        p.decay = borrowStateProcess_->expectation(t, borrowX0_ + 1.0, dt)
                - borrowStateProcess_->expectation(t, borrowX0_, dt);
        p.volstep = borrowStateProcess_->stdDeviation(t, borrowX0_, dt);
        return p;
    }

    Real BlackScholesHWBorrowProcess::equityVol(Time t, Real S) const {
        return volSurface_->blackVol(t, S);
    }

}
