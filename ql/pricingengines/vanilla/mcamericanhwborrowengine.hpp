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

/*! \file mcamericanhwborrowengine.hpp
    \brief Monte Carlo American option engine with stochastic rates
           and stochastic borrow (3-factor: S, r, b)
*/

#ifndef quantlib_mc_american_hw_borrow_engine_hpp
#define quantlib_mc_american_hw_borrow_engine_hpp

#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>
#include <ql/processes/blackscholeshwborrowprocess.hpp>
#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/math/matrixutilities/svd.hpp>
#include <vector>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    //! Monte Carlo American engine with stochastic rate and borrow
    /*! Uses Longstaff-Schwartz regression on (S, r, b) for early
        exercise decisions. Path generation uses exact OU transitions
        under the risk-neutral measure Q.

        \ingroup vanillaengines
    */
    template <class RNG = PseudoRandom>
    class MCAmericanHWBorrowEngine
        : public GenericEngine<VanillaOption::arguments,
                               VanillaOption::results> {
      public:
        MCAmericanHWBorrowEngine(
            ext::shared_ptr<BlackScholesHWBorrowProcess> process,
            Size timeSteps,
            Size requiredSamples,
            BigInteger seed = 0,
            bool antitheticVariate = false);

        void calculate() const override;

      private:
        struct TimeStepParams {
            Real dt, sqrt_dt;
            Real alpha_r, alpha_b;
            Real decay_r, volstep_r;
            Real decay_b, volstep_b;
        };

        std::vector<TimeStepParams> buildTimeGrid(Time T) const;

        void generatePaths(
            const std::vector<TimeStepParams>& grid,
            std::vector<std::vector<Real>>& S,
            std::vector<std::vector<Real>>& r,
            std::vector<std::vector<Real>>& b,
            std::vector<std::vector<Real>>& cumR) const;

        Real longstaffSchwartz(
            const std::vector<std::vector<Real>>& S,
            const std::vector<std::vector<Real>>& r,
            const std::vector<std::vector<Real>>& b,
            const std::vector<std::vector<Real>>& cumR,
            Real strike, Option::Type type) const;

        ext::shared_ptr<BlackScholesHWBorrowProcess> process_;
        Size timeSteps_;
        Size requiredSamples_;
        BigInteger seed_;
        bool antitheticVariate_;
    };

    typedef MCAmericanHWBorrowEngine<PseudoRandom>
        MCPRAmericanHWBorrowEngine;
    typedef MCAmericanHWBorrowEngine<LowDiscrepancy>
        MCLDAmericanHWBorrowEngine;


    // ====================================================================
    // Template implementation
    // ====================================================================

    template <class RNG>
    MCAmericanHWBorrowEngine<RNG>::MCAmericanHWBorrowEngine(
        ext::shared_ptr<BlackScholesHWBorrowProcess> process,
        Size timeSteps,
        Size requiredSamples,
        BigInteger seed,
        bool antitheticVariate)
        : process_(std::move(process)),
          timeSteps_(timeSteps),
          requiredSamples_(requiredSamples),
          seed_(seed),
          antitheticVariate_(antitheticVariate)
    {
        registerWith(process_);
    }


    template <class RNG>
    std::vector<typename MCAmericanHWBorrowEngine<RNG>::TimeStepParams>
    MCAmericanHWBorrowEngine<RNG>::buildTimeGrid(Time T) const
    {
        Size m = timeSteps_;
        Time dt = T / m;
        std::vector<TimeStepParams> grid(m + 1);

        for (Size i = 0; i <= m; ++i) {
            Time t = i * dt;
            grid[i].dt = dt;
            grid[i].sqrt_dt = std::sqrt(dt);
            grid[i].alpha_r = process_->rateAlpha(t);
            grid[i].alpha_b = process_->borrowAlpha(t);

            if (i < m) {
                auto rp = process_->rateOUParams(t, dt);
                grid[i].decay_r = rp.decay;
                grid[i].volstep_r = rp.volstep;
                auto bp = process_->borrowOUParams(t, dt);
                grid[i].decay_b = bp.decay;
                grid[i].volstep_b = bp.volstep;
            }
        }
        return grid;
    }


    template <class RNG>
    void MCAmericanHWBorrowEngine<RNG>::generatePaths(
        const std::vector<TimeStepParams>& grid,
        std::vector<std::vector<Real>>& S,
        std::vector<std::vector<Real>>& r,
        std::vector<std::vector<Real>>& b,
        std::vector<std::vector<Real>>& cumR) const
    {
        Size n = requiredSamples_;
        Size m = timeSteps_;
        Real S0 = process_->spot()->value();

        const Matrix& L = process_->choleskyDecomposition();
        Real L00 = L[0][0];
        Real L10 = L[1][0], L11 = L[1][1];
        Real L20 = L[2][0], L21 = L[2][1], L22 = L[2][2];

        typename RNG::rsg_type rsg =
            RNG::make_sequence_generator(3 * m, seed_);

        for (Size path = 0; path < n; ++path) {
            const auto& seq = rsg.nextSequence().value;

            Real dx_r = 0.0, dx_b = 0.0;
            Real logS = std::log(S0);
            Real cum_r = 0.0;

            S[path][0] = S0;
            r[path][0] = grid[0].alpha_r;
            b[path][0] = grid[0].alpha_b;
            cumR[path][0] = 0.0;

            for (Size i = 0; i < m; ++i) {
                const auto& p = grid[i];

                Real z0 = seq[3*i], z1 = seq[3*i+1], z2 = seq[3*i+2];
                Real w_S = L00*z0;
                Real w_r = L10*z0 + L11*z1;
                Real w_b = L20*z0 + L21*z1 + L22*z2;

                Real r_t = p.alpha_r + dx_r;
                Real b_t = p.alpha_b + dx_b;

                cum_r += r_t * p.dt;

                // Exact OU evolution of deviations under Q
                dx_r = p.decay_r * dx_r + p.volstep_r * w_r;
                dx_b = p.decay_b * dx_b + p.volstep_b * w_b;

                // Equity GBM step
                Real sigma = process_->equityVol(i * p.dt, std::exp(logS));
                logS += (r_t - b_t - 0.5*sigma*sigma) * p.dt
                      + sigma * p.sqrt_dt * w_S;

                S[path][i+1] = std::exp(logS);
                r[path][i+1] = grid[i+1].alpha_r + dx_r;
                b[path][i+1] = grid[i+1].alpha_b + dx_b;
                cumR[path][i+1] = cum_r;
            }
        }
    }


    template <class RNG>
    Real MCAmericanHWBorrowEngine<RNG>::longstaffSchwartz(
        const std::vector<std::vector<Real>>& S,
        const std::vector<std::vector<Real>>& r,
        const std::vector<std::vector<Real>>& b,
        const std::vector<std::vector<Real>>& cumR,
        Real strike, Option::Type type) const
    {
        Size n = requiredSamples_;
        Size m = timeSteps_;

        // Cashflow[path] = discounted payoff at the optimal stopping time
        std::vector<Real> cashflow(n, 0.0);
        std::vector<Size> exerciseTime(n, m);

        // Terminal payoff
        for (Size path = 0; path < n; ++path) {
            Real payoff = (type == Option::Call)
                ? std::max(S[path][m] - strike, 0.0)
                : std::max(strike - S[path][m], 0.0);
            cashflow[path] = payoff * std::exp(-cumR[path][m]);
        }

        // Backward induction
        for (Size step = m - 1; step >= 1; --step) {
            // Identify in-the-money paths at this step
            std::vector<Size> itm;
            itm.reserve(n);
            for (Size path = 0; path < n; ++path) {
                Real payoff = (type == Option::Call)
                    ? std::max(S[path][step] - strike, 0.0)
                    : std::max(strike - S[path][step], 0.0);
                if (payoff > 0.0)
                    itm.push_back(path);
            }

            if (itm.empty())
                continue;

            // Build regression matrix: polynomial basis on (S, r, b)
            // Basis: {1, s, r, b, s^2, sr, sb, r^2, rb, b^2}
            // where s = S/strike (normalized)
            Size nBasis = 10;
            Size nItm = itm.size();

            // Simple OLS via normal equations
            // A^T A beta = A^T y
            Matrix A(nItm, nBasis);
            Array y(nItm);

            for (Size j = 0; j < nItm; ++j) {
                Size path = itm[j];
                Real s = S[path][step] / strike;
                Real rv = r[path][step];
                Real bv = b[path][step];

                A[j][0] = 1.0;
                A[j][1] = s;
                A[j][2] = rv;
                A[j][3] = bv;
                A[j][4] = s * s;
                A[j][5] = s * rv;
                A[j][6] = s * bv;
                A[j][7] = rv * rv;
                A[j][8] = rv * bv;
                A[j][9] = bv * bv;

                // Continuation value: discounted future cashflow
                // relative to this step's discount factor
                Real discFactor = std::exp(-(cumR[path][exerciseTime[path]]
                                           - cumR[path][step]));
                Real futureCF = (type == Option::Call)
                    ? std::max(S[path][exerciseTime[path]] - strike, 0.0)
                    : std::max(strike - S[path][exerciseTime[path]], 0.0);
                y[j] = futureCF * discFactor;
            }

            // Solve via SVD for stability
            SVD svd(A);
            Array beta = svd.solveFor(y);

            // Exercise decision
            for (Size j = 0; j < nItm; ++j) {
                Size path = itm[j];
                Real s = S[path][step] / strike;
                Real rv = r[path][step];
                Real bv = b[path][step];

                Real continuation = beta[0]
                    + beta[1]*s + beta[2]*rv + beta[3]*bv
                    + beta[4]*s*s + beta[5]*s*rv + beta[6]*s*bv
                    + beta[7]*rv*rv + beta[8]*rv*bv + beta[9]*bv*bv;

                Real exerciseVal = (type == Option::Call)
                    ? std::max(S[path][step] - strike, 0.0)
                    : std::max(strike - S[path][step], 0.0);

                if (exerciseVal >= continuation) {
                    exerciseTime[path] = step;
                    cashflow[path] = exerciseVal
                                   * std::exp(-cumR[path][step]);
                }
            }
        }

        // Average discounted cashflows
        Real sum = 0.0;
        for (Size path = 0; path < n; ++path)
            sum += cashflow[path];

        return sum / n;
    }


    template <class RNG>
    void MCAmericanHWBorrowEngine<RNG>::calculate() const
    {
        ext::shared_ptr<StrikedTypePayoff> payoff =
            ext::dynamic_pointer_cast<StrikedTypePayoff>(
                this->arguments_.payoff);
        QL_REQUIRE(payoff, "StrikedTypePayoff required");

        ext::shared_ptr<Exercise> exercise = this->arguments_.exercise;
        QL_REQUIRE(exercise->type() == Exercise::American,
                   "American exercise required");

        Date today = Settings::instance().evaluationDate();
        Date maturity = exercise->lastDate();
        DayCounter dc = process_->volSurface()->dayCounter();
        Time T = dc.yearFraction(today, maturity);

        QL_REQUIRE(T > 0.0, "option has expired");

        auto grid = buildTimeGrid(T);

        Size n = requiredSamples_;
        Size m = timeSteps_;

        // Allocate path storage
        std::vector<std::vector<Real>> S(n, std::vector<Real>(m+1));
        std::vector<std::vector<Real>> r(n, std::vector<Real>(m+1));
        std::vector<std::vector<Real>> b(n, std::vector<Real>(m+1));
        std::vector<std::vector<Real>> cumR(n, std::vector<Real>(m+1));

        generatePaths(grid, S, r, b, cumR);

        Real npv = longstaffSchwartz(
            S, r, b, cumR,
            payoff->strike(), payoff->optionType());

        this->results_.value = std::max(0.0, npv);
    }

}

#endif
