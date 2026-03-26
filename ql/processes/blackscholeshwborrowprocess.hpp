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

/*! \file blackscholeshwborrowprocess.hpp
    \brief Black-Scholes process with Hull-White stochastic rate
           and stochastic borrow/repo rate
*/

#ifndef quantlib_black_scholes_hw_borrow_process_hpp
#define quantlib_black_scholes_hw_borrow_process_hpp

#include <ql/stochasticprocess.hpp>
#include <ql/models/shortrate/onefactormodel.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/math/matrixutilities/choleskydecomposition.hpp>

namespace QuantLib {

    //! 3D stochastic process: equity + stochastic rate + stochastic borrow
    /*! The dynamics are:
        \f[
        \begin{aligned}
        dS &= (r - b) S \, dt + \sigma(t,S) S \, dW_S \\
        dr &= [\theta_r(t) - a_r r] \, dt + \sigma_r \, dW_r \\
        db &= [\theta_b(t) - a_b b] \, dt + \sigma_b \, dW_b
        \end{aligned}
        \f]

        The rate and borrow factors accept any OneFactorAffineModel
        (HullWhite, Vasicek, Gsr, etc.). Internally, the process
        extracts OU parameters and evolves its own states under the
        risk-neutral measure Q. It does NOT delegate to the model
        process's evolve() method, which may be in a different measure.

        \ingroup processes
    */
    class BlackScholesHWBorrowProcess : public StochasticProcess {
      public:
        BlackScholesHWBorrowProcess(
            Handle<Quote> spot,
            Handle<BlackVolTermStructure> volSurface,
            ext::shared_ptr<OneFactorAffineModel> rateModel,
            ext::shared_ptr<OneFactorAffineModel> borrowModel,
            Handle<Quote> corrEquityRate,
            Handle<Quote> corrEquityBorrow,
            Handle<Quote> corrRateBorrow);

        Size size() const override { return 3; }
        Size factors() const override { return 3; }

        Array initialValues() const override;
        Array drift(Time t, const Array& x) const override;
        Matrix diffusion(Time t, const Array& x) const override;
        Array evolve(Time t0, const Array& x0,
                     Time dt, const Array& dw) const override;

        void update() override;

        // --- Accessors ---
        const Handle<Quote>& spot() const { return spot_; }
        const Handle<BlackVolTermStructure>& volSurface() const {
            return volSurface_;
        }
        const ext::shared_ptr<OneFactorAffineModel>& rateModel() const {
            return rateModel_;
        }
        const ext::shared_ptr<OneFactorAffineModel>& borrowModel() const {
            return borrowModel_;
        }
        const Matrix& choleskyDecomposition() const;

        //! Deterministic shift alpha(t) for the rate/borrow factor.
        //! Computed model-agnostically via discountBond.
        Real rateAlpha(Time t) const;
        Real borrowAlpha(Time t) const;

        //! OU transition parameters (measure-independent).
        struct OUParams {
            Real decay;    //!< mean-reversion decay factor over dt
            Real volstep;  //!< conditional std dev of state over dt
        };
        OUParams rateOUParams(Time t, Time dt) const;
        OUParams borrowOUParams(Time t, Time dt) const;

        //! Equity diffusion coefficient at (t, S).
        Real equityVol(Time t, Real S) const;

      private:
        Handle<Quote> spot_;
        Handle<BlackVolTermStructure> volSurface_;
        ext::shared_ptr<OneFactorAffineModel> rateModel_;
        ext::shared_ptr<OneFactorAffineModel> borrowModel_;

        ext::shared_ptr<StochasticProcess1D> rateStateProcess_;
        ext::shared_ptr<StochasticProcess1D> borrowStateProcess_;
        Real rateX0_, borrowX0_;
        Real rateR0_, borrowR0_;  //!< initial short rates via dynamics()

        Handle<Quote> corrSR_, corrSB_, corrRB_;

        mutable Matrix choleskyDecomp_;
        mutable bool choleskyValid_;
        void rebuildCholesky() const;

        static constexpr Real eps_ = 1e-6;

        static ext::shared_ptr<StochasticProcess1D>
        extractProcess(const ext::shared_ptr<OneFactorAffineModel>& model);
    };

}

#endif
