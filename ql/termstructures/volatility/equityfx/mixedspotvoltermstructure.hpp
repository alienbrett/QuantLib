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

/*! \file mixedspotvoltermstructure.hpp
    \brief Spot-mixture decorator over any Black volatility surface

    Implements the "random-coefficient" overlay from the mixture eSSVI
    spec (docs/mixture_essvi_spec.md), restricted to the single global
    forward-dispersion channel (sigma_F). Generic over any
    BlackVolTermStructure base: parametric (JW, eSSVI, SABR, AH, PWL)
    or interpolated.

    The construction is:

        F_i(T) = F(T) * exp(sigma_F * x_i) / cosh(sigma_F)
        C_mix(T,K) = sum_i w_i * C(T, K; F_i, vol_i)
        sigma_mix(T,K) = BS_inverse(C_mix)

    where vol_i = base.blackVol(T, K * F/F_i) — i.e. the base surface
    is queried at the synthetic strike that makes its internal
    log-moneyness equal to log(K/F_i). This avoids any need for the
    base to expose its forward or to be perturbed itself.

    Mix prices, never vols (spec section 1, 9). Default quadrature is
    the symmetric 2-point ±1 with weights 0.5/0.5.

    Optional exponential decay of sigma_F with maturity:

        sigma_F(T) = sigma_F * exp(-T / decay_tau)

    Useful for event-driven dispersion that fades at long tenors.
*/

#ifndef quantlib_mixed_spot_vol_term_structure_hpp
#define quantlib_mixed_spot_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <vector>

namespace QuantLib {

    //! Spot-mixture decorator over a base Black vol term structure
    /*! Generic over any BlackVolTermStructure. Builds a bimodal smile
        from a unimodal base via a 2-point spot-forward mixture with
        a single global parameter sigma_F.

        See class docstring (file comment) for the math.
    */
    class MixedSpotVolTermStructure : public BlackVolatilityTermStructure {
      public:
        //! Construct with default symmetric 2-point quadrature
        MixedSpotVolTermStructure(
            Handle<BlackVolTermStructure> base,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau = 0.0);

        //! Construct with custom quadrature nodes/weights
        MixedSpotVolTermStructure(
            Handle<BlackVolTermStructure> base,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau,
            std::vector<Real> quadNodes,
            std::vector<Real> quadWeights);

        //! \name TermStructure interface
        //@{
        Date maxDate() const override;
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Real minStrike() const override;
        Real maxStrike() const override;
        //@}

        //! \name Inspectors
        //@{
        Real sigmaF() const { return sigmaF_; }
        Real fwdDecayTau() const { return fwdDecayTau_; }
        const std::vector<Real>& quadNodes() const { return nodes_; }
        const std::vector<Real>& quadWeights() const { return weights_; }
        //@}

        //! \name Mutators
        //@{
        //! Update sigma_F in place; notifies observers
        void setSigmaF(Real sigmaF);
        //! Update fwd decay tau in place; notifies observers
        void setFwdDecayTau(Real tau);
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility blackVolImpl(Time t, Real strike) const override;

      private:
        Real forward_(Time t) const;
        Real sigmaFAt_(Time t) const;
        void validate_() const;

        Handle<BlackVolTermStructure> base_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        Real sigmaF_;
        Real fwdDecayTau_;
        std::vector<Real> nodes_;
        std::vector<Real> weights_;
    };

} // namespace QuantLib

#endif
