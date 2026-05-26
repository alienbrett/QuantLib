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

/*! \file mixedspotlocalvoltermstructure.hpp
    \brief Analytic local-vol surface for the spot-mixture decorator

    Density-weighted local vol for a price-mixture C_mix(K,T) =
    Σᵢ wᵢ Cᵢ(Fᵢ, K, σᵢ) where Fᵢ = F · exp(σ_F · x_i) / cosh(σ_F)
    and σᵢ = base.blackVol(T, K · F/Fᵢ).

    Identity (Dupire on a mixture; see Gatheral, "The Volatility
    Surface", §5):

        σ_LV,mix²(K, T) = E[σ_LV²(K, T) | S_T = K]
                       = (Σᵢ wᵢ ρᵢ(K) · σ_LV,base²(K · F/Fᵢ, T))
                         / (Σᵢ wᵢ ρᵢ(K))

    where ρᵢ(K) is the log-normal density at K under branch i with
    forward Fᵢ and BS vol σᵢ.

    Per-query cost: ``n_quad`` base-LV queries + ``n_quad`` base-Black
    queries + ``n_quad`` log-normal density evaluations + simple
    arithmetic.  Compared to numerical Dupire on the same mixture
    (NoExceptLocalVolSurface on MixedSpotVolTermStructure):
    ~10× faster.  Compared to base analytic LV alone (e.g.,
    ParametricLocalVolSurface): ~10× slower — the cost of the
    bimodality overlay.
*/

#ifndef quantlib_mixed_spot_local_vol_term_structure_hpp
#define quantlib_mixed_spot_local_vol_term_structure_hpp

#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/localvoltermstructure.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/quote.hpp>
#include <ql/handle.hpp>
#include <vector>

namespace QuantLib {

    //! Analytic local-vol surface for the spot-mixture overlay
    /*! Companion to MixedSpotVolTermStructure (the Black surface).
        Wraps a base LocalVolTermStructure (analytic Dupire of the
        base parametric Black surface, e.g. ParametricLocalVolSurface
        or EssviLocalVolSurface) and the base Black surface
        (needed for branch log-normal densities).

        See the file comment for the math.
    */
    class MixedSpotLocalVolTermStructure : public LocalVolTermStructure {
      public:
        //! Construct with default symmetric 2-point quadrature
        MixedSpotLocalVolTermStructure(
            Handle<LocalVolTermStructure> baseLV,
            Handle<BlackVolTermStructure> baseBlack,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Real sigmaF,
            Real fwdDecayTau = 0.0);

        //! Construct with custom quadrature nodes/weights
        MixedSpotLocalVolTermStructure(
            Handle<LocalVolTermStructure> baseLV,
            Handle<BlackVolTermStructure> baseBlack,
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
        //@}

        //! \name Mutators
        //@{
        void setSigmaF(Real sigmaF);
        void setFwdDecayTau(Real tau);
        //@}

        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&) override;
        //@}

      protected:
        Volatility localVolImpl(Time t, Real strike) const override;

      private:
        Real forward_(Time t) const;
        Real sigmaFAt_(Time t) const;
        void validate_() const;

        Handle<LocalVolTermStructure> baseLV_;
        Handle<BlackVolTermStructure> baseBlack_;
        Handle<Quote> spot_;
        Handle<YieldTermStructure> riskFreeRate_;
        Handle<YieldTermStructure> dividendYield_;
        Real sigmaF_;
        Real fwdDecayTau_;
        std::vector<Real> nodes_;
        std::vector<Real> weights_;
    };

}

#endif
