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

/*! \file essvihelpers.hpp
    \brief Global arbitrage-free eSSVI volatility surface helpers.

    Implements the parametrization from:
      Mingone (2022), "No arbitrage global parametrization for the eSSVI
      volatility surface", arXiv:2204.00312v1.

    The surface is parameterized per-slice by (theta_i, rho_i, psi_i) where:
      theta = ATM total implied variance
      rho   = correlation (skew), in (-1, 1)
      psi   = theta * phi, controls curvature
*/

#ifndef quantlib_essvi_helpers_hpp
#define quantlib_essvi_helpers_hpp

#include <ql/types.hpp>
#include <ql/errors.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace QuantLib {

    //! Native SSVI triple for a single expiry slice
    struct EssviSliceParams {
        Real theta;  //!< ATM total implied variance, > 0
        Real rho;    //!< correlation / skew, in (-1, 1)
        Real psi;    //!< theta * phi, curvature parameter, > 0

        Real phi() const { return theta > 0.0 ? psi / theta : 0.0; }
    };

    //! Global (arb-free) parameter set for the eSSVI surface
    struct EssviGlobalParams {
        std::vector<Real> rhos;    //!< rho_1 .. rho_N, each in (-1,1)
        Real              theta1;  //!< theta for the first slice, > 0
        std::vector<Real> as;      //!< a_2 .. a_N, each > 0
        std::vector<Real> cs;      //!< c_1 .. c_N, each in (0,1)

        Size numSlices() const { return rhos.size(); }
    };

    //! Butterfly arbitrage bound selector
    struct EssviButterflyCondition {
        enum Type {
            GatheralJacquier,   //!< sufficient (GJ): psi^2 <= 4*theta/(1+|rho|)
            MartiniMingone      //!< necessary & sufficient (MM): requires minimization
        };
    };

    // ---------------------------------------------------------------
    // Free functions: arbitrage bound helpers
    // ---------------------------------------------------------------

    //! GJ butterfly upper bound: f_GJ(theta, |rho|) = 4*theta / (1+|rho|)
    inline Real essviButterflyBoundGJ(Real theta, Real absRho) {
        return 4.0 * theta / (1.0 + absRho);
    }

    //! Roger Lee moment-formula bound on psi
    inline Real essviPsiUpperLee(Real absRho) {
        return 4.0 / (1.0 + absRho);
    }

    //! p_i = max((1+rho_{i-1})/(1+rho_i), (1-rho_{i-1})/(1-rho_i))
    inline Real essviCalendarP(Real rhoPrev, Real rhoCurr) {
        Real a = (1.0 + rhoPrev) / (1.0 + rhoCurr);
        Real b = (1.0 - rhoPrev) / (1.0 - rhoCurr);
        return std::max(a, b);
    }

    // ---------------------------------------------------------------
    // Martini-Mingone (MM) auxiliary functions
    // ---------------------------------------------------------------

    //! N(l, r) = sqrt(1 - r^2) + r*l + sqrt(l^2 + 1)
    inline Real essviMmN(Real l, Real r) {
        return std::sqrt(1.0 - r * r) + r * l + std::sqrt(l * l + 1.0);
    }

    //! N'(l, r) = r + l / sqrt(l^2 + 1)
    inline Real essviMmNprime(Real l, Real r) {
        return r + l / std::sqrt(l * l + 1.0);
    }

    //! N''(l) = 1 / (l^2 + 1)^{3/2}
    inline Real essviMmNpp(Real l) {
        Real u = l * l + 1.0;
        return 1.0 / (u * std::sqrt(u));
    }

    //! g(l, r) = N'(l, r) / 4
    inline Real essviMmG(Real l, Real r) {
        return essviMmNprime(l, r) / 4.0;
    }

    //! h(l, r) = 1 - (l - r / sqrt(1 - r^2)) * N'(l, r) / (2 * N(l, r))
    inline Real essviMmH(Real l, Real r) {
        Real s = std::sqrt(1.0 - r * r);
        Real Np = essviMmNprime(l, r);
        Real N  = essviMmN(l, r);
        return 1.0 - (l - r / s) * Np / (2.0 * N);
    }

    //! g2(l, r) = N''(l) - N'(l, r)^2 / (2 * N(l, r))
    inline Real essviMmG2(Real l, Real r) {
        Real Np  = essviMmNprime(l, r);
        Real Npp = essviMmNpp(l);
        Real N   = essviMmN(l, r);
        return Npp - Np * Np / (2.0 * N);
    }

    //! l2(|rho|) = tan(arccos(-|rho|) / 3)  (Fukasawa critical point)
    inline Real essviMmL2(Real absRho) {
        return std::tan(std::acos(-absRho) / 3.0);
    }

    //! R(l, theta, |rho|) ratio whose infimum gives f_MM
    inline Real essviMmRatio(Real l, Real theta, Real absRho) {
        Real s2   = 1.0 - absRho * absRho;
        Real srt  = std::sqrt(s2);
        Real g    = essviMmG(l, absRho);
        Real h    = essviMmH(l, absRho);
        Real g2   = essviMmG2(l, absRho);
        Real num  = 4.0 * theta * srt * h * h;
        Real den  = theta * srt * g * g - g2;
        if (den <= 0.0) return std::numeric_limits<Real>::max();
        return num / den;
    }

    //! MM butterfly bound: f_MM(theta, |rho|) = inf_{l > l2} R(l, theta, |rho|)
    Real essviButterflyBoundMM(Real theta, Real absRho);

    //! Unified butterfly bound: min(Lee bound, sqrt(f(theta, |rho|)))
    Real essviSliceButterflyBound(Real theta, Real absRho,
                                  EssviButterflyCondition::Type cond);

    // ---------------------------------------------------------------
    // Single-slice eSSVI total variance formula
    //   w(k) = 0.5 * (theta + rho*psi*k
    //           + sqrt((psi*k + theta*rho)^2 + theta^2*(1-rho^2)))
    // ---------------------------------------------------------------

    inline Real essviTotalVariance(const EssviSliceParams& s, Real k) {
        Real inner = s.psi * k + s.theta * s.rho;
        return 0.5 * (s.theta + s.rho * s.psi * k
                       + std::sqrt(inner * inner
                                   + s.theta * s.theta * (1.0 - s.rho * s.rho)));
    }

    // ---------------------------------------------------------------
    // Analytic gradient of eSSVI implied vol w.r.t. native slice params
    // ---------------------------------------------------------------

    //! Gradient of implied vol w.r.t. (theta, rho, psi) for one slice
    struct EssviSliceGradient {
        Real dSigma_dTheta;
        Real dSigma_dRho;
        Real dSigma_dPsi;
    };

    /*! Analytic gradient of total variance w(k) w.r.t. (theta, rho, psi).
        Let A = psi*k + theta*rho, D = A^2 + theta^2*(1-rho^2), sqrtD = sqrt(D).
        Then:
          dw/dtheta = 0.5 * (1 + (A*rho + theta*(1-rho^2)) / sqrtD)
          dw/drho   = 0.5 * psi*k * (1 + theta / sqrtD)
          dw/dpsi   = 0.5 * k * (rho + A / sqrtD)
    */
    inline void essviTotalVarianceGradient(const EssviSliceParams& s, Real k,
                                           Real& dw_dTheta,
                                           Real& dw_dRho,
                                           Real& dw_dPsi) {
        Real A = s.psi * k + s.theta * s.rho;
        Real D = A * A + s.theta * s.theta * (1.0 - s.rho * s.rho);
        Real sqrtD = std::sqrt(D);

        dw_dTheta = 0.5 * (1.0 + (A * s.rho + s.theta * (1.0 - s.rho * s.rho)) / sqrtD);
        dw_dRho   = 0.5 * s.psi * k * (1.0 + s.theta / sqrtD);
        dw_dPsi   = 0.5 * k * (s.rho + A / sqrtD);
    }

    /*! Analytic gradient of implied vol sigma = sqrt(w/T) w.r.t. (theta, rho, psi).
        dsigma/dparam = (dw/dparam) / (2 * T * sigma)
    */
    inline EssviSliceGradient essviImpliedVolGradient(
            const EssviSliceParams& s, Real k, Real T) {
        Real w = essviTotalVariance(s, k);
        Real sigma = std::sqrt(w / T);
        Real denom = 2.0 * T * sigma;

        Real dw_dTheta, dw_dRho, dw_dPsi;
        essviTotalVarianceGradient(s, k, dw_dTheta, dw_dRho, dw_dPsi);

        return { dw_dTheta / denom,
                 dw_dRho   / denom,
                 dw_dPsi   / denom };
    }

    /*! Analytic first and second derivatives of total variance w.r.t. log-moneyness k.
        dw/dk   = 0.5 * (rho*psi + psi*A/sqrt(D))
        d2w/dk2 = 0.5 * psi^2 * theta^2 * (1-rho^2) / D^{3/2}
        where A = psi*k + theta*rho, D = A^2 + theta^2*(1-rho^2).
    */
    inline void essviTotalVarianceStrikeDerivatives(
            const EssviSliceParams& s, Real k,
            Real& dwdk, Real& d2wdk2) {
        Real A = s.psi * k + s.theta * s.rho;
        Real D = A * A + s.theta * s.theta * (1.0 - s.rho * s.rho);
        Real sqrtD = std::sqrt(D);

        dwdk   = 0.5 * (s.rho * s.psi + s.psi * A / sqrtD);
        d2wdk2 = 0.5 * s.psi * s.psi * s.theta * s.theta
                 * (1.0 - s.rho * s.rho) / (D * sqrtD);
    }

    // ---------------------------------------------------------------
    // eSSVI core surface (standalone, no QuantLib term structure)
    // ---------------------------------------------------------------

    class EssviSurface {
    public:
        //! Construct from global (arb-free) parameters and maturities.
        EssviSurface(const std::vector<Real>& maturities,
                     const EssviGlobalParams& gparams,
                     EssviButterflyCondition::Type bflyCond
                         = EssviButterflyCondition::GatheralJacquier);

        //! Construct directly from per-slice native parameters.
        EssviSurface(const std::vector<Real>& maturities,
                     const std::vector<EssviSliceParams>& slices);

        //! Total implied variance w(k, T_i) for slice index i
        Real totalVariance(Size sliceIdx, Real logMoneyness) const;

        //! Total implied variance with maturity interpolation / extrapolation
        Real totalVariance(Real logMoneyness, Real maturity) const;

        //! Implied volatility sigma(k, T)
        Real impliedVol(Real logMoneyness, Real maturity) const;

        Size numSlices() const { return slices_.size(); }
        const std::vector<Real>& maturities() const { return T_; }
        const std::vector<EssviSliceParams>& slices() const { return slices_; }
        const EssviSliceParams& slice(Size i) const { return slices_.at(i); }

        //! Strike derivatives of total variance (analytic)
        void totalVarianceStrikeDerivatives(Real logMoneyness, Real maturity,
                                            Real& dwdk, Real& d2wdk2) const;

        //! Time derivative of total variance at (k, T) using param interpolation
        Real totalVarianceTimeDerivative(Real logMoneyness, Real maturity) const;

        //! Analytic Dupire local variance σ²_loc(k, T) from eSSVI params
        Real localVariance(Real logMoneyness, Real maturity) const;

        //! Analytic Dupire local volatility σ_loc(k, T) = sqrt(localVariance)
        Real localVol(Real logMoneyness, Real maturity) const;

        //! Gradient of implied vol w.r.t. native (theta_i, rho_i, psi_i)
        EssviSliceGradient impliedVolGradient(Size sliceIdx,
                                              Real logMoneyness) const;

        /*! Gradient of implied vol w.r.t. the full global parameter vector.
            Returns a vector of length 3N where N = numSlices().
            Layout: [drho_1..drho_N, dtheta1, da_2..da_N, dc_1..dc_N].
            Only the slice matching sliceIdx contributes non-trivially
            through the globalToSlice Jacobian.
        */
        std::vector<Real> impliedVolGlobalGradient(
            Size sliceIdx, Real logMoneyness,
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;

        //! Convert global params to native slice params (Proposition 3.1).
        static std::vector<EssviSliceParams>
        globalToSlice(const EssviGlobalParams& gp,
                      EssviButterflyCondition::Type cond
                          = EssviButterflyCondition::GatheralJacquier);

    private:
        EssviSliceParams interpolate(Real t) const;

        std::vector<Real>            T_;
        std::vector<EssviSliceParams> slices_;
        EssviButterflyCondition::Type bflyCond_;
    };

} // namespace QuantLib

#endif
