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

        /*! Chain Jacobian: d(theta_s, rho_s, psi_s) / d(global_param_j)
            for all slices.  Computed once, then reused for all observations.
            Returns flat row-major: N rows × nParams cols, where each row
            contains [dTheta_s/dg_0, ..., dTheta_s/dg_{P-1},
                      dPsi_s/dg_0, ..., dPsi_s/dg_{P-1}].
            Layout: first N*nParams entries are dTheta, next N*nParams are dPsi.
            Total size: 2 * N * nParams.  (rho is direct: drho_s/drho_s = 1)
        */
        std::vector<Real> chainJacobian(
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

    // ---------------------------------------------------------------
    // Dual-wing eSSVI: separate psi for k < 0 and k >= 0
    //
    // w(k) = essvi(k; theta, rho, psi_lo)  if k < 0
    //       = essvi(k; theta, rho, psi_hi)  if k >= 0
    //
    // C^0 at k=0 (w(0) = theta regardless of psi), but NOT C^1
    // unless psi_lo == psi_hi or rho == 0.
    // ---------------------------------------------------------------

    //! Native dual-wing SSVI params for a single expiry slice
    struct DualWingEssviSliceParams {
        Real theta;   //!< ATM total implied variance, > 0
        Real rho;     //!< correlation / skew, in (-1, 1)
        Real psi_lo;  //!< curvature for k < 0 (put wing), > 0
        Real psi_hi;  //!< curvature for k >= 0 (call wing), > 0
    };

    //! Global (arb-free) parameter set for dual-wing eSSVI
    /*! Layout: [rhos(N), theta1(1), as(N-1), cs_lo(N), cs_hi(N)] = 4N params.
        The theta chain and rho are shared; each wing has its own c-chain
        producing independent psi_lo and psi_hi sequences.
    */
    struct DualWingEssviGlobalParams {
        std::vector<Real> rhos;     //!< rho_1..rho_N, each in (-1,1)
        Real              theta1;   //!< theta for the first slice, > 0
        std::vector<Real> as;       //!< a_2..a_N, each > 0
        std::vector<Real> cs_lo;    //!< c_lo_1..c_lo_N, each in (0,1)
        std::vector<Real> cs_hi;    //!< c_hi_1..c_hi_N, each in (0,1)

        Size numSlices() const { return rhos.size(); }
    };

    //! Dual-wing eSSVI total variance: branches on sign(k)
    inline Real dualWingEssviTotalVariance(const DualWingEssviSliceParams& s, Real k) {
        Real psi = (k < 0.0) ? s.psi_lo : s.psi_hi;
        Real inner = psi * k + s.theta * s.rho;
        return 0.5 * (s.theta + s.rho * psi * k
                       + std::sqrt(inner * inner
                                   + s.theta * s.theta * (1.0 - s.rho * s.rho)));
    }

    //! Gradient of dual-wing implied vol w.r.t. (theta, rho, psi_lo, psi_hi)
    struct DualWingEssviSliceGradient {
        Real dSigma_dTheta;
        Real dSigma_dRho;
        Real dSigma_dPsiLo;
        Real dSigma_dPsiHi;
    };

    /*! Analytic gradient of dual-wing total variance w.r.t. (theta, rho, psi_lo, psi_hi).
        Only the active wing's psi gradient is non-zero.
    */
    inline void dualWingEssviTotalVarianceGradient(
            const DualWingEssviSliceParams& s, Real k,
            Real& dw_dTheta, Real& dw_dRho,
            Real& dw_dPsiLo, Real& dw_dPsiHi) {
        Real psi = (k < 0.0) ? s.psi_lo : s.psi_hi;
        Real A = psi * k + s.theta * s.rho;
        Real D = A * A + s.theta * s.theta * (1.0 - s.rho * s.rho);
        Real sqrtD = std::sqrt(D);

        dw_dTheta = 0.5 * (1.0 + (A * s.rho + s.theta * (1.0 - s.rho * s.rho)) / sqrtD);
        dw_dRho   = 0.5 * psi * k * (1.0 + s.theta / sqrtD);

        Real dw_dPsi = 0.5 * k * (s.rho + A / sqrtD);
        if (k < 0.0) {
            dw_dPsiLo = dw_dPsi;
            dw_dPsiHi = 0.0;
        } else {
            dw_dPsiLo = 0.0;
            dw_dPsiHi = dw_dPsi;
        }
    }

    /*! Analytic gradient of dual-wing implied vol sigma = sqrt(w/T). */
    inline DualWingEssviSliceGradient dualWingEssviImpliedVolGradient(
            const DualWingEssviSliceParams& s, Real k, Real T) {
        Real w = dualWingEssviTotalVariance(s, k);
        Real sigma = std::sqrt(w / T);
        Real denom = 2.0 * T * sigma;

        Real dw_dTheta, dw_dRho, dw_dPsiLo, dw_dPsiHi;
        dualWingEssviTotalVarianceGradient(s, k, dw_dTheta, dw_dRho, dw_dPsiLo, dw_dPsiHi);

        return { dw_dTheta / denom,
                 dw_dRho   / denom,
                 dw_dPsiLo / denom,
                 dw_dPsiHi / denom };
    }

    /*! Strike derivatives for dual-wing eSSVI (for local vol). */
    inline void dualWingEssviTotalVarianceStrikeDerivatives(
            const DualWingEssviSliceParams& s, Real k,
            Real& dwdk, Real& d2wdk2) {
        Real psi = (k < 0.0) ? s.psi_lo : s.psi_hi;
        Real A = psi * k + s.theta * s.rho;
        Real D = A * A + s.theta * s.theta * (1.0 - s.rho * s.rho);
        Real sqrtD = std::sqrt(D);

        dwdk   = 0.5 * (s.rho * psi + psi * A / sqrtD);
        d2wdk2 = 0.5 * psi * psi * s.theta * s.theta
                 * (1.0 - s.rho * s.rho) / (D * sqrtD);
    }

    // ---------------------------------------------------------------
    // Dual-wing eSSVI core surface
    // ---------------------------------------------------------------

    class DualWingEssviSurface {
    public:
        //! Construct from global (arb-free) dual-wing parameters.
        DualWingEssviSurface(const std::vector<Real>& maturities,
                             const DualWingEssviGlobalParams& gparams,
                             EssviButterflyCondition::Type bflyCond
                                 = EssviButterflyCondition::GatheralJacquier);

        //! Construct from per-slice native parameters.
        DualWingEssviSurface(const std::vector<Real>& maturities,
                             const std::vector<DualWingEssviSliceParams>& slices);

        //! Total implied variance w(k, T_i) for slice index i
        Real totalVariance(Size sliceIdx, Real logMoneyness) const;

        //! Total implied variance with maturity interpolation
        Real totalVariance(Real logMoneyness, Real maturity) const;

        //! Implied volatility sigma(k, T)
        Real impliedVol(Real logMoneyness, Real maturity) const;

        Size numSlices() const { return slices_.size(); }
        const std::vector<Real>& maturities() const { return T_; }
        const std::vector<DualWingEssviSliceParams>& slices() const { return slices_; }
        const DualWingEssviSliceParams& slice(Size i) const { return slices_.at(i); }

        //! Per-slice gradient w.r.t. native (theta, rho, psi_lo, psi_hi)
        DualWingEssviSliceGradient impliedVolGradient(Size sliceIdx,
                                                       Real logMoneyness) const;

        /*! Gradient of implied vol w.r.t. the full global parameter vector.
            Returns a vector of length 4N.
            Layout: [drho_1..drho_N, dtheta1, da_2..da_N, dc_lo_1..dc_lo_N, dc_hi_1..dc_hi_N].
        */
        std::vector<Real> impliedVolGlobalGradient(
            Size sliceIdx, Real logMoneyness,
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;

        /*! Chain Jacobian for dual-wing: d(theta_s, psi_lo_s, psi_hi_s) / d(global_param_j).
            Layout: [dTheta(N×P), dPsiLo(N×P), dPsiHi(N×P)]. Total: 3*N*P.
        */
        std::vector<Real> chainJacobian(
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;

        //! Convert global params to native slice params.
        static std::vector<DualWingEssviSliceParams>
        globalToSlice(const DualWingEssviGlobalParams& gp,
                      EssviButterflyCondition::Type cond
                          = EssviButterflyCondition::GatheralJacquier);

    private:
        DualWingEssviSliceParams interpolate(Real t) const;

        std::vector<Real>                     T_;
        std::vector<DualWingEssviSliceParams> slices_;
        EssviButterflyCondition::Type         bflyCond_;
    };

    // ---------------------------------------------------------------
    // Split-rho eSSVI: separate (rho, psi) for each wing
    //
    // w(k) = essvi(k; theta, rho_lo, psi_lo)  if k < 0
    //       = essvi(k; theta, rho_hi, psi_hi)  if k >= 0
    //
    // 5 native params per slice.  Each wing is independently arb-free
    // via its own Mingone c-chain.  Theta chain uses
    // p_i = max(p_lo_i, p_hi_i) to satisfy both wings.
    // ---------------------------------------------------------------

    //! Native split-rho SSVI params for a single expiry slice
    struct SplitRhoEssviSliceParams {
        Real theta;    //!< ATM total implied variance, > 0
        Real rho_lo;   //!< put wing correlation, in (-1, 1)
        Real psi_lo;   //!< put wing curvature, > 0
        Real rho_hi;   //!< call wing correlation, in (-1, 1)
        Real psi_hi;   //!< call wing curvature, > 0
    };

    //! Global (arb-free) parameter set for split-rho eSSVI
    /*! Layout: [rhos_lo(N), rhos_hi(N), theta1(1), as(N-1), cs_lo(N), cs_hi(N)] = 5N.
        Each wing has its own rho and c-chain; theta chain shared via max(p_lo, p_hi).
    */
    struct SplitRhoEssviGlobalParams {
        std::vector<Real> rhos_lo;   //!< put wing rho_1..rho_N
        std::vector<Real> rhos_hi;   //!< call wing rho_1..rho_N
        Real              theta1;    //!< theta for the first slice, > 0
        std::vector<Real> as;        //!< a_2..a_N, each > 0
        std::vector<Real> cs_lo;     //!< put wing c_1..c_N, each in (0,1)
        std::vector<Real> cs_hi;     //!< call wing c_1..c_N, each in (0,1)

        Size numSlices() const { return rhos_lo.size(); }
    };

    //! Split-rho eSSVI total variance: branches on sign(k)
    inline Real splitRhoEssviTotalVariance(const SplitRhoEssviSliceParams& s, Real k) {
        Real rho = (k < 0.0) ? s.rho_lo : s.rho_hi;
        Real psi = (k < 0.0) ? s.psi_lo : s.psi_hi;
        Real inner = psi * k + s.theta * rho;
        return 0.5 * (s.theta + rho * psi * k
                       + std::sqrt(inner * inner
                                   + s.theta * s.theta * (1.0 - rho * rho)));
    }

    //! Gradient of split-rho implied vol w.r.t. (theta, rho_lo, psi_lo, rho_hi, psi_hi)
    struct SplitRhoEssviSliceGradient {
        Real dSigma_dTheta;
        Real dSigma_dRhoLo;
        Real dSigma_dPsiLo;
        Real dSigma_dRhoHi;
        Real dSigma_dPsiHi;
    };

    /*! Analytic gradient of split-rho total variance. Only active wing contributes. */
    inline void splitRhoEssviTotalVarianceGradient(
            const SplitRhoEssviSliceParams& s, Real k,
            Real& dw_dTheta, Real& dw_dRhoLo, Real& dw_dPsiLo,
            Real& dw_dRhoHi, Real& dw_dPsiHi) {
        Real rho = (k < 0.0) ? s.rho_lo : s.rho_hi;
        Real psi = (k < 0.0) ? s.psi_lo : s.psi_hi;
        Real A = psi * k + s.theta * rho;
        Real D = A * A + s.theta * s.theta * (1.0 - rho * rho);
        Real sqrtD = std::sqrt(D);

        dw_dTheta = 0.5 * (1.0 + (A * rho + s.theta * (1.0 - rho * rho)) / sqrtD);
        Real dw_dRho = 0.5 * psi * k * (1.0 + s.theta / sqrtD);
        Real dw_dPsi = 0.5 * k * (rho + A / sqrtD);

        if (k < 0.0) {
            dw_dRhoLo = dw_dRho; dw_dPsiLo = dw_dPsi;
            dw_dRhoHi = 0.0;     dw_dPsiHi = 0.0;
        } else {
            dw_dRhoLo = 0.0;     dw_dPsiLo = 0.0;
            dw_dRhoHi = dw_dRho; dw_dPsiHi = dw_dPsi;
        }
    }

    /*! Analytic gradient of split-rho implied vol sigma = sqrt(w/T). */
    inline SplitRhoEssviSliceGradient splitRhoEssviImpliedVolGradient(
            const SplitRhoEssviSliceParams& s, Real k, Real T) {
        Real w = splitRhoEssviTotalVariance(s, k);
        Real sigma = std::sqrt(w / T);
        Real denom = 2.0 * T * sigma;

        Real dw_dTheta, dw_dRhoLo, dw_dPsiLo, dw_dRhoHi, dw_dPsiHi;
        splitRhoEssviTotalVarianceGradient(s, k,
            dw_dTheta, dw_dRhoLo, dw_dPsiLo, dw_dRhoHi, dw_dPsiHi);

        return { dw_dTheta / denom,
                 dw_dRhoLo / denom, dw_dPsiLo / denom,
                 dw_dRhoHi / denom, dw_dPsiHi / denom };
    }

    /*! Strike derivatives for split-rho eSSVI. */
    inline void splitRhoEssviTotalVarianceStrikeDerivatives(
            const SplitRhoEssviSliceParams& s, Real k,
            Real& dwdk, Real& d2wdk2) {
        Real rho = (k < 0.0) ? s.rho_lo : s.rho_hi;
        Real psi = (k < 0.0) ? s.psi_lo : s.psi_hi;
        Real A = psi * k + s.theta * rho;
        Real D = A * A + s.theta * s.theta * (1.0 - rho * rho);
        Real sqrtD = std::sqrt(D);

        dwdk   = 0.5 * (rho * psi + psi * A / sqrtD);
        d2wdk2 = 0.5 * psi * psi * s.theta * s.theta
                 * (1.0 - rho * rho) / (D * sqrtD);
    }

    // ---------------------------------------------------------------
    // Split-rho eSSVI core surface
    // ---------------------------------------------------------------

    class SplitRhoEssviSurface {
    public:
        SplitRhoEssviSurface(const std::vector<Real>& maturities,
                              const SplitRhoEssviGlobalParams& gparams,
                              EssviButterflyCondition::Type bflyCond
                                  = EssviButterflyCondition::GatheralJacquier);

        SplitRhoEssviSurface(const std::vector<Real>& maturities,
                              const std::vector<SplitRhoEssviSliceParams>& slices);

        Real totalVariance(Size sliceIdx, Real logMoneyness) const;
        Real totalVariance(Real logMoneyness, Real maturity) const;
        Real impliedVol(Real logMoneyness, Real maturity) const;

        Size numSlices() const { return slices_.size(); }
        const std::vector<Real>& maturities() const { return T_; }
        const std::vector<SplitRhoEssviSliceParams>& slices() const { return slices_; }
        const SplitRhoEssviSliceParams& slice(Size i) const { return slices_.at(i); }

        SplitRhoEssviSliceGradient impliedVolGradient(Size sliceIdx,
                                                       Real logMoneyness) const;

        /*! Gradient w.r.t. 5N global params.
            Layout: [rhos_lo(N), rhos_hi(N), theta1(1), as(N-1), cs_lo(N), cs_hi(N)].
        */
        std::vector<Real> impliedVolGlobalGradient(
            Size sliceIdx, Real logMoneyness,
            const SplitRhoEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond
                = EssviButterflyCondition::GatheralJacquier) const;

        static std::vector<SplitRhoEssviSliceParams>
        globalToSlice(const SplitRhoEssviGlobalParams& gp,
                      EssviButterflyCondition::Type cond
                          = EssviButterflyCondition::GatheralJacquier);

    private:
        SplitRhoEssviSliceParams interpolate(Real t) const;

        std::vector<Real>                      T_;
        std::vector<SplitRhoEssviSliceParams>  slices_;
        EssviButterflyCondition::Type          bflyCond_;
    };

} // namespace QuantLib

#endif
