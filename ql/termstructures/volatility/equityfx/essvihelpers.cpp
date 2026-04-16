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

#include <ql/termstructures/volatility/equityfx/essvihelpers.hpp>
#include <sstream>

namespace QuantLib {

    // =================================================================
    // Martini-Mingone butterfly bound: numerical infimum
    // =================================================================

    Real essviButterflyBoundMM(Real theta, Real absRho) {
        if (absRho >= 1.0 - 1e-12)
            absRho = 1.0 - 1e-12;

        const Real l2 = essviMmL2(absRho);
        const Real eps = 1e-8;

        Real searchWidth = 50.0;
        Real lLo = l2 + eps;
        Real lHi = lLo + searchWidth;

        // Phase 1: coarse grid (200 points)
        const Size nGrid = 200;
        Real bestL = lLo;
        Real bestR = essviMmRatio(lLo, theta, absRho);

        for (Size i = 1; i <= nGrid; ++i) {
            Real l = lLo + (lHi - lLo) * static_cast<Real>(i) / nGrid;
            Real R = essviMmRatio(l, theta, absRho);
            if (R < bestR) {
                bestR = R;
                bestL = l;
            }
        }

        // If minimum at right edge, extend search
        for (int ext = 0; ext < 3; ++ext) {
            if (bestL < lHi - (lHi - lLo) * 0.05) break;
            lLo = lHi;
            lHi = lLo + searchWidth;
            for (Size i = 0; i <= nGrid; ++i) {
                Real l = lLo + (lHi - lLo) * static_cast<Real>(i) / nGrid;
                Real R = essviMmRatio(l, theta, absRho);
                if (R < bestR) {
                    bestR = R;
                    bestL = l;
                }
            }
        }

        // Phase 2: golden-section refinement
        {
            const Real gr = 0.6180339887498949;
            Real a = std::max(l2 + eps, bestL - 2.0);
            Real b = bestL + 2.0;

            Real c = b - gr * (b - a);
            Real d = a + gr * (b - a);

            for (int iter = 0; iter < 80; ++iter) {
                Real fc = essviMmRatio(c, theta, absRho);
                Real fd = essviMmRatio(d, theta, absRho);
                if (fc < fd) {
                    b = d;
                } else {
                    a = c;
                }
                c = b - gr * (b - a);
                d = a + gr * (b - a);
                if (b - a < 1e-12) break;
            }
            bestR = essviMmRatio(0.5 * (a + b), theta, absRho);
        }

        return bestR;
    }

    // =================================================================
    // Unified sliceButterflyBound
    // =================================================================

    Real essviSliceButterflyBound(Real theta, Real absRho,
                                  EssviButterflyCondition::Type cond) {
        Real leeUpper = essviPsiUpperLee(absRho);
        Real fBfly;
        if (cond == EssviButterflyCondition::GatheralJacquier) {
            fBfly = essviButterflyBoundGJ(theta, absRho);
        } else {
            fBfly = essviButterflyBoundMM(theta, absRho);
        }
        return std::min(leeUpper, std::sqrt(fBfly));
    }

    // =================================================================
    // globalToSlice: Proposition 3.1 sequential parameter recovery
    // =================================================================

    std::vector<EssviSliceParams>
    EssviSurface::globalToSlice(const EssviGlobalParams& gp,
                                 EssviButterflyCondition::Type cond) {
        const Size N = gp.numSlices();
        QL_REQUIRE(N > 0, "EssviSurface: empty parameter set");
        QL_REQUIRE(gp.as.size() == N - 1,
                   "EssviSurface: as.size() (" << gp.as.size()
                   << ") must be N-1 (" << (N - 1) << ")");
        QL_REQUIRE(gp.cs.size() == N,
                   "EssviSurface: cs.size() (" << gp.cs.size()
                   << ") must be N (" << N << ")");

        // validate domains
        for (Size i = 0; i < N; ++i) {
            QL_REQUIRE(gp.rhos[i] > -1.0 && gp.rhos[i] < 1.0,
                       "rho[" << i << "] = " << gp.rhos[i]
                       << " must be in (-1, 1)");
            QL_REQUIRE(gp.cs[i] > 0.0 && gp.cs[i] < 1.0,
                       "c[" << i << "] = " << gp.cs[i]
                       << " must be in (0, 1)");
        }
        QL_REQUIRE(gp.theta1 > 0.0,
                   "theta1 = " << gp.theta1 << " must be > 0");
        for (Size i = 0; i < N - 1; ++i) {
            QL_REQUIRE(gp.as[i] > 0.0,
                       "a[" << i << "] = " << gp.as[i] << " must be > 0");
        }

        std::vector<EssviSliceParams> out(N);
        std::vector<Real> theta(N), psi(N), p(N), f(N);

        // Step 1: compute thetas
        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            p[i] = essviCalendarP(gp.rhos[i - 1], gp.rhos[i]);
            theta[i] = theta[i - 1] * p[i] + gp.as[i - 1];
        }

        // Step 2: compute f_i (butterfly bounds)
        for (Size i = 0; i < N; ++i) {
            f[i] = essviSliceButterflyBound(theta[i], std::abs(gp.rhos[i]), cond);
        }

        // Step 3: compute psi via c parameters
        std::vector<Real> Cpsi(N), Apsi(N);

        auto forwardTailMin = [&](Size i) -> Real {
            Real result = std::numeric_limits<Real>::max();
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= p[j];
                result = std::min(result, f[j] / prodP);
            }
            return result;
        };

        // Slice 0
        Apsi[0] = 0.0;
        {
            Real tailMin = forwardTailMin(0);
            Cpsi[0] = std::min(f[0], tailMin);
        }
        psi[0] = gp.cs[0] * (Cpsi[0] - Apsi[0]) + Apsi[0];

        // Slices 1 .. N-1
        for (Size i = 1; i < N; ++i) {
            Apsi[i] = psi[i - 1] * p[i];

            Real calBound = (psi[i - 1] / theta[i - 1]) * theta[i];
            Real tailMin  = forwardTailMin(i);
            Cpsi[i] = std::min({calBound, f[i], tailMin});

            QL_REQUIRE(Cpsi[i] > Apsi[i],
                       "EssviSurface::globalToSlice: infeasible at slice " << i
                       << " (Cpsi=" << Cpsi[i] << " <= Apsi=" << Apsi[i] << ")");

            psi[i] = gp.cs[i] * (Cpsi[i] - Apsi[i]) + Apsi[i];
        }

        // Assemble output
        for (Size i = 0; i < N; ++i) {
            out[i].theta = theta[i];
            out[i].rho   = gp.rhos[i];
            out[i].psi   = psi[i];
        }
        return out;
    }

    // =================================================================
    // Strike derivatives of total variance
    // =================================================================

    void EssviSurface::totalVarianceStrikeDerivatives(
            Real k, Real t, Real& dwdk, Real& d2wdk2) const {
        EssviSliceParams st = interpolate(t);
        essviTotalVarianceStrikeDerivatives(st, k, dwdk, d2wdk2);
    }

    // =================================================================
    // Time derivative of total variance via param interpolation
    // =================================================================

    Real EssviSurface::totalVarianceTimeDerivative(Real k, Real t) const {
        // ∂w/∂T = ∂w/∂θ · dθ/dT + ∂w/∂ρ · dρ/dT + ∂w/∂ψ · dψ/dT
        // where dθ/dT, dρ/dT, dψ/dT come from the piecewise-linear interpolation.
        const Size N = T_.size();

        Real dTheta_dT, dPsi_dT, dPsiRho_dT;

        if (t <= T_[0]) {
            // Left extrapolation: θ(t) = t/T[0] * θ[0], etc.
            dTheta_dT = slices_[0].theta / T_[0];
            dPsi_dT   = slices_[0].psi   / T_[0];
            dPsiRho_dT = slices_[0].psi * slices_[0].rho / T_[0];
        } else if (t >= T_[N-1]) {
            // Right extrapolation: linear theta extension
            Real slope = (N >= 2)
                ? (slices_[N-1].theta - slices_[N-2].theta) / (T_[N-1] - T_[N-2])
                : slices_[0].theta / T_[0];
            dTheta_dT = slope;
            dPsi_dT   = 0.0;  // psi flat in right extrapolation
            dPsiRho_dT = 0.0;
        } else {
            // Interior: find interval
            auto it = std::upper_bound(T_.begin(), T_.end(), t);
            Size hi = static_cast<Size>(it - T_.begin());
            Size lo = hi - 1;
            Real dT = T_[hi] - T_[lo];

            dTheta_dT  = (slices_[hi].theta - slices_[lo].theta) / dT;
            dPsi_dT    = (slices_[hi].psi   - slices_[lo].psi)   / dT;
            dPsiRho_dT = (slices_[hi].psi * slices_[hi].rho
                          - slices_[lo].psi * slices_[lo].rho) / dT;
        }

        // Interpolated params at t
        EssviSliceParams st = interpolate(t);

        // dρ/dT from the ψρ-blending: ρ = (ψρ)/ψ
        // dρ/dT = (d(ψρ)/dT - ρ · dψ/dT) / ψ
        Real dRho_dT = 0.0;
        if (st.psi > 1e-15) {
            dRho_dT = (dPsiRho_dT - st.rho * dPsi_dT) / st.psi;
        }

        // ∂w/∂(θ,ρ,ψ) at (k, st)
        Real dw_dTheta, dw_dRho, dw_dPsi;
        essviTotalVarianceGradient(st, k, dw_dTheta, dw_dRho, dw_dPsi);

        return dw_dTheta * dTheta_dT + dw_dRho * dRho_dT + dw_dPsi * dPsi_dT;
    }

    // =================================================================
    // Analytic Dupire local vol from eSSVI
    // =================================================================

    Real EssviSurface::localVariance(Real k, Real t) const {
        QL_REQUIRE(t > 0.0, "localVariance: maturity must be > 0");

        Real w = totalVariance(k, t);

        Real dwdk, d2wdk2;
        totalVarianceStrikeDerivatives(k, t, dwdk, d2wdk2);

        Real dwdt = totalVarianceTimeDerivative(k, t);

        // Dupire formula:
        // σ²_loc = dwdt / (1 - k/w*dwdk + ¼(-¼ - 1/w + k²/w²)*dwdk² + ½*d2wdk2)
        if (dwdk == 0.0 && d2wdk2 == 0.0) {
            return dwdt;
        }

        Real den1 = 1.0 - k / w * dwdk;
        Real den2 = 0.25 * (-0.25 - 1.0/w + k*k/(w*w)) * dwdk * dwdk;
        Real den3 = 0.5 * d2wdk2;
        Real den = den1 + den2 + den3;

        QL_REQUIRE(den > 0.0,
                   "EssviSurface::localVariance: non-positive denominator "
                   << den << " at k=" << k << ", T=" << t);

        return dwdt / den;
    }

    Real EssviSurface::localVol(Real k, Real t) const {
        Real lv2 = localVariance(k, t);
        QL_REQUIRE(lv2 >= 0.0,
                   "EssviSurface::localVol: negative local variance "
                   << lv2 << " at k=" << k << ", T=" << t);
        return std::sqrt(lv2);
    }

    // =================================================================
    // Per-slice gradient
    // =================================================================

    EssviSliceGradient EssviSurface::impliedVolGradient(Size idx, Real k) const {
        return essviImpliedVolGradient(slices_.at(idx), k, T_.at(idx));
    }

    // =================================================================
    // Chain Jacobian: d(theta_s, psi_s) / d(global_param_j) for all slices.
    // Computed once, reused for all observations.
    // Layout: [dTheta[0..N-1][0..P-1], dPsi[0..N-1][0..P-1]]
    // Total: 2 * N * P entries, row-major.
    // =================================================================

    std::vector<Real> EssviSurface::chainJacobian(
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {

        const Size N = gp.numSlices();
        const Size P = 3 * N;

        // ── Replay globalToSlice recording intermediates ──
        std::vector<Real> theta(N), psi(N), p(N), f(N);
        std::vector<Real> Apsi(N), Cpsi(N);
        std::vector<int> p_branch(N, -1), Cpsi_branch(N, -1), f_branch(N, -1);

        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            Real a = (1.0 + gp.rhos[i-1]) / (1.0 + gp.rhos[i]);
            Real b = (1.0 - gp.rhos[i-1]) / (1.0 - gp.rhos[i]);
            p[i] = std::max(a, b);
            p_branch[i] = (a >= b) ? 0 : 1;
            theta[i] = theta[i-1] * p[i] + gp.as[i-1];
        }

        for (Size i = 0; i < N; ++i) {
            Real absRho = std::abs(gp.rhos[i]);
            Real lee = essviPsiUpperLee(absRho);
            Real fBfly = (bflyCond == EssviButterflyCondition::GatheralJacquier)
                ? essviButterflyBoundGJ(theta[i], absRho)
                : essviButterflyBoundMM(theta[i], absRho);
            Real sqrtF = std::sqrt(fBfly);
            if (lee <= sqrtF) { f[i] = lee; f_branch[i] = 0; }
            else              { f[i] = sqrtF; f_branch[i] = 1; }
        }

        auto forwardTailMinIdx = [&](Size i) -> std::pair<Real, int> {
            Real result = std::numeric_limits<Real>::max();
            int bestJ = -1;
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= p[j];
                Real val = f[j] / prodP;
                if (val < result) { result = val; bestJ = static_cast<int>(j); }
            }
            return {result, bestJ};
        };

        Apsi[0] = 0.0;
        {
            auto [tailMin, tailJ] = forwardTailMinIdx(0);
            if (f[0] <= tailMin) { Cpsi[0] = f[0]; Cpsi_branch[0] = 0; }
            else                 { Cpsi[0] = tailMin; Cpsi_branch[0] = tailJ; }
        }
        psi[0] = gp.cs[0] * (Cpsi[0] - Apsi[0]) + Apsi[0];

        for (Size i = 1; i < N; ++i) {
            Apsi[i] = psi[i-1] * p[i];
            Real calBound = (psi[i-1] / theta[i-1]) * theta[i];
            auto [tailMin, tailJ] = forwardTailMinIdx(i);
            Real Cval = calBound; int branch = -1;
            if (f[i] < Cval) { Cval = f[i]; branch = static_cast<int>(i); }
            if (tailMin < Cval) { Cval = tailMin; branch = tailJ; }
            Cpsi[i] = Cval; Cpsi_branch[i] = branch;
            psi[i] = gp.cs[i] * (Cpsi[i] - Apsi[i]) + Apsi[i];
        }

        // ── Forward-mode AD: compute dTheta_s/dg and dPsi_s/dg for all slices ──
        Size rhoOff = 0, th1Off = N, aOff = N + 1, cOff = 2 * N;

        // Output: dTheta[s * P + g] and dPsi[s * P + g]
        std::vector<Real> result(2 * N * P, 0.0);
        Real* dTheta = result.data();          // first N*P entries
        Real* dPsi   = result.data() + N * P;  // next N*P entries

        for (Size gi = 0; gi < P; ++gi) {
            auto dp_dg = [&](Size i) -> Real {
                if (i == 0) return 0.0;
                Real dp_drho_prev, dp_drho_curr;
                if (p_branch[i] == 0) {
                    dp_drho_prev = 1.0 / (1.0 + gp.rhos[i]);
                    dp_drho_curr = -(1.0 + gp.rhos[i-1]) / ((1.0 + gp.rhos[i]) * (1.0 + gp.rhos[i]));
                } else {
                    dp_drho_prev = -1.0 / (1.0 - gp.rhos[i]);
                    dp_drho_curr = (1.0 - gp.rhos[i-1]) / ((1.0 - gp.rhos[i]) * (1.0 - gp.rhos[i]));
                }
                Real r = 0.0;
                if (gi == rhoOff + i - 1) r += dp_drho_prev;
                if (gi == rhoOff + i)     r += dp_drho_curr;
                return r;
            };

            // Theta tangents for all slices
            std::vector<Real> dTh(N, 0.0);
            dTh[0] = (gi == th1Off) ? 1.0 : 0.0;
            for (Size j = 1; j < N; ++j) {
                dTh[j] = dTh[j-1] * p[j] + theta[j-1] * dp_dg(j);
                if (gi == aOff + (j - 1)) dTh[j] += 1.0;
            }
            for (Size s = 0; s < N; ++s)
                dTheta[s * P + gi] = dTh[s];

            auto df_dg = [&](Size j) -> Real {
                Real absRho = std::abs(gp.rhos[j]);
                Real signRho = (gp.rhos[j] >= 0.0) ? 1.0 : -1.0;
                if (f_branch[j] == 0) {
                    if (gi == rhoOff + j)
                        return -4.0 * signRho / ((1.0 + absRho) * (1.0 + absRho));
                    return 0.0;
                } else {
                    Real df_dt = 2.0 / (f[j] * (1.0 + absRho));
                    Real df_dr = -2.0 * theta[j] * signRho
                                 / (f[j] * (1.0 + absRho) * (1.0 + absRho));
                    Real r = df_dt * dTh[j];
                    if (gi == rhoOff + j) r += df_dr;
                    return r;
                }
            };

            // Psi tangents for all slices
            Real dPsi_prev = 0.0;
            for (Size s = 0; s < N; ++s) {
                Real dApsi_s = (s == 0) ? 0.0
                    : dPsi_prev * p[s] + psi[s-1] * dp_dg(s);

                Real dCpsi_s = 0.0;
                if (Cpsi_branch[s] == -1 && s > 0) {
                    dCpsi_s = (dPsi_prev * theta[s] + psi[s-1] * dTh[s]
                               - psi[s-1] * theta[s] / theta[s-1] * dTh[s-1])
                              / theta[s-1];
                } else if (Cpsi_branch[s] >= 0) {
                    int j = Cpsi_branch[s];
                    if (j == static_cast<int>(s)) {
                        dCpsi_s = df_dg(static_cast<Size>(j));
                    } else {
                        Real prodP = 1.0;
                        for (Size m = s + 1; m <= static_cast<Size>(j); ++m)
                            prodP *= p[m];
                        Real dfj = df_dg(static_cast<Size>(j));
                        Real dProdP = 0.0;
                        for (Size m = s + 1; m <= static_cast<Size>(j); ++m)
                            dProdP += dp_dg(m) * prodP / p[m];
                        dCpsi_s = dfj / prodP - f[j] * dProdP / (prodP * prodP);
                    }
                }

                Real dPsi_s = gp.cs[s] * dCpsi_s + (1.0 - gp.cs[s]) * dApsi_s;
                if (gi == cOff + s) dPsi_s += (Cpsi[s] - Apsi[s]);

                dPsi[s * P + gi] = dPsi_s;
                dPsi_prev = dPsi_s;
            }
        }

        return result;
    }

    // =================================================================
    // Global gradient — Jacobian of globalToSlice via forward-mode AD
    // =================================================================

    std::vector<Real> EssviSurface::impliedVolGlobalGradient(
            Size idx, Real k,
            const EssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {

        const Size N = gp.numSlices();
        QL_REQUIRE(N == slices_.size(),
                   "global params N=" << N << " != slices " << slices_.size());
        QL_REQUIRE(idx < N, "sliceIdx " << idx << " >= N=" << N);

        // ── Step 1: replay globalToSlice, recording intermediates ────
        std::vector<Real> theta(N), psi(N), p(N), f(N);
        std::vector<Real> Apsi(N), Cpsi(N);
        std::vector<int> p_branch(N, -1);
        std::vector<int> Cpsi_branch(N, -1);  // -1=calBound, >=0 = f[j] branch

        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            Real a = (1.0 + gp.rhos[i-1]) / (1.0 + gp.rhos[i]);
            Real b = (1.0 - gp.rhos[i-1]) / (1.0 - gp.rhos[i]);
            p[i] = std::max(a, b);
            p_branch[i] = (a >= b) ? 0 : 1;
            theta[i] = theta[i-1] * p[i] + gp.as[i-1];
        }

        // f_i = min(Lee, sqrt(f_bfly)) and track which branch is active
        // f_branch: 0 = Lee bound (4/(1+|rho|)), 1 = sqrt(GJ or MM)
        std::vector<int> f_branch(N, -1);
        for (Size i = 0; i < N; ++i) {
            Real absRho = std::abs(gp.rhos[i]);
            Real lee = essviPsiUpperLee(absRho);
            Real fBfly;
            if (bflyCond == EssviButterflyCondition::GatheralJacquier)
                fBfly = essviButterflyBoundGJ(theta[i], absRho);
            else
                fBfly = essviButterflyBoundMM(theta[i], absRho);
            Real sqrtF = std::sqrt(fBfly);
            if (lee <= sqrtF) {
                f[i] = lee;
                f_branch[i] = 0;  // Lee bound active
            } else {
                f[i] = sqrtF;
                f_branch[i] = 1;  // sqrt(GJ/MM) active
            }
        }

        // df_i/dtheta_i and df_i/drho_i (for GJ bound, which is the common case)
        // When f_branch[i]==0 (Lee): f = 4/(1+|rho|), df/dtheta = 0
        // When f_branch[i]==1 (GJ):  f = sqrt(4*theta/(1+|rho|))
        //   df/dtheta = 1 / (2*f*(1+|rho|)/4) = 2 / (f*(1+|rho|))  -- actually:
        //   f = sqrt(4θ/(1+|ρ|)), df/dθ = (1/(2f)) * 4/(1+|ρ|) = 2/(f*(1+|ρ|))
        //   df/d|ρ| = (1/(2f)) * (-4θ/(1+|ρ|)^2) = -2θ/(f*(1+|ρ|)^2)
        //   df/dρ = df/d|ρ| * sign(ρ) when ρ != 0

        auto forwardTailMinIdx = [&](Size i) -> std::pair<Real, int> {
            Real result = std::numeric_limits<Real>::max();
            int bestJ = -1;
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= p[j];
                Real val = f[j] / prodP;
                if (val < result) { result = val; bestJ = static_cast<int>(j); }
            }
            return {result, bestJ};
        };

        Apsi[0] = 0.0;
        {
            auto [tailMin, tailJ] = forwardTailMinIdx(0);
            if (f[0] <= tailMin) {
                Cpsi[0] = f[0]; Cpsi_branch[0] = 0;
            } else {
                Cpsi[0] = tailMin; Cpsi_branch[0] = tailJ;
            }
        }
        psi[0] = gp.cs[0] * (Cpsi[0] - Apsi[0]) + Apsi[0];

        for (Size i = 1; i < N; ++i) {
            Apsi[i] = psi[i-1] * p[i];
            Real calBound = (psi[i-1] / theta[i-1]) * theta[i];
            auto [tailMin, tailJ] = forwardTailMinIdx(i);

            Real Cval = calBound; int branch = -1;
            if (f[i] < Cval) { Cval = f[i]; branch = static_cast<int>(i); }
            if (tailMin < Cval) { Cval = tailMin; branch = tailJ; }
            Cpsi[i] = Cval; Cpsi_branch[i] = branch;
            psi[i] = gp.cs[i] * (Cpsi[i] - Apsi[i]) + Apsi[i];
        }

        // ── Step 2: per-slice gradient ∂σ/∂(native theta_idx, rho_idx, psi_idx) ──
        EssviSliceGradient sliceGrad = essviImpliedVolGradient(slices_.at(idx), k, T_.at(idx));

        // ── Step 3: Jacobian via forward-mode AD ─────────────────────
        // For each global param g, we compute dtheta_idx/dg and dpsi_idx/dg
        // by replaying the sequential recovery with tangent propagation.
        //
        // Global param layout: [rho_0..rho_{N-1}, theta1, a_1..a_{N-1}, c_0..c_{N-1}]
        Size nParams = 3 * N;
        Size rhoOff = 0, th1Off = N, aOff = N + 1, cOff = 2 * N;

        std::vector<Real> grad(nParams, 0.0);

        // ρ_idx is direct
        grad[rhoOff + idx] += sliceGrad.dSigma_dRho;

        // For each global param, propagate tangent through theta and psi chains
        for (Size gi = 0; gi < nParams; ++gi) {
            Real dPsi_prev = 0.0;    // dpsi[i-1]/dg
            Real dTheta_i = 0.0;     // dtheta[i]/dg (set from dTheta_all)

            // p[i] derivatives w.r.t. global param gi
            auto dp_dg = [&](Size i) -> Real {
                if (i == 0) return 0.0;
                Real dp_drho_prev, dp_drho_curr;
                if (p_branch[i] == 0) {
                    dp_drho_prev = 1.0 / (1.0 + gp.rhos[i]);
                    dp_drho_curr = -(1.0 + gp.rhos[i-1]) / ((1.0 + gp.rhos[i]) * (1.0 + gp.rhos[i]));
                } else {
                    dp_drho_prev = -1.0 / (1.0 - gp.rhos[i]);
                    dp_drho_curr = (1.0 - gp.rhos[i-1]) / ((1.0 - gp.rhos[i]) * (1.0 - gp.rhos[i]));
                }
                Real result = 0.0;
                if (gi == rhoOff + i - 1) result += dp_drho_prev;
                if (gi == rhoOff + i)     result += dp_drho_curr;
                return result;
            };

            // We need theta tangents for all slices (not just up to idx)
            // because f[j] for j > idx may determine Cpsi via tail bounds.
            // Pre-compute all dTheta[j]/dg.
            std::vector<Real> dTheta_all(N, 0.0);
            dTheta_all[0] = (gi == th1Off) ? 1.0 : 0.0;
            for (Size j = 1; j < N; ++j) {
                Real dp_j = dp_dg(j);
                dTheta_all[j] = dTheta_all[j-1] * p[j] + theta[j-1] * dp_j;
                if (gi == aOff + (j - 1))
                    dTheta_all[j] += 1.0;
            }

            // f[j] derivative: df[j]/dg
            // f[j] = min(Lee, sqrt(GJ))
            // Lee = 4/(1+|rho_j|) → only depends on rho_j
            // sqrt(GJ) = sqrt(4*theta_j/(1+|rho_j|)) → depends on theta_j and rho_j
            auto df_dg = [&](Size j) -> Real {
                Real absRho = std::abs(gp.rhos[j]);
                Real signRho = (gp.rhos[j] >= 0.0) ? 1.0 : -1.0;
                if (f_branch[j] == 0) {
                    // Lee: f = 4/(1+|rho|)
                    // df/drho_j = -4*sign(rho) / (1+|rho|)^2
                    if (gi == rhoOff + j) {
                        Real denom = (1.0 + absRho) * (1.0 + absRho);
                        return -4.0 * signRho / denom;
                    }
                    return 0.0;
                } else {
                    // sqrt(GJ): f = sqrt(4*theta/(1+|rho|))
                    // df/dtheta = 2 / (f * (1+|rho|))
                    // df/drho   = -2*theta*sign(rho) / (f * (1+|rho|)^2)
                    Real df_dtheta = 2.0 / (f[j] * (1.0 + absRho));
                    Real df_drho   = -2.0 * theta[j] * signRho
                                     / (f[j] * (1.0 + absRho) * (1.0 + absRho));
                    Real result = df_dtheta * dTheta_all[j];
                    if (gi == rhoOff + j)
                        result += df_drho;
                    return result;
                }
            };

            // Propagate through slices 0..idx
            for (Size i = 0; i <= idx; ++i) {
                // ── theta[i] tangent ──
                dTheta_i = dTheta_all[i];

                // ── psi[i] tangent ──
                Real dPsi_i = 0.0;

                // Apsi[i] tangent
                Real dApsi_i;
                if (i == 0) {
                    dApsi_i = 0.0;
                } else {
                    // Apsi[i] = psi[i-1] * p[i]
                    dApsi_i = dPsi_prev * p[i] + psi[i-1] * dp_dg(i);
                }

                // Cpsi[i] tangent (depends on active branch)
                Real dCpsi_i = 0.0;
                if (Cpsi_branch[i] == -1 && i > 0) {
                    // calBound = psi[i-1] / theta[i-1] * theta[i]
                    dCpsi_i = (dPsi_prev * theta[i] + psi[i-1] * dTheta_i
                               - psi[i-1] * theta[i] / theta[i-1] * dTheta_all[i-1])
                              / theta[i-1];
                } else if (Cpsi_branch[i] >= 0) {
                    // f[j] branch (could be f[i] itself or tail f[j]/prod(p))
                    int j = Cpsi_branch[i];
                    if (j == static_cast<int>(i)) {
                        // Cpsi = f[i]
                        dCpsi_i = df_dg(static_cast<Size>(j));
                    } else {
                        // Cpsi = f[j] / prod(p_{i+1}..p_j)
                        Real prodP = 1.0;
                        for (Size m = i + 1; m <= static_cast<Size>(j); ++m)
                            prodP *= p[m];
                        // d/dg (f[j] / prodP) = df[j]/dg / prodP
                        //   - f[j] / prodP^2 * d(prodP)/dg
                        Real dfj = df_dg(static_cast<Size>(j));
                        Real dProdP = 0.0;
                        for (Size m = i + 1; m <= static_cast<Size>(j); ++m) {
                            Real dp_m = dp_dg(m);
                            // d(prodP)/dp_m = prodP / p[m]
                            dProdP += dp_m * prodP / p[m];
                        }
                        dCpsi_i = dfj / prodP - f[j] * dProdP / (prodP * prodP);
                    }
                }

                // psi[i] = c_i * (Cpsi[i] - Apsi[i]) + Apsi[i]
                //        = c_i * Cpsi[i] + (1 - c_i) * Apsi[i]
                dPsi_i = gp.cs[i] * dCpsi_i + (1.0 - gp.cs[i]) * dApsi_i;
                if (gi == cOff + i)
                    dPsi_i += (Cpsi[i] - Apsi[i]);

                // Save for next iteration
                dPsi_prev = dPsi_i;
            }

            // Accumulate: dσ/dg = dσ/dθ_idx * dθ_idx/dg + dσ/dψ_idx * dψ_idx/dg
            grad[gi] += sliceGrad.dSigma_dTheta * dTheta_all[idx]
                      + sliceGrad.dSigma_dPsi   * dPsi_prev;
        }

        // Note: ρ_idx contribution from the eSSVI formula was already added above.
        // The rho also affects theta and psi through p_i, which is captured in
        // the forward-mode loop. But the direct ∂σ/∂ρ_idx (from the eSSVI formula
        // itself) is separate and was added at the start. We need to avoid
        // double-counting: the loop above computes dσ/dg via dθ and dψ only,
        // so the direct ρ term is correctly separate.

        return grad;
    }

    // =================================================================
    // Constructors
    // =================================================================

    EssviSurface::EssviSurface(const std::vector<Real>& maturities,
                                const EssviGlobalParams& gparams,
                                EssviButterflyCondition::Type bflyCond)
        : T_(maturities),
          slices_(globalToSlice(gparams, bflyCond)),
          bflyCond_(bflyCond)
    {
        QL_REQUIRE(T_.size() == slices_.size(),
                   "maturity count (" << T_.size()
                   << ") != slice count (" << slices_.size() << ")");
        for (Size i = 0; i < T_.size(); ++i) {
            QL_REQUIRE(T_[i] > 0.0,
                       "maturity[" << i << "] = " << T_[i] << " must be > 0");
        }
        for (Size i = 1; i < T_.size(); ++i) {
            QL_REQUIRE(T_[i] > T_[i - 1],
                       "maturities must be strictly increasing: T["
                       << i - 1 << "]=" << T_[i-1] << " >= T["
                       << i << "]=" << T_[i]);
        }
    }

    EssviSurface::EssviSurface(const std::vector<Real>& maturities,
                                const std::vector<EssviSliceParams>& slices)
        : T_(maturities), slices_(slices),
          bflyCond_(EssviButterflyCondition::GatheralJacquier)
    {
        QL_REQUIRE(T_.size() == slices_.size(),
                   "maturity count (" << T_.size()
                   << ") != slice count (" << slices_.size() << ")");
        QL_REQUIRE(!T_.empty(), "maturities must not be empty");
        for (Size i = 0; i < T_.size(); ++i) {
            QL_REQUIRE(T_[i] > 0.0,
                       "maturity[" << i << "] = " << T_[i] << " must be > 0");
        }
        for (Size i = 1; i < T_.size(); ++i) {
            QL_REQUIRE(T_[i] > T_[i - 1],
                       "maturities must be strictly increasing: T["
                       << i - 1 << "]=" << T_[i-1] << " >= T["
                       << i << "]=" << T_[i]);
        }
        // Validate native params
        for (Size i = 0; i < slices_.size(); ++i) {
            QL_REQUIRE(slices_[i].theta > 0.0,
                       "theta[" << i << "] = " << slices_[i].theta
                       << " must be > 0");
            QL_REQUIRE(slices_[i].rho > -1.0 && slices_[i].rho < 1.0,
                       "rho[" << i << "] = " << slices_[i].rho
                       << " must be in (-1, 1)");
            QL_REQUIRE(slices_[i].psi > 0.0,
                       "psi[" << i << "] = " << slices_[i].psi
                       << " must be > 0");
        }
    }

    // =================================================================
    // In-place global param update
    // =================================================================

    void EssviSurface::setGlobalParams(const EssviGlobalParams& gp,
                                        EssviButterflyCondition::Type bflyCond) {
        auto newSlices = globalToSlice(gp, bflyCond);
        QL_REQUIRE(newSlices.size() == T_.size(),
                   "setGlobalParams: slice count (" << newSlices.size()
                   << ") != maturity count (" << T_.size() << ")");
        slices_ = std::move(newSlices);
        bflyCond_ = bflyCond;
    }

    // =================================================================
    // Evaluation: slice index
    // =================================================================

    Real EssviSurface::totalVariance(Size idx, Real k) const {
        return essviTotalVariance(slices_.at(idx), k);
    }

    // =================================================================
    // Interpolation (Section 5 of Mingone 2022)
    // =================================================================

    EssviSliceParams EssviSurface::interpolate(Real t) const {
        const Size N = T_.size();

        // Left extrapolation
        if (t <= T_[0]) {
            Real lambda = t / T_[0];
            if (lambda < 0.0) lambda = 0.0;
            return { lambda * slices_[0].theta,
                     slices_[0].rho,
                     lambda * slices_[0].psi };
        }

        // Right extrapolation
        if (t >= T_[N - 1]) {
            Real slope = 0.0;
            if (N >= 2) {
                slope = (slices_[N - 1].theta - slices_[N - 2].theta)
                        / (T_[N - 1] - T_[N - 2]);
            } else {
                slope = slices_[0].theta / T_[0];
            }
            Real thetaT = slices_[N - 1].theta + slope * (t - T_[N - 1]);
            return { thetaT,
                     slices_[N - 1].rho,
                     slices_[N - 1].psi };
        }

        // Interpolation between benchmarks
        auto it = std::upper_bound(T_.begin(), T_.end(), t);
        Size hi = static_cast<Size>(it - T_.begin());
        Size lo = hi - 1;

        Real lambda = (t - T_[lo]) / (T_[hi] - T_[lo]);

        const auto& s1 = slices_[lo];
        const auto& s2 = slices_[hi];

        EssviSliceParams st;
        st.theta = (1.0 - lambda) * s1.theta + lambda * s2.theta;
        st.psi   = (1.0 - lambda) * s1.psi   + lambda * s2.psi;
        Real psiRho = (1.0 - lambda) * s1.psi * s1.rho
                      + lambda * s2.psi * s2.rho;
        st.rho = (st.psi > 1e-15) ? psiRho / st.psi : 0.0;

        return st;
    }

    // =================================================================
    // Evaluation: arbitrary maturity
    // =================================================================

    Real EssviSurface::totalVariance(Real k, Real t) const {
        EssviSliceParams st = interpolate(t);
        return essviTotalVariance(st, k);
    }

    Real EssviSurface::impliedVol(Real k, Real t) const {
        QL_REQUIRE(t > 0.0, "impliedVol: maturity must be > 0");
        Real w = totalVariance(k, t);
        QL_REQUIRE(w >= 0.0, "impliedVol: negative total variance " << w);
        return std::sqrt(w / t);
    }

    // =================================================================
    // DualWingEssviSurface: globalToSlice
    // =================================================================

    std::vector<DualWingEssviSliceParams>
    DualWingEssviSurface::globalToSlice(const DualWingEssviGlobalParams& gp,
                                         EssviButterflyCondition::Type cond) {
        const Size N = gp.numSlices();
        QL_REQUIRE(N > 0, "DualWingEssviSurface: empty parameter set");
        QL_REQUIRE(gp.as.size() == N - 1,
                   "DualWingEssviSurface: as.size() (" << gp.as.size()
                   << ") must be N-1 (" << (N - 1) << ")");
        QL_REQUIRE(gp.cs_lo.size() == N,
                   "DualWingEssviSurface: cs_lo.size() (" << gp.cs_lo.size()
                   << ") must be N (" << N << ")");
        QL_REQUIRE(gp.cs_hi.size() == N,
                   "DualWingEssviSurface: cs_hi.size() (" << gp.cs_hi.size()
                   << ") must be N (" << N << ")");

        for (Size i = 0; i < N; ++i) {
            QL_REQUIRE(gp.rhos[i] > -1.0 && gp.rhos[i] < 1.0,
                       "rho[" << i << "] = " << gp.rhos[i] << " must be in (-1, 1)");
            QL_REQUIRE(gp.cs_lo[i] > 0.0 && gp.cs_lo[i] < 1.0,
                       "cs_lo[" << i << "] = " << gp.cs_lo[i] << " must be in (0, 1)");
            QL_REQUIRE(gp.cs_hi[i] > 0.0 && gp.cs_hi[i] < 1.0,
                       "cs_hi[" << i << "] = " << gp.cs_hi[i] << " must be in (0, 1)");
        }
        QL_REQUIRE(gp.theta1 > 0.0, "theta1 must be > 0");
        for (Size i = 0; i < N - 1; ++i)
            QL_REQUIRE(gp.as[i] > 0.0, "a[" << i << "] must be > 0");

        // Step 1: compute theta chain (shared between wings)
        std::vector<Real> theta(N), p(N, 0.0), f(N);
        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            p[i] = essviCalendarP(gp.rhos[i - 1], gp.rhos[i]);
            theta[i] = theta[i - 1] * p[i] + gp.as[i - 1];
        }

        // Step 2: butterfly bounds (same for both wings — same theta, rho)
        for (Size i = 0; i < N; ++i)
            f[i] = essviSliceButterflyBound(theta[i], std::abs(gp.rhos[i]), cond);

        // Step 3: forward tail helper
        auto forwardTailMin = [&](Size i) -> Real {
            Real result = std::numeric_limits<Real>::max();
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= p[j];
                result = std::min(result, f[j] / prodP);
            }
            return result;
        };

        // Step 4: run two parallel c-chains for psi_lo and psi_hi
        auto runPsiChain = [&](const std::vector<Real>& cs) -> std::vector<Real> {
            std::vector<Real> psi(N), Apsi(N), Cpsi(N);
            Apsi[0] = 0.0;
            Cpsi[0] = std::min(f[0], forwardTailMin(0));
            psi[0] = cs[0] * (Cpsi[0] - Apsi[0]) + Apsi[0];

            for (Size i = 1; i < N; ++i) {
                Apsi[i] = psi[i - 1] * p[i];
                Real calBound = (psi[i - 1] / theta[i - 1]) * theta[i];
                Real tailMin = forwardTailMin(i);
                Cpsi[i] = std::min({calBound, f[i], tailMin});

                QL_REQUIRE(Cpsi[i] > Apsi[i],
                           "DualWingEssviSurface::globalToSlice: infeasible at slice "
                           << i << " (Cpsi=" << Cpsi[i] << " <= Apsi=" << Apsi[i] << ")");
                psi[i] = cs[i] * (Cpsi[i] - Apsi[i]) + Apsi[i];
            }
            return psi;
        };

        std::vector<Real> psi_lo = runPsiChain(gp.cs_lo);
        std::vector<Real> psi_hi = runPsiChain(gp.cs_hi);

        // Assemble output
        std::vector<DualWingEssviSliceParams> out(N);
        for (Size i = 0; i < N; ++i) {
            out[i].theta  = theta[i];
            out[i].rho    = gp.rhos[i];
            out[i].psi_lo = psi_lo[i];
            out[i].psi_hi = psi_hi[i];
        }
        return out;
    }

    // =================================================================
    // DualWingEssviSurface: constructors
    // =================================================================

    DualWingEssviSurface::DualWingEssviSurface(
            const std::vector<Real>& maturities,
            const DualWingEssviGlobalParams& gparams,
            EssviButterflyCondition::Type bflyCond)
        : T_(maturities),
          slices_(globalToSlice(gparams, bflyCond)),
          bflyCond_(bflyCond)
    {
        QL_REQUIRE(T_.size() == slices_.size(),
                   "maturity count != slice count");
        for (Size i = 0; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > 0.0, "maturity[" << i << "] must be > 0");
        for (Size i = 1; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > T_[i - 1], "maturities must be strictly increasing");
    }

    DualWingEssviSurface::DualWingEssviSurface(
            const std::vector<Real>& maturities,
            const std::vector<DualWingEssviSliceParams>& slices)
        : T_(maturities), slices_(slices),
          bflyCond_(EssviButterflyCondition::GatheralJacquier)
    {
        QL_REQUIRE(T_.size() == slices_.size(),
                   "maturity count != slice count");
        QL_REQUIRE(!T_.empty(), "maturities must not be empty");
        for (Size i = 0; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > 0.0, "maturity[" << i << "] must be > 0");
        for (Size i = 1; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > T_[i - 1], "maturities must be strictly increasing");
    }

    // =================================================================
    // DualWingEssviSurface: evaluation
    // =================================================================

    Real DualWingEssviSurface::totalVariance(Size idx, Real k) const {
        return dualWingEssviTotalVariance(slices_.at(idx), k);
    }

    DualWingEssviSliceParams DualWingEssviSurface::interpolate(Real t) const {
        const Size N = T_.size();

        if (t <= T_[0]) {
            Real lambda = std::max(t / T_[0], 0.0);
            return { lambda * slices_[0].theta,
                     slices_[0].rho,
                     lambda * slices_[0].psi_lo,
                     lambda * slices_[0].psi_hi };
        }

        if (t >= T_[N - 1]) {
            Real slope = (N >= 2)
                ? (slices_[N-1].theta - slices_[N-2].theta) / (T_[N-1] - T_[N-2])
                : slices_[0].theta / T_[0];
            Real thetaT = slices_[N - 1].theta + slope * (t - T_[N - 1]);
            return { thetaT,
                     slices_[N - 1].rho,
                     slices_[N - 1].psi_lo,
                     slices_[N - 1].psi_hi };
        }

        auto it = std::upper_bound(T_.begin(), T_.end(), t);
        Size hi = static_cast<Size>(it - T_.begin());
        Size lo = hi - 1;
        Real lambda = (t - T_[lo]) / (T_[hi] - T_[lo]);

        const auto& s1 = slices_[lo];
        const auto& s2 = slices_[hi];

        DualWingEssviSliceParams st;
        st.theta  = (1.0 - lambda) * s1.theta  + lambda * s2.theta;
        st.psi_lo = (1.0 - lambda) * s1.psi_lo + lambda * s2.psi_lo;
        st.psi_hi = (1.0 - lambda) * s1.psi_hi + lambda * s2.psi_hi;
        // Blend rho via psi-weighted average (using avg of lo/hi psi)
        Real avgPsi1 = 0.5 * (s1.psi_lo + s1.psi_hi);
        Real avgPsi2 = 0.5 * (s2.psi_lo + s2.psi_hi);
        Real psiBlend = (1.0 - lambda) * avgPsi1 + lambda * avgPsi2;
        Real psiRho = (1.0 - lambda) * avgPsi1 * s1.rho
                      + lambda * avgPsi2 * s2.rho;
        st.rho = (psiBlend > 1e-15) ? psiRho / psiBlend : 0.0;

        return st;
    }

    Real DualWingEssviSurface::totalVariance(Real k, Real t) const {
        DualWingEssviSliceParams st = interpolate(t);
        return dualWingEssviTotalVariance(st, k);
    }

    Real DualWingEssviSurface::impliedVol(Real k, Real t) const {
        QL_REQUIRE(t > 0.0, "impliedVol: maturity must be > 0");
        Real w = totalVariance(k, t);
        QL_REQUIRE(w >= 0.0, "impliedVol: negative total variance " << w);
        return std::sqrt(w / t);
    }

    // =================================================================
    // DualWingEssviSurface: per-slice gradient
    // =================================================================

    DualWingEssviSliceGradient
    DualWingEssviSurface::impliedVolGradient(Size idx, Real k) const {
        return dualWingEssviImpliedVolGradient(slices_.at(idx), k, T_.at(idx));
    }

    // =================================================================
    // DualWingEssviSurface: factored chain Jacobian
    // Layout: [dTheta(N×P), dPsiLo(N×P), dPsiHi(N×P)]. Total: 3*N*P.
    // =================================================================

    std::vector<Real> DualWingEssviSurface::chainJacobian(
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {

        const Size N = gp.numSlices();
        const Size P = 4 * N;

        // Replay globalToSlice recording intermediates
        std::vector<Real> theta(N), p(N, 0.0), f(N);
        std::vector<int> p_branch(N, -1), f_branch(N, -1);

        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            Real a = (1.0 + gp.rhos[i-1]) / (1.0 + gp.rhos[i]);
            Real b = (1.0 - gp.rhos[i-1]) / (1.0 - gp.rhos[i]);
            p[i] = std::max(a, b);
            p_branch[i] = (a >= b) ? 0 : 1;
            theta[i] = theta[i-1] * p[i] + gp.as[i-1];
        }

        for (Size i = 0; i < N; ++i) {
            Real absRho = std::abs(gp.rhos[i]);
            Real lee = essviPsiUpperLee(absRho);
            Real fBfly = (bflyCond == EssviButterflyCondition::GatheralJacquier)
                ? essviButterflyBoundGJ(theta[i], absRho)
                : essviButterflyBoundMM(theta[i], absRho);
            Real sqrtF = std::sqrt(fBfly);
            if (lee <= sqrtF) { f[i] = lee; f_branch[i] = 0; }
            else              { f[i] = sqrtF; f_branch[i] = 1; }
        }

        auto forwardTailMinIdx = [&](Size i) -> std::pair<Real, int> {
            Real result = std::numeric_limits<Real>::max();
            int bestJ = -1;
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= p[j];
                Real val = f[j] / prodP;
                if (val < result) { result = val; bestJ = static_cast<int>(j); }
            }
            return {result, bestJ};
        };

        // Replay both c-chains to get psi values and branch info
        auto replayChain = [&](const std::vector<Real>& cs)
            -> std::tuple<std::vector<Real>, std::vector<Real>, std::vector<Real>,
                          std::vector<int>> {
            std::vector<Real> psi_out(N), Apsi_out(N), Cpsi_out(N);
            std::vector<int> Cpsi_br(N, -1);
            Apsi_out[0] = 0.0;
            {
                auto [tailMin, tailJ] = forwardTailMinIdx(0);
                if (f[0] <= tailMin) { Cpsi_out[0] = f[0]; Cpsi_br[0] = 0; }
                else                 { Cpsi_out[0] = tailMin; Cpsi_br[0] = tailJ; }
            }
            psi_out[0] = cs[0] * (Cpsi_out[0] - Apsi_out[0]) + Apsi_out[0];
            for (Size i = 1; i < N; ++i) {
                Apsi_out[i] = psi_out[i-1] * p[i];
                Real calBound = (psi_out[i-1] / theta[i-1]) * theta[i];
                auto [tailMin, tailJ] = forwardTailMinIdx(i);
                Real Cval = calBound; int branch = -1;
                if (f[i] < Cval) { Cval = f[i]; branch = static_cast<int>(i); }
                if (tailMin < Cval) { Cval = tailMin; branch = tailJ; }
                Cpsi_out[i] = Cval; Cpsi_br[i] = branch;
                psi_out[i] = cs[i] * (Cpsi_out[i] - Apsi_out[i]) + Apsi_out[i];
            }
            return {psi_out, Apsi_out, Cpsi_out, Cpsi_br};
        };

        auto [psi_lo, Apsi_lo, Cpsi_lo, Cpsi_br_lo] =
            replayChain(std::vector<Real>(gp.cs_lo.begin(), gp.cs_lo.end()));
        auto [psi_hi, Apsi_hi, Cpsi_hi, Cpsi_br_hi] =
            replayChain(std::vector<Real>(gp.cs_hi.begin(), gp.cs_hi.end()));

        // Global param layout: [rho(N), theta1(1), a(N-1), cs_lo(N), cs_hi(N)]
        Size rhoOff = 0, th1Off = N, aOff = N + 1;
        Size cLoOff = 2 * N, cHiOff = 3 * N;

        // Output: [dTheta(N×P), dPsiLo(N×P), dPsiHi(N×P)]
        std::vector<Real> result(3 * N * P, 0.0);
        Real* dTheta = result.data();
        Real* dPsiLo = result.data() + N * P;
        Real* dPsiHi = result.data() + 2 * N * P;

        // Helper to propagate one c-chain's psi tangents
        auto propagatePsiChain = [&](const std::vector<Real>& cs,
                                      const std::vector<Real>& psi_v,
                                      const std::vector<Real>& Apsi_v,
                                      const std::vector<Real>& Cpsi_v,
                                      const std::vector<int>& Cpsi_br_v,
                                      Size cOff_v,
                                      Real* dPsi_out,
                                      Size gi,
                                      const std::vector<Real>& dTh,
                                      auto& dp_dg_fn, auto& df_dg_fn) {
            Real dPsi_prev = 0.0;
            for (Size s = 0; s < N; ++s) {
                Real dApsi = (s == 0) ? 0.0
                    : dPsi_prev * p[s] + psi_v[s-1] * dp_dg_fn(s);

                Real dCpsi = 0.0;
                if (Cpsi_br_v[s] == -1 && s > 0) {
                    dCpsi = (dPsi_prev * theta[s] + psi_v[s-1] * dTh[s]
                             - psi_v[s-1] * theta[s] / theta[s-1] * dTh[s-1])
                            / theta[s-1];
                } else if (Cpsi_br_v[s] >= 0) {
                    int j = Cpsi_br_v[s];
                    if (j == static_cast<int>(s)) {
                        dCpsi = df_dg_fn(static_cast<Size>(j));
                    } else {
                        Real prodP_v = 1.0;
                        for (Size m = s + 1; m <= static_cast<Size>(j); ++m)
                            prodP_v *= p[m];
                        Real dfj = df_dg_fn(static_cast<Size>(j));
                        Real dProdP = 0.0;
                        for (Size m = s + 1; m <= static_cast<Size>(j); ++m)
                            dProdP += dp_dg_fn(m) * prodP_v / p[m];
                        dCpsi = dfj / prodP_v - f[j] * dProdP / (prodP_v * prodP_v);
                    }
                }

                Real dPsi_s = cs[s] * dCpsi + (1.0 - cs[s]) * dApsi;
                if (gi == cOff_v + s) dPsi_s += (Cpsi_v[s] - Apsi_v[s]);

                dPsi_out[s * P + gi] = dPsi_s;
                dPsi_prev = dPsi_s;
            }
        };

        std::vector<Real> cs_lo_v(gp.cs_lo.begin(), gp.cs_lo.end());
        std::vector<Real> cs_hi_v(gp.cs_hi.begin(), gp.cs_hi.end());

        for (Size gi = 0; gi < P; ++gi) {
            auto dp_dg = [&](Size i) -> Real {
                if (i == 0) return 0.0;
                Real dp_drho_prev, dp_drho_curr;
                if (p_branch[i] == 0) {
                    dp_drho_prev = 1.0 / (1.0 + gp.rhos[i]);
                    dp_drho_curr = -(1.0 + gp.rhos[i-1]) / ((1.0 + gp.rhos[i]) * (1.0 + gp.rhos[i]));
                } else {
                    dp_drho_prev = -1.0 / (1.0 - gp.rhos[i]);
                    dp_drho_curr = (1.0 - gp.rhos[i-1]) / ((1.0 - gp.rhos[i]) * (1.0 - gp.rhos[i]));
                }
                Real r = 0.0;
                if (gi == rhoOff + i - 1) r += dp_drho_prev;
                if (gi == rhoOff + i)     r += dp_drho_curr;
                return r;
            };

            // Theta tangents
            std::vector<Real> dTh(N, 0.0);
            dTh[0] = (gi == th1Off) ? 1.0 : 0.0;
            for (Size j = 1; j < N; ++j) {
                dTh[j] = dTh[j-1] * p[j] + theta[j-1] * dp_dg(j);
                if (gi == aOff + (j - 1)) dTh[j] += 1.0;
            }
            for (Size s = 0; s < N; ++s)
                dTheta[s * P + gi] = dTh[s];

            auto df_dg = [&](Size j) -> Real {
                Real absRho = std::abs(gp.rhos[j]);
                Real signRho = (gp.rhos[j] >= 0.0) ? 1.0 : -1.0;
                if (f_branch[j] == 0) {
                    if (gi == rhoOff + j)
                        return -4.0 * signRho / ((1.0 + absRho) * (1.0 + absRho));
                    return 0.0;
                } else {
                    Real df_dt = 2.0 / (f[j] * (1.0 + absRho));
                    Real df_dr = -2.0 * theta[j] * signRho
                                 / (f[j] * (1.0 + absRho) * (1.0 + absRho));
                    Real r = df_dt * dTh[j];
                    if (gi == rhoOff + j) r += df_dr;
                    return r;
                }
            };

            // Skip hi c-params for lo chain and vice versa
            bool skip_lo = (gi >= cHiOff && gi < cHiOff + N);
            bool skip_hi = (gi >= cLoOff && gi < cLoOff + N);

            if (!skip_lo)
                propagatePsiChain(cs_lo_v, psi_lo, Apsi_lo, Cpsi_lo, Cpsi_br_lo,
                                  cLoOff, dPsiLo, gi, dTh, dp_dg, df_dg);
            if (!skip_hi)
                propagatePsiChain(cs_hi_v, psi_hi, Apsi_hi, Cpsi_hi, Cpsi_br_hi,
                                  cHiOff, dPsiHi, gi, dTh, dp_dg, df_dg);
        }

        return result;
    }

    // =================================================================
    // DualWingEssviSurface: global gradient (forward-mode AD)
    // =================================================================

    std::vector<Real> DualWingEssviSurface::impliedVolGlobalGradient(
            Size idx, Real k,
            const DualWingEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {

        const Size N = gp.numSlices();
        QL_REQUIRE(N == slices_.size(), "global params N != slices count");
        QL_REQUIRE(idx < N, "sliceIdx >= N");

        // Determine which wing is active
        bool loWing = (k < 0.0);

        // ── Replay globalToSlice recording intermediates ──
        std::vector<Real> theta(N), p(N, 0.0), f(N);
        std::vector<int> p_branch(N, -1), f_branch(N, -1);

        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            Real a = (1.0 + gp.rhos[i-1]) / (1.0 + gp.rhos[i]);
            Real b = (1.0 - gp.rhos[i-1]) / (1.0 - gp.rhos[i]);
            p[i] = std::max(a, b);
            p_branch[i] = (a >= b) ? 0 : 1;
            theta[i] = theta[i-1] * p[i] + gp.as[i-1];
        }

        for (Size i = 0; i < N; ++i) {
            Real absRho = std::abs(gp.rhos[i]);
            Real lee = essviPsiUpperLee(absRho);
            Real fBfly = (bflyCond == EssviButterflyCondition::GatheralJacquier)
                ? essviButterflyBoundGJ(theta[i], absRho)
                : essviButterflyBoundMM(theta[i], absRho);
            Real sqrtF = std::sqrt(fBfly);
            if (lee <= sqrtF) { f[i] = lee; f_branch[i] = 0; }
            else              { f[i] = sqrtF; f_branch[i] = 1; }
        }

        auto forwardTailMinIdx = [&](Size i) -> std::pair<Real, int> {
            Real result = std::numeric_limits<Real>::max();
            int bestJ = -1;
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= p[j];
                Real val = f[j] / prodP;
                if (val < result) { result = val; bestJ = static_cast<int>(j); }
            }
            return {result, bestJ};
        };

        // Replay the active wing's c-chain
        const std::vector<Real>& cs = loWing ? gp.cs_lo : gp.cs_hi;

        std::vector<Real> psi(N), Apsi(N), Cpsi(N);
        std::vector<int> Cpsi_branch(N, -1);

        Apsi[0] = 0.0;
        {
            auto [tailMin, tailJ] = forwardTailMinIdx(0);
            if (f[0] <= tailMin) { Cpsi[0] = f[0]; Cpsi_branch[0] = 0; }
            else                 { Cpsi[0] = tailMin; Cpsi_branch[0] = tailJ; }
        }
        psi[0] = cs[0] * (Cpsi[0] - Apsi[0]) + Apsi[0];

        for (Size i = 1; i < N; ++i) {
            Apsi[i] = psi[i-1] * p[i];
            Real calBound = (psi[i-1] / theta[i-1]) * theta[i];
            auto [tailMin, tailJ] = forwardTailMinIdx(i);
            Real Cval = calBound; int branch = -1;
            if (f[i] < Cval) { Cval = f[i]; branch = static_cast<int>(i); }
            if (tailMin < Cval) { Cval = tailMin; branch = tailJ; }
            Cpsi[i] = Cval; Cpsi_branch[i] = branch;
            psi[i] = cs[i] * (Cpsi[i] - Apsi[i]) + Apsi[i];
        }

        // ── Per-slice gradient ──
        DualWingEssviSliceGradient sliceGrad =
            dualWingEssviImpliedVolGradient(slices_.at(idx), k, T_.at(idx));

        // ── Forward-mode AD through global params ──
        // Layout: [rho(N), theta1(1), a(N-1), cs_lo(N), cs_hi(N)]
        Size nParams = 4 * N;
        Size rhoOff = 0, th1Off = N, aOff = N + 1;
        Size cLoOff = 2 * N, cHiOff = 3 * N;
        Size cOff = loWing ? cLoOff : cHiOff;

        std::vector<Real> grad(nParams, 0.0);

        // Direct rho contribution from the eSSVI formula
        grad[rhoOff + idx] += sliceGrad.dSigma_dRho;

        // Active psi sensitivity
        Real dSigma_dPsi = loWing ? sliceGrad.dSigma_dPsiLo : sliceGrad.dSigma_dPsiHi;

        for (Size gi = 0; gi < nParams; ++gi) {
            // Skip inactive wing's c params — they don't affect this observation
            if (loWing && gi >= cHiOff && gi < cHiOff + N) continue;
            if (!loWing && gi >= cLoOff && gi < cLoOff + N) continue;

            auto dp_dg = [&](Size i) -> Real {
                if (i == 0) return 0.0;
                Real dp_drho_prev, dp_drho_curr;
                if (p_branch[i] == 0) {
                    dp_drho_prev = 1.0 / (1.0 + gp.rhos[i]);
                    dp_drho_curr = -(1.0 + gp.rhos[i-1]) / ((1.0 + gp.rhos[i]) * (1.0 + gp.rhos[i]));
                } else {
                    dp_drho_prev = -1.0 / (1.0 - gp.rhos[i]);
                    dp_drho_curr = (1.0 - gp.rhos[i-1]) / ((1.0 - gp.rhos[i]) * (1.0 - gp.rhos[i]));
                }
                Real result = 0.0;
                if (gi == rhoOff + i - 1) result += dp_drho_prev;
                if (gi == rhoOff + i)     result += dp_drho_curr;
                return result;
            };

            std::vector<Real> dTheta_all(N, 0.0);
            dTheta_all[0] = (gi == th1Off) ? 1.0 : 0.0;
            for (Size j = 1; j < N; ++j) {
                dTheta_all[j] = dTheta_all[j-1] * p[j] + theta[j-1] * dp_dg(j);
                if (gi == aOff + (j - 1))
                    dTheta_all[j] += 1.0;
            }

            auto df_dg = [&](Size j) -> Real {
                Real absRho = std::abs(gp.rhos[j]);
                Real signRho = (gp.rhos[j] >= 0.0) ? 1.0 : -1.0;
                if (f_branch[j] == 0) {
                    if (gi == rhoOff + j)
                        return -4.0 * signRho / ((1.0 + absRho) * (1.0 + absRho));
                    return 0.0;
                } else {
                    Real df_dtheta = 2.0 / (f[j] * (1.0 + absRho));
                    Real df_drho = -2.0 * theta[j] * signRho
                                   / (f[j] * (1.0 + absRho) * (1.0 + absRho));
                    Real result = df_dtheta * dTheta_all[j];
                    if (gi == rhoOff + j) result += df_drho;
                    return result;
                }
            };

            Real dPsi_prev = 0.0;
            for (Size i = 0; i <= idx; ++i) {
                Real dApsi_i = (i == 0) ? 0.0
                    : dPsi_prev * p[i] + psi[i-1] * dp_dg(i);

                Real dCpsi_i = 0.0;
                if (Cpsi_branch[i] == -1 && i > 0) {
                    dCpsi_i = (dPsi_prev * theta[i] + psi[i-1] * dTheta_all[i]
                               - psi[i-1] * theta[i] / theta[i-1] * dTheta_all[i-1])
                              / theta[i-1];
                } else if (Cpsi_branch[i] >= 0) {
                    int j = Cpsi_branch[i];
                    if (j == static_cast<int>(i)) {
                        dCpsi_i = df_dg(static_cast<Size>(j));
                    } else {
                        Real prodP = 1.0;
                        for (Size m = i + 1; m <= static_cast<Size>(j); ++m)
                            prodP *= p[m];
                        Real dfj = df_dg(static_cast<Size>(j));
                        Real dProdP = 0.0;
                        for (Size m = i + 1; m <= static_cast<Size>(j); ++m)
                            dProdP += dp_dg(m) * prodP / p[m];
                        dCpsi_i = dfj / prodP - f[j] * dProdP / (prodP * prodP);
                    }
                }

                Real dPsi_i = cs[i] * dCpsi_i + (1.0 - cs[i]) * dApsi_i;
                if (gi == cOff + i)
                    dPsi_i += (Cpsi[i] - Apsi[i]);

                dPsi_prev = dPsi_i;
            }

            grad[gi] += sliceGrad.dSigma_dTheta * dTheta_all[idx]
                      + dSigma_dPsi * dPsi_prev;
        }

        return grad;
    }

    // =================================================================
    // SplitRhoEssviSurface: globalToSlice
    // =================================================================

    std::vector<SplitRhoEssviSliceParams>
    SplitRhoEssviSurface::globalToSlice(const SplitRhoEssviGlobalParams& gp,
                                          EssviButterflyCondition::Type cond) {
        const Size N = gp.numSlices();
        QL_REQUIRE(N > 0, "SplitRhoEssviSurface: empty parameter set");
        QL_REQUIRE(gp.rhos_hi.size() == N, "rhos_hi size mismatch");
        QL_REQUIRE(gp.as.size() == N - 1, "as size mismatch");
        QL_REQUIRE(gp.cs_lo.size() == N && gp.cs_hi.size() == N, "cs size mismatch");
        QL_REQUIRE(gp.theta1 > 0.0, "theta1 must be > 0");

        for (Size i = 0; i < N; ++i) {
            QL_REQUIRE(gp.rhos_lo[i] > -1.0 && gp.rhos_lo[i] < 1.0,
                       "rhos_lo[" << i << "] out of range");
            QL_REQUIRE(gp.rhos_hi[i] > -1.0 && gp.rhos_hi[i] < 1.0,
                       "rhos_hi[" << i << "] out of range");
            QL_REQUIRE(gp.cs_lo[i] > 0.0 && gp.cs_lo[i] < 1.0,
                       "cs_lo[" << i << "] out of range");
            QL_REQUIRE(gp.cs_hi[i] > 0.0 && gp.cs_hi[i] < 1.0,
                       "cs_hi[" << i << "] out of range");
        }
        for (Size i = 0; i < N - 1; ++i)
            QL_REQUIRE(gp.as[i] > 0.0, "a[" << i << "] must be > 0");

        // Per-wing p values and combined p for theta chain
        std::vector<Real> p_lo(N, 0.0), p_hi(N, 0.0), p(N, 0.0);
        std::vector<Real> theta(N), f_lo(N), f_hi(N);

        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            p_lo[i] = essviCalendarP(gp.rhos_lo[i - 1], gp.rhos_lo[i]);
            p_hi[i] = essviCalendarP(gp.rhos_hi[i - 1], gp.rhos_hi[i]);
            p[i] = std::max(p_lo[i], p_hi[i]);
            theta[i] = theta[i - 1] * p[i] + gp.as[i - 1];
        }

        for (Size i = 0; i < N; ++i) {
            f_lo[i] = essviSliceButterflyBound(theta[i], std::abs(gp.rhos_lo[i]), cond);
            f_hi[i] = essviSliceButterflyBound(theta[i], std::abs(gp.rhos_hi[i]), cond);
        }

        // Run independent c-chain for each wing
        auto runChain = [&](const std::vector<Real>& rhos,
                            const std::vector<Real>& cs,
                            const std::vector<Real>& pw,
                            const std::vector<Real>& fw) {
            std::vector<Real> psi(N);

            auto fwdTail = [&](Size i) -> Real {
                Real result = std::numeric_limits<Real>::max();
                Real prodP = 1.0;
                for (Size j = i + 1; j < N; ++j) {
                    prodP *= pw[j];
                    result = std::min(result, fw[j] / prodP);
                }
                return result;
            };

            Real Apsi = 0.0;
            Real Cpsi = std::min(fw[0], fwdTail(0));
            psi[0] = cs[0] * (Cpsi - Apsi) + Apsi;

            for (Size i = 1; i < N; ++i) {
                Apsi = psi[i - 1] * pw[i];
                Real calBound = (psi[i - 1] / theta[i - 1]) * theta[i];
                Real tail = fwdTail(i);
                Cpsi = std::min({calBound, fw[i], tail});
                QL_REQUIRE(Cpsi > Apsi,
                           "SplitRho: infeasible at slice " << i);
                psi[i] = cs[i] * (Cpsi - Apsi) + Apsi;
            }
            return psi;
        };

        std::vector<Real> psi_lo = runChain(gp.rhos_lo, gp.cs_lo, p_lo, f_lo);
        std::vector<Real> psi_hi = runChain(gp.rhos_hi, gp.cs_hi, p_hi, f_hi);

        std::vector<SplitRhoEssviSliceParams> out(N);
        for (Size i = 0; i < N; ++i) {
            out[i].theta  = theta[i];
            out[i].rho_lo = gp.rhos_lo[i];
            out[i].psi_lo = psi_lo[i];
            out[i].rho_hi = gp.rhos_hi[i];
            out[i].psi_hi = psi_hi[i];
        }
        return out;
    }

    // =================================================================
    // SplitRhoEssviSurface: constructors
    // =================================================================

    SplitRhoEssviSurface::SplitRhoEssviSurface(
            const std::vector<Real>& maturities,
            const SplitRhoEssviGlobalParams& gparams,
            EssviButterflyCondition::Type bflyCond)
        : T_(maturities),
          slices_(globalToSlice(gparams, bflyCond)),
          bflyCond_(bflyCond)
    {
        QL_REQUIRE(T_.size() == slices_.size(), "maturity/slice count mismatch");
        for (Size i = 0; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > 0.0, "maturity[" << i << "] must be > 0");
        for (Size i = 1; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > T_[i - 1], "maturities must be strictly increasing");
    }

    SplitRhoEssviSurface::SplitRhoEssviSurface(
            const std::vector<Real>& maturities,
            const std::vector<SplitRhoEssviSliceParams>& slices)
        : T_(maturities), slices_(slices),
          bflyCond_(EssviButterflyCondition::GatheralJacquier)
    {
        QL_REQUIRE(T_.size() == slices_.size(), "maturity/slice count mismatch");
        QL_REQUIRE(!T_.empty(), "maturities must not be empty");
        for (Size i = 0; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > 0.0, "maturity[" << i << "] must be > 0");
        for (Size i = 1; i < T_.size(); ++i)
            QL_REQUIRE(T_[i] > T_[i - 1], "maturities must be strictly increasing");
    }

    // =================================================================
    // SplitRhoEssviSurface: evaluation
    // =================================================================

    Real SplitRhoEssviSurface::totalVariance(Size idx, Real k) const {
        return splitRhoEssviTotalVariance(slices_.at(idx), k);
    }

    SplitRhoEssviSliceParams SplitRhoEssviSurface::interpolate(Real t) const {
        const Size N = T_.size();
        if (t <= T_[0]) {
            Real lam = std::max(t / T_[0], 0.0);
            return { lam * slices_[0].theta,
                     slices_[0].rho_lo, lam * slices_[0].psi_lo,
                     slices_[0].rho_hi, lam * slices_[0].psi_hi };
        }
        if (t >= T_[N - 1]) {
            Real slope = (N >= 2)
                ? (slices_[N-1].theta - slices_[N-2].theta) / (T_[N-1] - T_[N-2])
                : slices_[0].theta / T_[0];
            return { slices_[N-1].theta + slope * (t - T_[N-1]),
                     slices_[N-1].rho_lo, slices_[N-1].psi_lo,
                     slices_[N-1].rho_hi, slices_[N-1].psi_hi };
        }
        auto it = std::upper_bound(T_.begin(), T_.end(), t);
        Size hi = static_cast<Size>(it - T_.begin());
        Size lo = hi - 1;
        Real lam = (t - T_[lo]) / (T_[hi] - T_[lo]);

        const auto& s1 = slices_[lo];
        const auto& s2 = slices_[hi];

        SplitRhoEssviSliceParams st;
        st.theta  = (1.0 - lam) * s1.theta  + lam * s2.theta;
        st.psi_lo = (1.0 - lam) * s1.psi_lo + lam * s2.psi_lo;
        st.psi_hi = (1.0 - lam) * s1.psi_hi + lam * s2.psi_hi;
        // Blend rho per wing via psi-weighted average
        Real psiRhoLo = (1.0 - lam) * s1.psi_lo * s1.rho_lo + lam * s2.psi_lo * s2.rho_lo;
        st.rho_lo = (st.psi_lo > 1e-15) ? psiRhoLo / st.psi_lo : 0.0;
        Real psiRhoHi = (1.0 - lam) * s1.psi_hi * s1.rho_hi + lam * s2.psi_hi * s2.rho_hi;
        st.rho_hi = (st.psi_hi > 1e-15) ? psiRhoHi / st.psi_hi : 0.0;
        return st;
    }

    Real SplitRhoEssviSurface::totalVariance(Real k, Real t) const {
        return splitRhoEssviTotalVariance(interpolate(t), k);
    }

    Real SplitRhoEssviSurface::impliedVol(Real k, Real t) const {
        QL_REQUIRE(t > 0.0, "impliedVol: maturity must be > 0");
        Real w = totalVariance(k, t);
        QL_REQUIRE(w >= 0.0, "impliedVol: negative total variance " << w);
        return std::sqrt(w / t);
    }

    SplitRhoEssviSliceGradient
    SplitRhoEssviSurface::impliedVolGradient(Size idx, Real k) const {
        return splitRhoEssviImpliedVolGradient(slices_.at(idx), k, T_.at(idx));
    }

    // =================================================================
    // SplitRhoEssviSurface: global gradient (forward-mode AD)
    // =================================================================

    std::vector<Real> SplitRhoEssviSurface::impliedVolGlobalGradient(
            Size idx, Real k,
            const SplitRhoEssviGlobalParams& gp,
            EssviButterflyCondition::Type bflyCond) const {

        const Size N = gp.numSlices();
        QL_REQUIRE(N == slices_.size(), "global params N != slices count");
        QL_REQUIRE(idx < N, "sliceIdx >= N");

        bool loWing = (k < 0.0);

        // Replay globalToSlice with intermediates
        std::vector<Real> theta(N), p_lo(N, 0.0), p_hi(N, 0.0), p(N, 0.0);
        std::vector<Real> f_w(N);  // butterfly bound for active wing
        std::vector<int> p_lo_branch(N, -1), p_hi_branch(N, -1), p_max_branch(N, -1);
        std::vector<int> f_branch(N, -1);

        const std::vector<Real>& rhos = loWing ? gp.rhos_lo : gp.rhos_hi;
        const std::vector<Real>& cs   = loWing ? gp.cs_lo   : gp.cs_hi;

        theta[0] = gp.theta1;
        for (Size i = 1; i < N; ++i) {
            Real a_lo = (1.0 + gp.rhos_lo[i-1]) / (1.0 + gp.rhos_lo[i]);
            Real b_lo = (1.0 - gp.rhos_lo[i-1]) / (1.0 - gp.rhos_lo[i]);
            p_lo[i] = std::max(a_lo, b_lo);
            p_lo_branch[i] = (a_lo >= b_lo) ? 0 : 1;

            Real a_hi = (1.0 + gp.rhos_hi[i-1]) / (1.0 + gp.rhos_hi[i]);
            Real b_hi = (1.0 - gp.rhos_hi[i-1]) / (1.0 - gp.rhos_hi[i]);
            p_hi[i] = std::max(a_hi, b_hi);
            p_hi_branch[i] = (a_hi >= b_hi) ? 0 : 1;

            p[i] = std::max(p_lo[i], p_hi[i]);
            p_max_branch[i] = (p_lo[i] >= p_hi[i]) ? 0 : 1;
            theta[i] = theta[i - 1] * p[i] + gp.as[i - 1];
        }

        // Active wing's butterfly bounds and p
        const std::vector<Real>& pw = loWing ? p_lo : p_hi;

        for (Size i = 0; i < N; ++i) {
            Real absRho = std::abs(rhos[i]);
            Real lee = essviPsiUpperLee(absRho);
            Real fBfly = (bflyCond == EssviButterflyCondition::GatheralJacquier)
                ? essviButterflyBoundGJ(theta[i], absRho)
                : essviButterflyBoundMM(theta[i], absRho);
            Real sqrtF = std::sqrt(fBfly);
            if (lee <= sqrtF) { f_w[i] = lee; f_branch[i] = 0; }
            else              { f_w[i] = sqrtF; f_branch[i] = 1; }
        }

        // Replay active wing c-chain
        auto fwdTailMinIdx = [&](Size i) -> std::pair<Real, int> {
            Real result = std::numeric_limits<Real>::max();
            int bestJ = -1;
            Real prodP = 1.0;
            for (Size j = i + 1; j < N; ++j) {
                prodP *= pw[j];
                Real val = f_w[j] / prodP;
                if (val < result) { result = val; bestJ = static_cast<int>(j); }
            }
            return {result, bestJ};
        };

        std::vector<Real> psi(N), Apsi(N), Cpsi(N);
        std::vector<int> Cpsi_branch(N, -1);

        Apsi[0] = 0.0;
        {
            auto [tailMin, tailJ] = fwdTailMinIdx(0);
            if (f_w[0] <= tailMin) { Cpsi[0] = f_w[0]; Cpsi_branch[0] = 0; }
            else                   { Cpsi[0] = tailMin; Cpsi_branch[0] = tailJ; }
        }
        psi[0] = cs[0] * (Cpsi[0] - Apsi[0]) + Apsi[0];

        for (Size i = 1; i < N; ++i) {
            Apsi[i] = psi[i - 1] * pw[i];
            Real calBound = (psi[i - 1] / theta[i - 1]) * theta[i];
            auto [tailMin, tailJ] = fwdTailMinIdx(i);
            Real Cval = calBound; int branch = -1;
            if (f_w[i] < Cval) { Cval = f_w[i]; branch = static_cast<int>(i); }
            if (tailMin < Cval) { Cval = tailMin; branch = tailJ; }
            Cpsi[i] = Cval; Cpsi_branch[i] = branch;
            psi[i] = cs[i] * (Cpsi[i] - Apsi[i]) + Apsi[i];
        }

        // Per-slice gradient
        SplitRhoEssviSliceGradient sliceGrad =
            splitRhoEssviImpliedVolGradient(slices_.at(idx), k, T_.at(idx));

        // Global param layout: [rhos_lo(N), rhos_hi(N), theta1(1), as(N-1), cs_lo(N), cs_hi(N)]
        Size nParams = 5 * N;
        Size rLoOff = 0, rHiOff = N, th1Off = 2*N, aOff = 2*N + 1;
        Size cLoOff = 3*N, cHiOff = 4*N;

        // Active wing offsets
        Size rOff = loWing ? rLoOff : rHiOff;
        Size cOff = loWing ? cLoOff : cHiOff;

        std::vector<Real> grad(nParams, 0.0);

        // Direct rho contribution from the eSSVI formula
        Real dSigma_dRho = loWing ? sliceGrad.dSigma_dRhoLo : sliceGrad.dSigma_dRhoHi;
        Real dSigma_dPsi = loWing ? sliceGrad.dSigma_dPsiLo : sliceGrad.dSigma_dPsiHi;
        grad[rOff + idx] += dSigma_dRho;

        // Forward-mode AD: for each global param, propagate tangent
        for (Size gi = 0; gi < nParams; ++gi) {
            // Skip inactive wing's c params and rho params that don't affect theta
            bool isInactiveC = loWing ? (gi >= cHiOff && gi < cHiOff + N)
                                      : (gi >= cLoOff && gi < cLoOff + N);
            if (isInactiveC) continue;

            // Inactive wing's rho only affects theta through p_max
            // (it affects p_hi or p_lo which may determine p[i] via max)

            // dp_lo/dg and dp_hi/dg
            auto dp_wing_dg = [&](Size i, const std::vector<Real>& wing_rhos,
                                   const std::vector<int>& wing_pb, Size wing_rOff) -> Real {
                if (i == 0) return 0.0;
                Real dp_dr_prev, dp_dr_curr;
                if (wing_pb[i] == 0) {
                    dp_dr_prev = 1.0 / (1.0 + wing_rhos[i]);
                    dp_dr_curr = -(1.0 + wing_rhos[i-1]) / ((1.0 + wing_rhos[i]) * (1.0 + wing_rhos[i]));
                } else {
                    dp_dr_prev = -1.0 / (1.0 - wing_rhos[i]);
                    dp_dr_curr = (1.0 - wing_rhos[i-1]) / ((1.0 - wing_rhos[i]) * (1.0 - wing_rhos[i]));
                }
                Real result = 0.0;
                if (gi == wing_rOff + i - 1) result += dp_dr_prev;
                if (gi == wing_rOff + i)     result += dp_dr_curr;
                return result;
            };

            // dp[i]/dg = dp_lo/dg or dp_hi/dg depending on which is max
            auto dp_dg = [&](Size i) -> Real {
                if (i == 0) return 0.0;
                if (p_max_branch[i] == 0)
                    return dp_wing_dg(i, gp.rhos_lo, p_lo_branch, rLoOff);
                else
                    return dp_wing_dg(i, gp.rhos_hi, p_hi_branch, rHiOff);
            };

            // Active wing's p derivative (for the c-chain)
            auto dpw_dg = [&](Size i) -> Real {
                if (i == 0) return 0.0;
                if (loWing)
                    return dp_wing_dg(i, gp.rhos_lo, p_lo_branch, rLoOff);
                else
                    return dp_wing_dg(i, gp.rhos_hi, p_hi_branch, rHiOff);
            };

            // Theta tangents (theta chain uses combined p)
            std::vector<Real> dTheta(N, 0.0);
            dTheta[0] = (gi == th1Off) ? 1.0 : 0.0;
            for (Size j = 1; j < N; ++j) {
                dTheta[j] = dTheta[j-1] * p[j] + theta[j-1] * dp_dg(j);
                if (gi == aOff + (j - 1))
                    dTheta[j] += 1.0;
            }

            // Active wing's f derivative
            auto df_dg = [&](Size j) -> Real {
                Real absRho = std::abs(rhos[j]);
                Real signRho = (rhos[j] >= 0.0) ? 1.0 : -1.0;
                if (f_branch[j] == 0) {
                    if (gi == rOff + j)
                        return -4.0 * signRho / ((1.0 + absRho) * (1.0 + absRho));
                    return 0.0;
                } else {
                    Real df_dt = 2.0 / (f_w[j] * (1.0 + absRho));
                    Real df_dr = -2.0 * theta[j] * signRho
                                 / (f_w[j] * (1.0 + absRho) * (1.0 + absRho));
                    Real result = df_dt * dTheta[j];
                    if (gi == rOff + j) result += df_dr;
                    return result;
                }
            };

            // Propagate through c-chain
            Real dPsi_prev = 0.0;
            for (Size i = 0; i <= idx; ++i) {
                Real dApsi = (i == 0) ? 0.0
                    : dPsi_prev * pw[i] + psi[i-1] * dpw_dg(i);

                Real dCpsi = 0.0;
                if (Cpsi_branch[i] == -1 && i > 0) {
                    dCpsi = (dPsi_prev * theta[i] + psi[i-1] * dTheta[i]
                             - psi[i-1] * theta[i] / theta[i-1] * dTheta[i-1])
                            / theta[i-1];
                } else if (Cpsi_branch[i] >= 0) {
                    int j = Cpsi_branch[i];
                    if (j == static_cast<int>(i)) {
                        dCpsi = df_dg(static_cast<Size>(j));
                    } else {
                        Real prodP = 1.0;
                        for (Size m = i + 1; m <= static_cast<Size>(j); ++m)
                            prodP *= pw[m];
                        Real dfj = df_dg(static_cast<Size>(j));
                        Real dProdP = 0.0;
                        for (Size m = i + 1; m <= static_cast<Size>(j); ++m)
                            dProdP += dpw_dg(m) * prodP / pw[m];
                        dCpsi = dfj / prodP - f_w[j] * dProdP / (prodP * prodP);
                    }
                }

                Real dPsi_i = cs[i] * dCpsi + (1.0 - cs[i]) * dApsi;
                if (gi == cOff + i)
                    dPsi_i += (Cpsi[i] - Apsi[i]);

                dPsi_prev = dPsi_i;
            }

            grad[gi] += sliceGrad.dSigma_dTheta * dTheta[idx]
                      + dSigma_dPsi * dPsi_prev;
        }

        return grad;
    }

} // namespace QuantLib
