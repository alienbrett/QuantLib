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

} // namespace QuantLib
