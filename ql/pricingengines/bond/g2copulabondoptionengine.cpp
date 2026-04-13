/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 chloride contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
*/

#include <ql/pricingengines/bond/g2copulabondoptionengine.hpp>
#include <ql/math/integrals/gaussianquadratures.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/termstructures/volatility/sabr.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace QuantLib {

    namespace {

        // ── SABR CDF builder (per-tenor) ─────────────────────────────

        class SabrCdf {
          public:
            SabrCdf(Real forward, Real atmVol, Real T,
                    Real alpha, Real beta, Real nu, Real rho, Real shift,
                    Size nStrikes, Real nStd)
            : forward_(forward), T_(T)
            {
                Real volT = atmVol * std::sqrt(T);
                Real lo = forward - nStd * volT;
                Real hi = forward + nStd * volT;

                // Strike grid
                K_.resize(nStrikes);
                for (Size i = 0; i < nStrikes; ++i)
                    K_[i] = lo + (hi - lo) * i / (nStrikes - 1);

                // Call prices across strike grid
                std::vector<Real> calls(nStrikes);
                bool useSabr = alpha > 0.0 && T > 0.0;
                Real sqrtT = std::sqrt(T);
                Real shiftedFwd = forward + shift;

                for (Size i = 0; i < nStrikes; ++i) {
                    if (useSabr) {
                        Real shiftedK = K_[i] + shift;
                        if (shiftedK > 0.0 && shiftedFwd > 0.0) {
                            Real slnVol = sabrFlochKennedyVolatility(
                                shiftedK, shiftedFwd, T, alpha, beta, nu, rho);
                            calls[i] = blackFormula(
                                Option::Call, shiftedK, shiftedFwd, slnVol * sqrtT);
                        } else {
                            calls[i] = std::max(forward - K_[i], 0.0);
                        }
                    } else {
                        // Flat Bachelier
                        if (volT < 1e-15) {
                            calls[i] = std::max(forward - K_[i], 0.0);
                        } else {
                            Real d = (forward - K_[i]) / volT;
                            CumulativeNormalDistribution Phi;
                            NormalDistribution phi;
                            calls[i] = (forward - K_[i]) * Phi(d) + volT * phi(d);
                        }
                    }
                }

                // Breeden-Litzenberger: CDF = 1 + dC/dK
                Size nMid = nStrikes - 1;
                Kmid_.resize(nMid);
                cdf_.resize(nMid);

                for (Size i = 0; i < nMid; ++i) {
                    Real dK = K_[i + 1] - K_[i];
                    Real dC = calls[i + 1] - calls[i];
                    Kmid_[i] = 0.5 * (K_[i] + K_[i + 1]);
                    cdf_[i] = std::max(0.0, std::min(1.0, 1.0 + dC / dK));
                }

                // Enforce monotonicity
                for (Size i = 1; i < nMid; ++i)
                    cdf_[i] = std::max(cdf_[i], cdf_[i - 1]);
            }

            //! Inverse CDF: uniform → swap rate
            Real invCdf(Real u) const {
                u = std::max(1e-10, std::min(1.0 - 1e-10, u));
                // Linear interpolation on (cdf_, Kmid_)
                if (u <= cdf_.front()) return Kmid_.front();
                if (u >= cdf_.back()) return Kmid_.back();

                // Binary search
                auto it = std::lower_bound(cdf_.begin(), cdf_.end(), u);
                Size idx = std::distance(cdf_.begin(), it);
                if (idx == 0) return Kmid_[0];
                Real t = (u - cdf_[idx - 1]) / (cdf_[idx] - cdf_[idx - 1]);
                return Kmid_[idx - 1] + t * (Kmid_[idx] - Kmid_[idx - 1]);
            }

          private:
            Real forward_, T_;
            std::vector<Real> K_, Kmid_, cdf_;
        };

        // ── G2 factor loadings and Cholesky ──────────────────────────

        //! Co-initial par swap rate at G2 state (x1, x2)
        Real swapRateAt(const ext::shared_ptr<G2>& model,
                        Time T, Real tenor, Real x1, Real x2,
                        Real freq = 1.0) {
            if (tenor <= freq + 1e-10) {
                // Sub-annual: single payment
                Real df = model->discountBond(T, T + tenor, x1, x2);
                return (std::abs(df) > 1e-15) ? (1.0 - df) / (tenor * df) : 0.0;
            }
            Size nPeriods = std::max<Size>(1, (Size)std::round(tenor / freq));
            Real annuity = 0.0;
            for (Size k = 1; k <= nPeriods; ++k) {
                Real tk = T + k * freq;
                if (tk > T + tenor + 1e-10) break;
                annuity += freq * model->discountBond(T, tk, x1, x2);
            }
            Real dfFinal = model->discountBond(T, T + tenor, x1, x2);
            return (std::abs(annuity) > 1e-15) ? (1.0 - dfFinal) / annuity : 0.0;
        }

        //! Forward swap rates at (0, 0)
        std::vector<Real> forwardSwapRates(
            const ext::shared_ptr<G2>& model,
            Time T, const std::vector<Real>& tenors) {
            std::vector<Real> rates(tenors.size());
            for (Size i = 0; i < tenors.size(); ++i)
                rates[i] = swapRateAt(model, T, tenors[i], 0.0, 0.0);
            return rates;
        }

        //! Factor loadings (alpha_i, beta_i) via central differences
        void factorLoadings(
            const ext::shared_ptr<G2>& model,
            Time T, const std::vector<Real>& tenors,
            std::vector<Real>& alphas, std::vector<Real>& betas) {
            const Real eps = 1e-4;
            Size n = tenors.size();
            alphas.resize(n);
            betas.resize(n);
            for (Size i = 0; i < n; ++i) {
                Real sUp = swapRateAt(model, T, tenors[i], eps, 0.0);
                Real sDn = swapRateAt(model, T, tenors[i], -eps, 0.0);
                alphas[i] = (sUp - sDn) / (2 * eps);

                sUp = swapRateAt(model, T, tenors[i], 0.0, eps);
                sDn = swapRateAt(model, T, tenors[i], 0.0, -eps);
                betas[i] = (sUp - sDn) / (2 * eps);
            }
        }

        //! G2 state covariance at time T
        void stateCovariance(const ext::shared_ptr<G2>& model, Time T,
                             Real& var1, Real& var2, Real& cov12) {
            Real a = model->a(), sigma = model->sigma();
            Real b = model->b(), eta = model->eta();
            Real rhoVal = model->rho();
            var1 = sigma * sigma / (2 * a) * (1.0 - std::exp(-2 * a * T));
            var2 = eta * eta / (2 * b) * (1.0 - std::exp(-2 * b * T));
            cov12 = rhoVal * sigma * eta / (a + b) * (1.0 - std::exp(-(a + b) * T));
        }

        //! Rank-2 Cholesky: N×2 matrix, rows normalized to unit variance
        //! L[i] = (alpha_i, beta_i) @ chol(V), then normalized
        Matrix rank2Cholesky(
            const ext::shared_ptr<G2>& model,
            Time T, const std::vector<Real>& tenors) {
            std::vector<Real> alphas, betas;
            factorLoadings(model, T, tenors, alphas, betas);

            Real var1, var2, cov12;
            stateCovariance(model, T, var1, var2, cov12);

            // Cholesky of 2x2 covariance
            Real L00 = std::sqrt(var1);
            Real L10 = cov12 / L00;
            Real L11 = std::sqrt(std::max(var2 - L10 * L10, 0.0));

            Size n = tenors.size();
            Matrix L(n, 2);
            for (Size i = 0; i < n; ++i) {
                Real r0 = alphas[i] * L00 + betas[i] * L10;
                Real r1 = betas[i] * L11;
                Real norm = std::sqrt(r0 * r0 + r1 * r1);
                norm = std::max(norm, 1e-15);
                L[i][0] = r0 / norm;
                L[i][1] = r1 / norm;
            }
            return L;
        }

        // ── DF bootstrap from co-initial swap rates ──────────────────

        //! Bootstrap DFs from par swap rates at sparse tenors
        std::vector<Real> bootstrapDFs(
            const std::vector<Real>& tenors,
            const std::vector<Real>& swapRates) {
            Size n = tenors.size();
            std::vector<Real> dfs(n);

            // Map: tenor → DF (for interpolation of intermediate annual nodes)
            std::vector<std::pair<Real, Real>> known;

            for (Size i = 0; i < n; ++i) {
                Real tau = tenors[i];
                Real S = swapRates[i];

                if (tau <= 1.0 + 1e-10) {
                    // Sub-annual or 1Y: DF = 1 / (1 + S * tau)
                    Real denom = 1.0 + S * tau;
                    dfs[i] = (std::abs(denom) > 1e-15) ? 1.0 / denom : 1.0;
                } else {
                    // Multi-year: annuity from known DFs
                    Size nPeriods = (Size)std::round(tau);
                    Real annuity = 0.0;
                    for (Size k = 1; k < nPeriods; ++k) {
                        Real tk = (Real)k;
                        // Log-linear interpolation from known DFs
                        Real df_k = 1.0;
                        if (!known.empty()) {
                            // Find bracketing known DFs
                            if (tk <= known.front().first) {
                                df_k = known.front().second;
                            } else if (tk >= known.back().first) {
                                df_k = known.back().second;
                            } else {
                                for (Size j = 1; j < known.size(); ++j) {
                                    if (known[j].first >= tk) {
                                        Real t0 = known[j-1].first, t1 = known[j].first;
                                        Real ld0 = std::log(known[j-1].second);
                                        Real ld1 = std::log(known[j].second);
                                        Real w = (tk - t0) / (t1 - t0);
                                        df_k = std::exp(ld0 + w * (ld1 - ld0));
                                        break;
                                    }
                                }
                            }
                        }
                        annuity += df_k;
                    }
                    Real denom = 1.0 + S;
                    dfs[i] = (std::abs(denom) > 1e-15) ?
                        (1.0 - S * annuity) / denom : 1.0;
                }
                known.push_back({tau, dfs[i]});
            }
            return dfs;
        }

        // ── Bond price from DFs ──────────────────────────────────────

        //! Log-linear interpolation of DFs at cashflow dates
        Real bondPriceFromDFs(
            const std::vector<CashflowEntry>& cashflows,
            const std::vector<Real>& tenors,
            const std::vector<Real>& dfs) {
            Size nT = tenors.size();
            if (nT == 0) return 0.0;

            // Precompute log DFs
            std::vector<Real> logDFs(nT);
            for (Size i = 0; i < nT; ++i)
                logDFs[i] = std::log(std::max(dfs[i], 1e-15));

            Real price = 0.0;
            for (const auto& cf : cashflows) {
                Real t = cf.timeFromExpiry;
                // np.interp equivalent: linear on log DFs
                Real logDF;
                if (t <= tenors.front()) {
                    logDF = logDFs.front();
                } else if (t >= tenors.back()) {
                    logDF = logDFs.back();
                } else {
                    auto it = std::lower_bound(tenors.begin(), tenors.end(), t);
                    Size idx = std::distance(tenors.begin(), it);
                    if (idx == 0) idx = 1;
                    Real w = (t - tenors[idx-1]) / (tenors[idx] - tenors[idx-1]);
                    logDF = logDFs[idx-1] + w * (logDFs[idx] - logDFs[idx-1]);
                }
                price += cf.amount * std::exp(logDF);
            }
            return price;
        }

    }  // anonymous namespace


    // ── Public API ───────────────────────────────────────────────────

    Real g2CopulaBondOptionPV(
        const ext::shared_ptr<G2>& model,
        const std::vector<CashflowEntry>& cashflows,
        const std::vector<Real>& tenors,
        const std::vector<SabrTenorParams>& sabrParams,
        Real strike,
        Time expiryTime,
        Option::Type type,
        Size quadOrder,
        Size nCdfStrikes,
        Real nCdfStd)
    {
        QL_REQUIRE(tenors.size() == sabrParams.size(),
                   "tenors and sabrParams must have same size");
        QL_REQUIRE(!cashflows.empty(), "no cashflows provided");
        QL_REQUIRE(!tenors.empty(), "no tenors provided");
        QL_REQUIRE(expiryTime > 0.0, "expiryTime must be positive");

        Size nTenors = tenors.size();

        // Forward swap rates at (0, 0)
        std::vector<Real> fwdRates = forwardSwapRates(model, expiryTime, tenors);

        // Build SABR CDFs
        std::vector<SabrCdf> cdfs;
        cdfs.reserve(nTenors);
        for (Size i = 0; i < nTenors; ++i) {
            const auto& sp = sabrParams[i];
            cdfs.emplace_back(
                fwdRates[i], sp.atmNormalVol, expiryTime,
                sp.alpha, sp.beta, sp.nu, sp.rho, sp.shift,
                nCdfStrikes, nCdfStd);
        }

        // Rank-2 Cholesky
        Matrix L = rank2Cholesky(model, expiryTime, tenors);

        // Gauss-Hermite quadrature
        GaussHermiteIntegration ghq(quadOrder);
        const Array& nodes = ghq.x();
        const Array& weights = ghq.weights();

        // QL's GaussHermiteIntegration weights do NOT include exp(-x²).
        // numpy's hermgauss weights DO.  We need physicist's Hermite:
        //   ∫f(x)exp(-x²)dx ≈ Σ w_i*exp(-x_i²) * f(x_i)
        // Then convert to probabilist's (standard normal expectation):
        //   E[g(Z)] = (1/√π) Σ [w_i*exp(-x_i²)] * g(x_i*√2)
        Real sqrt2 = std::sqrt(2.0);
        Real invSqrtPi = 1.0 / std::sqrt(M_PI);

        // Precompute corrected weights: w_i * exp(-x_i²) / √π
        Array corrWeights(quadOrder);
        for (Size i = 0; i < quadOrder; ++i)
            corrWeights[i] = weights[i] * std::exp(-nodes[i] * nodes[i]) * invSqrtPi;

        CumulativeNormalDistribution Phi;

        Real putCall = (type == Option::Call) ? 1.0 : -1.0;
        Real optionValue = 0.0;

        for (Size i = 0; i < quadOrder; ++i) {
            for (Size j = 0; j < quadOrder; ++j) {
                Real w = corrWeights[i] * corrWeights[j];
                Real z0 = nodes[i] * sqrt2;
                Real z1 = nodes[j] * sqrt2;

                // Correlated standard normals: xi[k] = L[k][0]*z0 + L[k][1]*z1
                // → uniform via Phi → swap rates via inv CDF
                std::vector<Real> swapRates(nTenors);
                for (Size k = 0; k < nTenors; ++k) {
                    Real xi = L[k][0] * z0 + L[k][1] * z1;
                    Real u = Phi(xi);
                    swapRates[k] = cdfs[k].invCdf(u);
                }

                // Bootstrap DFs → bond price
                std::vector<Real> dfs = bootstrapDFs(tenors, swapRates);
                Real bondPrice = bondPriceFromDFs(cashflows, tenors, dfs);

                // Payoff
                Real payoff = std::max(putCall * (bondPrice - strike), 0.0);
                optionValue += w * payoff;
            }
        }

        return optionValue;
    }

}
