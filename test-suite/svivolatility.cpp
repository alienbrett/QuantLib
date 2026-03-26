/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2022 Skandinaviska Enskilda Banken AB (publ)

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

#include "toplevelfixture.hpp"
#include "utilities.hpp"
#include <ql/experimental/volatility/svismilesection.hpp>
#include <ql/termstructures/volatility/equityfx/essvihelpers.hpp>
#include <ql/termstructures/volatility/equityfx/essvivoltermstructure.hpp>
#include <ql/termstructures/volatility/equityfx/essvilocalvolsurface.hpp>
#include <ql/termstructures/volatility/equityfx/localvolsurface.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <cmath>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(SviVolatilityTests)

BOOST_AUTO_TEST_CASE(testSviSmileSection) {

    BOOST_TEST_MESSAGE("Testing SviSmileSection construction...");

    Date today = Settings::instance().evaluationDate();

    // Test time based constructor
    Time tte = 11.0 / 365;
    Real forward = 123.45;
    Real a = -0.0666;
    Real b = 0.229;
    Real sigma = 0.337;
    Real rho = 0.439;
    Real m = 0.193;
    std::vector<Real> sviParameters = {a, b, sigma, rho, m};
    Real strike = forward * std::exp(m);
    ext::shared_ptr<SviSmileSection> time_section;

    BOOST_CHECK_NO_THROW(time_section =
                             ext::make_shared<SviSmileSection>(tte, forward, sviParameters));
    BOOST_CHECK_EQUAL(time_section->atmLevel(), forward);
    QL_CHECK_CLOSE(time_section->variance(strike), a + b * sigma, 1E-10);

    Date date = today + Period(11, Days);
    ext::shared_ptr<SviSmileSection> date_section;

    BOOST_CHECK_NO_THROW(date_section =
                             ext::make_shared<SviSmileSection>(date, forward, sviParameters));

    BOOST_CHECK_EQUAL(date_section->atmLevel(), forward);
    QL_CHECK_CLOSE(date_section->variance(strike), a + b * sigma, 1E-10);
}

// =====================================================================
// eSSVI per-slice gradient: analytic vs finite difference
// =====================================================================

namespace {

    // Rebuild surface with one native param bumped and return implied vol
    Real bumpedImpliedVol(const std::vector<Real>& maturities,
                          std::vector<EssviSliceParams> slices,
                          Size sliceIdx,
                          int paramIdx,   // 0=theta, 1=rho, 2=psi
                          Real bump,
                          Real k) {
        if (paramIdx == 0) slices[sliceIdx].theta += bump;
        if (paramIdx == 1) slices[sliceIdx].rho   += bump;
        if (paramIdx == 2) slices[sliceIdx].psi   += bump;
        EssviSurface surface(maturities, slices);
        return surface.impliedVol(k, maturities[sliceIdx]);
    }

    Real fdSliceGradient(const std::vector<Real>& maturities,
                         const std::vector<EssviSliceParams>& slices,
                         Size sliceIdx, int paramIdx, Real k,
                         Real eps = 1e-6) {
        Real up   = bumpedImpliedVol(maturities, slices, sliceIdx, paramIdx, +eps, k);
        Real down = bumpedImpliedVol(maturities, slices, sliceIdx, paramIdx, -eps, k);
        return (up - down) / (2.0 * eps);
    }

}

BOOST_AUTO_TEST_CASE(testEssviSliceGradient) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI per-slice analytic gradient vs finite differences...");

    std::vector<Real> maturities = {0.25, 1.0, 2.0};
    std::vector<EssviSliceParams> slices = {
        {0.01,  -0.30, 0.05},
        {0.04,  -0.25, 0.10},
        {0.10,  -0.15, 0.18},
    };

    std::vector<Real> kValues = {-0.50, -0.20, -0.05, 0.0, 0.05, 0.20, 0.50};
    Real eps = 1e-6;
    Real tol = 1e-5;

    for (Size si = 0; si < slices.size(); ++si) {
        for (Real k : kValues) {
            EssviSurface surface(maturities, slices);
            EssviSliceGradient grad = surface.impliedVolGradient(si, k);

            Real fd_theta = fdSliceGradient(maturities, slices, si, 0, k, eps);
            Real fd_rho   = fdSliceGradient(maturities, slices, si, 1, k, eps);
            Real fd_psi   = fdSliceGradient(maturities, slices, si, 2, k, eps);

            BOOST_CHECK_SMALL(grad.dSigma_dTheta - fd_theta, tol);
            BOOST_CHECK_SMALL(grad.dSigma_dRho   - fd_rho,   tol);
            BOOST_CHECK_SMALL(grad.dSigma_dPsi   - fd_psi,   tol);
        }
    }
}

BOOST_AUTO_TEST_CASE(testEssviGradientATM) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI gradient at ATM (k=0) has known structure...");

    // At k=0: w = theta, sigma = sqrt(theta/T)
    // dw/dtheta = 1, dsigma/dtheta = 1/(2*sqrt(theta*T))
    // dw/drho = 0, dw/dpsi = 0

    std::vector<Real> maturities = {0.5};
    std::vector<EssviSliceParams> slices = {{0.02, -0.30, 0.08}};
    EssviSurface surface(maturities, slices);

    EssviSliceGradient grad = surface.impliedVolGradient(0, 0.0);

    Real expected_dTheta = 1.0 / (2.0 * std::sqrt(0.02 * 0.5));
    BOOST_CHECK_SMALL(grad.dSigma_dTheta - expected_dTheta, 1e-12);
    BOOST_CHECK_SMALL(grad.dSigma_dRho, 1e-12);
    BOOST_CHECK_SMALL(grad.dSigma_dPsi, 1e-12);
}

BOOST_AUTO_TEST_CASE(testEssviGradientSymmetryRhoZero) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI gradient symmetry when rho=0...");

    // When rho=0: w(k) = 0.5*(theta + sqrt(psi^2*k^2 + theta^2))
    // w is symmetric (even) in k.
    // dw/dtheta is even in k (depends on k^2 only)
    // dw/drho = 0.5*psi*k*(1+theta/sqrtD), which is ODD in k
    // dw/dpsi = 0.5*k*(A/sqrtD) = 0.5*k*(psi*k/sqrtD) = 0.5*psi*k^2/sqrtD, EVEN in k
    // So dsigma/dtheta and dsigma/dpsi are symmetric, dsigma/drho is antisymmetric.

    std::vector<Real> maturities = {1.0};
    std::vector<EssviSliceParams> slices = {{0.04, 0.0, 0.10}};
    EssviSurface surface(maturities, slices);

    std::vector<Real> kValues = {0.05, 0.10, 0.30};
    Real tol = 1e-12;

    for (Real k : kValues) {
        EssviSliceGradient gp = surface.impliedVolGradient(0, +k);
        EssviSliceGradient gm = surface.impliedVolGradient(0, -k);

        // theta and psi gradients symmetric (even)
        BOOST_CHECK_SMALL(gp.dSigma_dTheta - gm.dSigma_dTheta, tol);
        BOOST_CHECK_SMALL(gp.dSigma_dPsi   - gm.dSigma_dPsi,   tol);
        // rho gradient antisymmetric (odd)
        BOOST_CHECK_SMALL(gp.dSigma_dRho   + gm.dSigma_dRho,   tol);
    }
}

BOOST_AUTO_TEST_CASE(testEssviGradientExtremeParams) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI gradient with extreme parameters...");

    std::vector<Real> maturities = {0.1, 2.0};
    std::vector<EssviSliceParams> slices = {
        {0.002, -0.80, 0.01},
        {0.10,  -0.10, 0.20},
    };

    std::vector<Real> kValues = {-0.20, -0.05, 0.0, 0.05, 0.20};
    Real eps = 1e-6;
    Real tol = 1e-4;

    for (Size si = 0; si < slices.size(); ++si) {
        for (Real k : kValues) {
            EssviSurface surface(maturities, slices);
            EssviSliceGradient grad = surface.impliedVolGradient(si, k);

            Real fd_theta = fdSliceGradient(maturities, slices, si, 0, k, eps);
            Real fd_rho   = fdSliceGradient(maturities, slices, si, 1, k, eps);
            Real fd_psi   = fdSliceGradient(maturities, slices, si, 2, k, eps);

            BOOST_CHECK_SMALL(grad.dSigma_dTheta - fd_theta, tol);
            BOOST_CHECK_SMALL(grad.dSigma_dRho   - fd_rho,   tol);
            BOOST_CHECK_SMALL(grad.dSigma_dPsi   - fd_psi,   tol);
        }
    }
}

// =====================================================================
// eSSVI global gradient: analytic vs finite difference
// =====================================================================

namespace {

    // Build surface from global params with one param bumped
    Real globalBumpedVol(const std::vector<Real>& maturities,
                         EssviGlobalParams gp,
                         Size globalIdx,
                         Real bump,
                         Size sliceIdx,
                         Real k,
                         EssviButterflyCondition::Type cond) {
        Size N = gp.numSlices();
        // Layout: [rho_0..rho_{N-1}, theta1, a_0..a_{N-2}, c_0..c_{N-1}]
        if (globalIdx < N) {
            gp.rhos[globalIdx] += bump;
        } else if (globalIdx == N) {
            gp.theta1 += bump;
        } else if (globalIdx < 2 * N) {
            gp.as[globalIdx - N - 1] += bump;
        } else {
            gp.cs[globalIdx - 2 * N] += bump;
        }
        EssviSurface surface(maturities, gp, cond);
        return surface.impliedVol(k, maturities[sliceIdx]);
    }

}

BOOST_AUTO_TEST_CASE(testEssviGlobalGradient) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI global-param analytic gradient vs finite differences...");

    std::vector<Real> maturities = {0.25, 0.5, 1.0};

    EssviGlobalParams gp;
    gp.rhos   = {-0.30, -0.25, -0.20};
    gp.theta1 = 0.01;
    gp.as     = {0.005, 0.010};
    gp.cs     = {0.50, 0.50, 0.50};

    auto cond = EssviButterflyCondition::GatheralJacquier;
    EssviSurface surface(maturities, gp, cond);

    Size N = gp.numSlices();
    Size nParams = 3 * N;  // 9 params total

    std::vector<Real> kValues = {-0.10, 0.0, 0.10};
    Real eps = 1e-6;
    Real tol = 1e-4;

    for (Size si = 0; si < N; ++si) {
        for (Real k : kValues) {
            std::vector<Real> grad = surface.impliedVolGlobalGradient(
                si, k, gp, cond);

            BOOST_REQUIRE_EQUAL(grad.size(), nParams);

            for (Size gi = 0; gi < nParams; ++gi) {
                // Check if bump keeps params in valid domain
                EssviGlobalParams gp_test = gp;
                bool valid = true;
                if (gi < N) {
                    if (gp.rhos[gi] + eps >= 1.0 || gp.rhos[gi] - eps <= -1.0)
                        valid = false;
                } else if (gi == N) {
                    if (gp.theta1 - eps <= 0.0)
                        valid = false;
                } else if (gi < 2 * N) {
                    if (gp.as[gi - N - 1] - eps <= 0.0)
                        valid = false;
                } else {
                    if (gp.cs[gi - 2 * N] + eps >= 1.0 || gp.cs[gi - 2 * N] - eps <= 0.0)
                        valid = false;
                }
                if (!valid) continue;

                Real vol_up = globalBumpedVol(
                    maturities, gp, gi, +eps, si, k, cond);
                Real vol_dn = globalBumpedVol(
                    maturities, gp, gi, -eps, si, k, cond);
                Real fd = (vol_up - vol_dn) / (2.0 * eps);

                BOOST_CHECK_SMALL(grad[gi] - fd, tol);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(testEssviTermStructureGradient) {

    BOOST_TEST_MESSAGE(
        "Testing EssviVolatilityTermStructure gradient vs finite differences...");

    Date refDate(17, March, 2026);
    Settings::instance().evaluationDate() = refDate;
    DayCounter dc = Actual365Fixed();

    Real spot = 100.0;
    Real rf = 0.04;
    Real div = 0.01;

    auto spotH = Handle<Quote>(ext::make_shared<SimpleQuote>(spot));
    auto rfH = Handle<YieldTermStructure>(
        ext::make_shared<FlatForward>(refDate, rf, dc));
    auto divH = Handle<YieldTermStructure>(
        ext::make_shared<FlatForward>(refDate, div, dc));

    std::vector<Date> dates = {
        refDate + Period(3, Months),
        refDate + Period(1, Years),
    };
    std::vector<Real> thetas = {0.01, 0.04};
    std::vector<Real> rhos   = {-0.30, -0.25};
    std::vector<Real> psis   = {0.05, 0.10};

    EssviVolatilityTermStructure surface(
        refDate, dates, thetas, rhos, psis, spotH, rfH, divH, dc);

    std::vector<Real> strikes = {85.0, 95.0, 100.0, 105.0, 115.0};
    Real eps = 1e-6;
    Real tol = 1e-5;

    for (Size si = 0; si < dates.size(); ++si) {
        for (Real K : strikes) {
            EssviSliceGradient grad = surface.impliedVolGradient(si, K);

            for (int param = 0; param < 3; ++param) {
                std::vector<Real> th_up = thetas, th_dn = thetas;
                std::vector<Real> rh_up = rhos,   rh_dn = rhos;
                std::vector<Real> ps_up = psis,   ps_dn = psis;

                if (param == 0) { th_up[si] += eps; th_dn[si] -= eps; }
                if (param == 1) { rh_up[si] += eps; rh_dn[si] -= eps; }
                if (param == 2) { ps_up[si] += eps; ps_dn[si] -= eps; }

                EssviVolatilityTermStructure s_up(
                    refDate, dates, th_up, rh_up, ps_up, spotH, rfH, divH, dc);
                EssviVolatilityTermStructure s_dn(
                    refDate, dates, th_dn, rh_dn, ps_dn, spotH, rfH, divH, dc);

                Real vol_up = s_up.blackVol(dates[si], K);
                Real vol_dn = s_dn.blackVol(dates[si], K);
                Real fd = (vol_up - vol_dn) / (2.0 * eps);

                Real analytic = (param == 0) ? grad.dSigma_dTheta
                              : (param == 1) ? grad.dSigma_dRho
                              :                grad.dSigma_dPsi;

                BOOST_CHECK_SMALL(analytic - fd, tol);
            }
        }
    }
}

// =====================================================================
// eSSVI analytic Dupire local vol vs QL numeric LocalVolSurface
// =====================================================================

BOOST_AUTO_TEST_CASE(testEssviLocalVolVsNumericDupire) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI analytic local vol vs QL numeric Dupire...");

    Date refDate(17, March, 2026);
    Settings::instance().evaluationDate() = refDate;
    DayCounter dc = Actual365Fixed();

    Real spot = 100.0;
    Real rf = 0.04;
    Real div = 0.01;

    auto spotQ = ext::make_shared<SimpleQuote>(spot);
    auto spotH = Handle<Quote>(spotQ);
    auto rfH = Handle<YieldTermStructure>(
        ext::make_shared<FlatForward>(refDate, rf, dc));
    auto divH = Handle<YieldTermStructure>(
        ext::make_shared<FlatForward>(refDate, div, dc));

    // 3-slice surface with realistic params
    std::vector<Date> dates = {
        refDate + Period(3, Months),
        refDate + Period(6, Months),
        refDate + Period(1, Years),
    };
    std::vector<Real> thetas = {0.01, 0.02, 0.04};
    std::vector<Real> rhos   = {-0.30, -0.25, -0.20};
    std::vector<Real> psis   = {0.05, 0.08, 0.12};

    auto essviSurface = ext::make_shared<EssviVolatilityTermStructure>(
        refDate, dates, thetas, rhos, psis, spotH, rfH, divH, dc);

    // Analytic local vol
    EssviLocalVolSurface analyticLV(essviSurface, rfH, divH, spotH);

    // Numeric local vol (QL's built-in Dupire)
    auto blackH = Handle<BlackVolTermStructure>(essviSurface);
    LocalVolSurface numericLV(blackH, rfH, divH, spotH);

    // Compare at a grid of (T, K) points
    // Use maturities between benchmark dates (interpolation region)
    std::vector<Time> times = {0.15, 0.25, 0.35, 0.50, 0.75, 1.0};
    std::vector<Real> strikes = {85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0};

    // QL's numeric Dupire uses finite differences with relative bump sizes.
    // Interior points (between benchmark maturities) agree very well.
    // At boundary/extrapolation regions, the FD time derivative differs
    // because QL adjusts strikes for changing forwards. Use 5% tolerance
    // to cover these edge effects.
    Real relTol = 0.05;

    for (Time t : times) {
        for (Real K : strikes) {
            Real analytic = analyticLV.localVol(t, K, true);
            Real numeric  = numericLV.localVol(t, K, true);

            Real diff = std::abs(analytic - numeric);
            Real scale = std::max(std::abs(analytic), 1e-6);

            BOOST_CHECK_MESSAGE(
                diff / scale < relTol,
                "local vol mismatch at T=" << t << ", K=" << K
                << ": analytic=" << analytic << ", numeric=" << numeric
                << ", rel diff=" << diff / scale);
        }
    }
}

BOOST_AUTO_TEST_CASE(testEssviLocalVolATM) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI local vol at ATM has expected structure...");

    Date refDate(17, March, 2026);
    Settings::instance().evaluationDate() = refDate;
    DayCounter dc = Actual365Fixed();

    Real spot = 100.0;
    auto spotH = Handle<Quote>(ext::make_shared<SimpleQuote>(spot));
    auto rfH = Handle<YieldTermStructure>(
        ext::make_shared<FlatForward>(refDate, 0.0, dc));
    auto divH = Handle<YieldTermStructure>(
        ext::make_shared<FlatForward>(refDate, 0.0, dc));

    // Flat smile (rho=0, uniform theta growth) → local vol ≈ implied vol
    std::vector<Date> dates = {
        refDate + Period(6, Months),
        refDate + Period(1, Years),
    };
    // Constant ATM vol: theta proportional to T → sigma = sqrt(theta/T) = const
    std::vector<Real> thetas = {0.02, 0.04};  // both give sigma = sqrt(0.04) = 0.2
    std::vector<Real> rhos   = {0.0, 0.0};
    std::vector<Real> psis   = {0.04, 0.06};

    auto essviSurface = ext::make_shared<EssviVolatilityTermStructure>(
        refDate, dates, thetas, rhos, psis, spotH, rfH, divH, dc);

    EssviLocalVolSurface lv(essviSurface, rfH, divH, spotH);

    // At ATM with zero rates, forward = spot, k = 0
    // For flat vol (theta linear in T), local vol at ATM ≈ implied vol
    Real impliedAtm = essviSurface->blackVol(0.75, spot);
    Real localAtm   = lv.localVol(0.75, spot, true);

    // Local vol should be close to implied vol for near-flat surface
    // (not exact because psi creates curvature)
    BOOST_CHECK_SMALL(localAtm - impliedAtm, 0.02);

    // Local vol should be positive everywhere
    BOOST_CHECK(localAtm > 0.0);
}

BOOST_AUTO_TEST_CASE(testEssviStrikeDerivatives) {

    BOOST_TEST_MESSAGE(
        "Testing eSSVI total variance strike derivatives vs finite differences...");

    std::vector<Real> maturities = {0.5, 1.0};
    std::vector<EssviSliceParams> slices = {
        {0.02, -0.30, 0.08},
        {0.04, -0.25, 0.12},
    };
    EssviSurface surface(maturities, slices);

    std::vector<Real> kValues = {-0.30, -0.10, 0.0, 0.10, 0.30};
    Real eps = 1e-6;
    Real tol = 1e-5;

    for (Real T : maturities) {
        for (Real k : kValues) {
            Real dwdk_a, d2wdk2_a;
            surface.totalVarianceStrikeDerivatives(k, T, dwdk_a, d2wdk2_a);

            // FD first derivative
            Real wp = surface.totalVariance(k + eps, T);
            Real wm = surface.totalVariance(k - eps, T);
            Real dwdk_fd = (wp - wm) / (2.0 * eps);

            // FD second derivative (use wider bump for better FD accuracy)
            Real w0 = surface.totalVariance(k, T);
            Real d2wdk2_fd = (wp - 2.0 * w0 + wm) / (eps * eps);

            BOOST_CHECK_SMALL(dwdk_a - dwdk_fd, tol);
            BOOST_CHECK_SMALL(d2wdk2_a - d2wdk2_fd, tol * 10.0);  // wider for 2nd deriv FD
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
