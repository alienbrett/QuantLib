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

#include "toplevelfixture.hpp"
#include "utilities.hpp"
#include <ql/termstructures/volatility/equityfx/parametricvoltermstructure.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <cmath>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(ParametricVolTermStructureTests)

namespace {

    // Constant shape: f(z) = 1, dfdz = 0.  Black vol == atmIv everywhere.
    class FlatShape : public ParametricVolShape {
      public:
        Real f(Real /*z*/, const std::vector<Real>& /*p*/) const override {
            return 1.0;
        }
        Real dfdz(Real /*z*/, const std::vector<Real>& /*p*/) const override {
            return 0.0;
        }
    };

    // Pure-quadratic shape in z: f(z) = 1 + a*z + b*z^2.
    // Used to verify the surface delegates correctly to a custom C++
    // subclass (analogue of a Python override).
    class QuadShape : public ParametricVolShape {
      public:
        Real f(Real z, const std::vector<Real>& p) const override {
            return 1.0 + p[0] * z + p[1] * z * z;
        }
        Real dfdz(Real z, const std::vector<Real>& p) const override {
            return p[0] + 2.0 * p[1] * z;
        }
    };

    struct Fixture {
        Date today;
        DayCounter dc;
        Handle<Quote> spot;
        Handle<YieldTermStructure> rTs;
        Handle<YieldTermStructure> qTs;

        Fixture()
        : today(Settings::instance().evaluationDate()),
          dc(Actual365Fixed()),
          spot(ext::make_shared<SimpleQuote>(100.0)),
          rTs(ext::make_shared<FlatForward>(today, 0.03, dc)),
          qTs(ext::make_shared<FlatForward>(today, 0.01, dc)) {}
    };

}  // namespace

// S3 shape arithmetic at z = 0 must give f = 1 and dfdz = s2.
BOOST_AUTO_TEST_CASE(testS3ShapeAtAtm) {
    BOOST_TEST_MESSAGE("Testing S3Shape values at z = 0...");
    S3Shape s3;
    std::vector<Real> p = {-0.3, 0.4};  // (s2, c2)
    QL_CHECK_CLOSE(s3.f(0.0, p), 1.0, 1e-12);
    QL_CHECK_CLOSE(s3.dfdz(0.0, p), p[0], 1e-12);
}

// Strike derivative agrees with central FD for several (z, params).
BOOST_AUTO_TEST_CASE(testS3ShapeDerivativeFiniteDifference) {
    BOOST_TEST_MESSAGE("Testing S3Shape dfdz against central FD...");
    S3Shape s3;
    std::vector<std::vector<Real>> paramSets = {
        {-0.3, 0.4},
        {-0.5, 0.8},
        { 0.1, 0.05},
        { 0.0, 0.2},
    };
    const Real h = 1e-5;
    for (const auto& p : paramSets) {
        for (Real z : {-2.0, -0.5, 0.0, 0.5, 2.0}) {
            Real fd = (s3.f(z + h, p) - s3.f(z - h, p)) / (2.0 * h);
            Real analytic = s3.dfdz(z, p);
            QL_CHECK_CLOSE(analytic, fd, 1e-6);
        }
    }
}

// Flat shape gives blackVol == atmIv across strikes and times.
BOOST_AUTO_TEST_CASE(testFlatShapeReproducesAtmIv) {
    BOOST_TEST_MESSAGE("Testing flat shape reproduces ATM iv at all (T,K)...");
    Fixture f;
    std::vector<Date> dates = {f.today + 60, f.today + 180, f.today + 365};
    std::vector<ParametricVolSlice> slices = {
        {0.22, {}},
        {0.20, {}},
        {0.18, {}},
    };
    auto shape = ext::make_shared<FlatShape>();
    ParametricVolTermStructure surface(
        f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);

    // ATM == per-slice atmIv at the pillar.
    for (Size i = 0; i < dates.size(); ++i) {
        Real F = surface.forward(f.dc.yearFraction(f.today, dates[i]));
        Real v = surface.blackVol(dates[i], F);
        QL_CHECK_CLOSE(v, slices[i].atmIv, 1e-10);
    }

    // Off-ATM also flat (FlatShape says f == 1 anywhere).
    for (Real K : {80.0, 95.0, 105.0, 120.0}) {
        Real v = surface.blackVol(dates[1], K);
        QL_CHECK_CLOSE(v, slices[1].atmIv, 1e-10);
    }
}

// S3 shape: at the pillar, ATM vol equals atmIv (f(0) = 1).
BOOST_AUTO_TEST_CASE(testS3ShapeAtmIdentity) {
    BOOST_TEST_MESSAGE("Testing S3 shape preserves ATM iv at pillar...");
    Fixture f;
    std::vector<Date> dates = {f.today + 30, f.today + 90, f.today + 365};
    std::vector<ParametricVolSlice> slices = {
        {0.30, {-0.5, 0.6}},
        {0.25, {-0.4, 0.4}},
        {0.20, {-0.3, 0.3}},
    };
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);

    for (Size i = 0; i < dates.size(); ++i) {
        Time t = f.dc.yearFraction(f.today, dates[i]);
        Real F = surface.forward(t);
        Real v = surface.blackVol(dates[i], F);
        QL_CHECK_CLOSE(v, slices[i].atmIv, 1e-10);
    }
}

// At a pillar, blackVol matches sqrt(atmIv^2 * f(z)) when computed by hand.
BOOST_AUTO_TEST_CASE(testS3OffAtmAgainstClosedForm) {
    BOOST_TEST_MESSAGE("Testing S3 off-ATM vol matches closed form...");
    Fixture f;
    Date pillar = f.today + 365;
    std::vector<Date> dates = {pillar};
    ParametricVolSlice slc{0.25, {-0.4, 0.3}};
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, dates, {slc}, shape, f.spot, f.rTs, f.qTs, f.dc);

    Time T = f.dc.yearFraction(f.today, pillar);
    Real F = surface.forward(T);
    Real sigmaHat = slc.atmIv * std::sqrt(T);
    for (Real K : {70.0, 90.0, 110.0, 130.0}) {
        Real k = std::log(K / F);
        Real z = k / sigmaHat;
        Real expected = slc.atmIv * std::sqrt(shape->f(z, slc.params));
        Real got = surface.blackVol(pillar, K);
        QL_CHECK_CLOSE(got, expected, 1e-10);
    }
}

// Linear-in-T total-variance interpolation between two pillars.
BOOST_AUTO_TEST_CASE(testTotalVarianceTimeInterpolation) {
    BOOST_TEST_MESSAGE("Testing total variance interpolates linearly in T...");
    Fixture f;
    Date d1 = f.today + 90;
    Date d2 = f.today + 365;
    std::vector<Date> dates = {d1, d2};
    std::vector<ParametricVolSlice> slices = {
        {0.30, {-0.4, 0.3}},
        {0.20, {-0.3, 0.2}},
    };
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);

    Time T1 = f.dc.yearFraction(f.today, d1);
    Time T2 = f.dc.yearFraction(f.today, d2);
    Time tMid = 0.5 * (T1 + T2);
    Real k = -0.05;  // arbitrary log-moneyness

    Real w1 = surface.totalVariance(k, T1);
    Real w2 = surface.totalVariance(k, T2);
    Real wMid = surface.totalVariance(k, tMid);
    Real expected = 0.5 * (w1 + w2);
    QL_CHECK_CLOSE(wMid, expected, 1e-12);
}

// Polymorphism: surface delegates to a custom C++ subclass.
BOOST_AUTO_TEST_CASE(testCustomShapeSubclass) {
    BOOST_TEST_MESSAGE("Testing surface delegates to custom shape subclass...");
    Fixture f;
    Date pillar = f.today + 180;
    std::vector<Date> dates = {pillar};
    ParametricVolSlice slc{0.20, {0.10, 0.30}};  // (a, b) for QuadShape
    auto shape = ext::make_shared<QuadShape>();
    ParametricVolTermStructure surface(
        f.today, dates, {slc}, shape, f.spot, f.rTs, f.qTs, f.dc);

    Time T = f.dc.yearFraction(f.today, pillar);
    Real F = surface.forward(T);
    Real sigmaHat = slc.atmIv * std::sqrt(T);
    for (Real K : {85.0, 100.0, 115.0}) {
        Real k = std::log(K / F);
        Real z = k / sigmaHat;
        Real fz = 1.0 + slc.params[0] * z + slc.params[1] * z * z;
        Real expected = slc.atmIv * std::sqrt(fz);
        Real got = surface.blackVol(pillar, K);
        QL_CHECK_CLOSE(got, expected, 1e-10);
    }
}

// JW reduces to S3 when c_minus = c_plus.
BOOST_AUTO_TEST_CASE(testJWShapeReducesToS3) {
    BOOST_TEST_MESSAGE("Testing JW with c_minus = c_plus matches S3...");
    S3Shape s3;
    JWShape jw;
    std::vector<Real> s3p = {-0.4, 0.3};
    std::vector<Real> jwp = {-0.4, 0.3, 0.3};
    for (Real z : {-2.0, -0.5, 0.0, 0.5, 2.0}) {
        QL_CHECK_CLOSE(jw.f(z, jwp), s3.f(z, s3p), 1e-12);
        QL_CHECK_CLOSE(jw.dfdz(z, jwp), s3.dfdz(z, s3p), 1e-12);
    }
}

// JW with c_plus != c_minus uses the right kernel on each side.
BOOST_AUTO_TEST_CASE(testJWShapeAsymmetricWings) {
    BOOST_TEST_MESSAGE("Testing JW with asymmetric wing curvatures...");
    JWShape jw;
    std::vector<Real> p = {-0.3, 0.5, 0.2};  // s2, c_minus, c_plus
    // At z = 0: f = 1, df = s2 (independent of c+/c-).
    QL_CHECK_CLOSE(jw.f(0.0, p), 1.0, 1e-12);
    QL_CHECK_CLOSE(jw.dfdz(0.0, p), p[0], 1e-12);
    // At z = +1: should match S3 with c2 = c_plus.
    S3Shape s3;
    std::vector<Real> s3plus = {p[0], p[2]};
    QL_CHECK_CLOSE(jw.f(1.0, p), s3.f(1.0, s3plus), 1e-12);
    // At z = -1: should match S3 with c2 = c_minus.
    std::vector<Real> s3minus = {p[0], p[1]};
    QL_CHECK_CLOSE(jw.f(-1.0, p), s3.f(-1.0, s3minus), 1e-12);
}

// Butterfly density: BS limit is g(z) = 1 everywhere.
BOOST_AUTO_TEST_CASE(testButterflyDensityFlatShape) {
    BOOST_TEST_MESSAGE("Testing butterfly density g(z) for flat smile = 1...");
    Fixture f;
    Date pillar = f.today + 365;
    ParametricVolSlice slc{0.20, {}};
    auto shape = ext::make_shared<FlatShape>();
    ParametricVolTermStructure surface(
        f.today, {pillar}, {slc}, shape, f.spot, f.rTs, f.qTs, f.dc);

    for (Real z : {-2.0, -0.5, 0.0, 0.5, 2.0}) {
        Real g = surface.butterflyDensity(0, z);
        QL_CHECK_CLOSE(g, 1.0, 1e-6);
    }
    std::vector<Real> zGrid = {-3.0, -1.0, 0.0, 1.0, 3.0};
    QL_CHECK_CLOSE(surface.butterflyArbViolation(0, zGrid), 0.0, 1e-12);
}

// S3 with reasonable params: density should be > 0 across wide z range.
BOOST_AUTO_TEST_CASE(testS3ArbFreeOnTypicalParams) {
    BOOST_TEST_MESSAGE("Testing S3 with typical equity smile is arb-free...");
    Fixture f;
    Date pillar = f.today + 180;
    ParametricVolSlice slc{0.20, {-0.4, 0.3}};
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, {pillar}, {slc}, shape, f.spot, f.rTs, f.qTs, f.dc);

    std::vector<Real> zGrid;
    for (Real z = -4.0; z <= 4.01; z += 0.25) zGrid.push_back(z);
    Real viol = surface.butterflyArbViolation(0, zGrid);
    BOOST_CHECK_MESSAGE(viol == 0.0,
                       "S3 with typical params should be arb-free, got "
                       << viol);
}

// Detect a butterfly violation: pushing s2 too negative breaks ATF bound.
BOOST_AUTO_TEST_CASE(testButterflyViolationDetected) {
    BOOST_TEST_MESSAGE("Testing butterfly violation detected for extreme s2...");
    Fixture f;
    Date pillar = f.today + 365;
    // ATF bound: s2^2 <= (4 + 2*c2) / (1 + sigma_hat0^2 / 4).
    // With c2 = 0, sigma_hat = atm_iv * sqrt(T) ~ 0.30 * 1 = 0.30:
    //   bound = 4 / (1 + 0.0225) = ~3.91 => |s2| ≲ 1.98.
    // Pick s2 = -2.5 to push past it.
    ParametricVolSlice slc{0.30, {-2.5, 0.0}};
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, {pillar}, {slc}, shape, f.spot, f.rTs, f.qTs, f.dc);
    Real g0 = surface.butterflyDensity(0, 0.0);
    BOOST_CHECK_MESSAGE(g0 < 0.0,
                       "expected ATF butterfly violation (g(0) < 0), got "
                       << g0);
}

// Calendar arb check: monotone w(T) => 0; flipped slices => positive.
BOOST_AUTO_TEST_CASE(testCalendarArbCheck) {
    BOOST_TEST_MESSAGE("Testing calendar arb detection between two pillars...");
    Fixture f;
    std::vector<Date> dates = {f.today + 90, f.today + 365};
    std::vector<Real> kGrid = {-0.5, -0.25, 0.0, 0.25, 0.5};

    // Monotone: T=90d atmIv=0.18, T=365d atmIv=0.22 => w_lo < w_hi.
    {
        std::vector<ParametricVolSlice> slices = {
            {0.18, {-0.2, 0.3}}, {0.22, {-0.2, 0.3}},
        };
        auto shape = ext::make_shared<S3Shape>();
        ParametricVolTermStructure surface(
            f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);
        QL_CHECK_CLOSE(surface.calendarArbViolation(kGrid), 0.0, 1e-12);
    }

    // Reverse: T=90d atmIv=0.40 (high), T=365d atmIv=0.10 (low) => violation.
    {
        std::vector<ParametricVolSlice> slices = {
            {0.40, {-0.2, 0.3}}, {0.10, {-0.2, 0.3}},
        };
        auto shape = ext::make_shared<S3Shape>();
        ParametricVolTermStructure surface(
            f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);
        Real v = surface.calendarArbViolation(kGrid);
        BOOST_CHECK_MESSAGE(v > 0.0,
                           "expected calendar violation > 0, got " << v);
    }
}

// Analytic dfdParams must match central FD on f for S3.
BOOST_AUTO_TEST_CASE(testS3DfdParamsMatchesFD) {
    BOOST_TEST_MESSAGE("Testing S3 dfdParams vs central FD...");
    S3Shape s3;
    std::vector<Real> p = {-0.4, 0.3};
    const Real h = 1e-6;
    for (Real z : {-2.0, -0.5, 0.0, 0.5, 2.0}) {
        auto analytic = s3.dfdParams(z, p);
        for (Size i = 0; i < p.size(); ++i) {
            auto pp = p; pp[i] += h;
            Real fp = s3.f(z, pp);
            pp[i] = p[i] - h;
            Real fm = s3.f(z, pp);
            Real fd = (fp - fm) / (2.0 * h);
            QL_CHECK_CLOSE(analytic[i], fd, 1e-6);
        }
    }
}

// Analytic dfdParams must match central FD on f for JW (away from z=0).
BOOST_AUTO_TEST_CASE(testJWDfdParamsMatchesFD) {
    BOOST_TEST_MESSAGE("Testing JW dfdParams vs central FD...");
    JWShape jw;
    std::vector<Real> p = {-0.3, 0.5, 0.2};
    const Real h = 1e-6;
    for (Real z : {-2.0, -1.0, -0.5, 0.5, 1.0, 2.0}) {
        auto analytic = jw.dfdParams(z, p);
        for (Size i = 0; i < p.size(); ++i) {
            auto pp = p; pp[i] += h;
            Real fp = jw.f(z, pp);
            pp[i] = p[i] - h;
            Real fm = jw.f(z, pp);
            Real fd = (fp - fm) / (2.0 * h);
            QL_CHECK_CLOSE(analytic[i], fd, 1e-5);
        }
    }
}

// Batch: batchBlackVol matches scalar blackVol at each (sliceIdx, K).
BOOST_AUTO_TEST_CASE(testBatchBlackVolMatchesScalar) {
    BOOST_TEST_MESSAGE("Testing batchBlackVol matches scalar blackVol...");
    Fixture f;
    std::vector<Date> dates = {f.today + 60, f.today + 180, f.today + 365};
    std::vector<ParametricVolSlice> slices = {
        {0.25, {-0.4, 0.30}}, {0.22, {-0.3, 0.25}}, {0.20, {-0.2, 0.20}},
    };
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);

    std::vector<Size> sl;
    std::vector<Real> Ks;
    for (Size i = 0; i < dates.size(); ++i) {
        Real F = surface.forward(f.dc.yearFraction(f.today, dates[i]));
        for (Real km : {-0.2, -0.05, 0.0, 0.05, 0.2}) {
            sl.push_back(i);
            Ks.push_back(F * std::exp(km));
        }
    }
    std::vector<Real> batch = surface.batchBlackVol(sl, Ks);
    for (Size j = 0; j < sl.size(); ++j) {
        Real scalar = surface.blackVol(dates[sl[j]], Ks[j]);
        QL_CHECK_CLOSE(batch[j], scalar, 1e-10);
    }
}

// Batch arb sanity: batched eval doesn't perturb arb checks.
BOOST_AUTO_TEST_CASE(testBatchEvalArbInvariant) {
    BOOST_TEST_MESSAGE("Testing arb checks invariant under batch eval...");
    Fixture f;
    std::vector<Date> dates = {f.today + 60, f.today + 180, f.today + 365};
    std::vector<ParametricVolSlice> slices = {
        {0.25, {-0.4, 0.30}}, {0.22, {-0.3, 0.25}}, {0.20, {-0.2, 0.20}},
    };
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, dates, slices, shape, f.spot, f.rTs, f.qTs, f.dc);

    std::vector<Real> zGrid;
    for (Real z = -4.0; z <= 4.01; z += 0.25) zGrid.push_back(z);
    std::vector<Real> kGrid = {-0.5, -0.25, 0.0, 0.25, 0.5};

    Real bf0 = surface.butterflyArbViolation(1, zGrid);
    Real cal0 = surface.calendarArbViolation(kGrid);

    // Exercise the batch path.
    std::vector<Size> sl(zGrid.size(), 1);
    std::vector<Real> Ks(zGrid.size());
    Real F = surface.forward(f.dc.yearFraction(f.today, dates[1]));
    for (Size j = 0; j < zGrid.size(); ++j)
        Ks[j] = F * std::exp(zGrid[j] * slices[1].atmIv
                              * std::sqrt(f.dc.yearFraction(f.today, dates[1])));
    surface.batchBlackVol(sl, Ks);
    surface.batchBlackVolAtTimes(
        std::vector<Real>(Ks.size(),
                          f.dc.yearFraction(f.today, dates[1])),
        Ks);

    Real bf1 = surface.butterflyArbViolation(1, zGrid);
    Real cal1 = surface.calendarArbViolation(kGrid);
    QL_CHECK_CLOSE(bf0, bf1, 1e-12);
    QL_CHECK_CLOSE(cal0, cal1, 1e-12);
}

// Mutator path: setSlice updates observers and the next query reflects it.
BOOST_AUTO_TEST_CASE(testSetSliceUpdates) {
    BOOST_TEST_MESSAGE("Testing setSlice mutator updates blackVol...");
    Fixture f;
    Date pillar = f.today + 90;
    std::vector<Date> dates = {pillar};
    ParametricVolSlice slc{0.20, {-0.3, 0.2}};
    auto shape = ext::make_shared<S3Shape>();
    ParametricVolTermStructure surface(
        f.today, dates, {slc}, shape, f.spot, f.rTs, f.qTs, f.dc);

    Time T = f.dc.yearFraction(f.today, pillar);
    Real F = surface.forward(T);
    Real v0 = surface.blackVol(pillar, F);
    QL_CHECK_CLOSE(v0, 0.20, 1e-10);

    surface.setSlice(0, 0.35, std::vector<Real>{-0.5, 0.4});
    Real v1 = surface.blackVol(pillar, F);
    QL_CHECK_CLOSE(v1, 0.35, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
