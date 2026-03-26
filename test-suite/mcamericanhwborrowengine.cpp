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

#include "toplevelfixture.hpp"
#include "utilities.hpp"
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/mcamericanhwborrowengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/processes/blackscholeshwborrowprocess.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/calendars/nullcalendar.hpp>
#include <cmath>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(MCAmericanHWBorrowEngineTests)

BOOST_AUTO_TEST_CASE(testZeroVolLimitConvergence) {
    BOOST_TEST_MESSAGE("Testing MC American HW Borrow engine convergence"
                       " to FD Black-Scholes in zero rate-vol limit...");

    DayCounter dc = Actual365Fixed();
    Date today(15, January, 2025);
    Settings::instance().evaluationDate() = today;
    Date maturity = today + Period(1, Years);

    // Market data
    Real spotVal = 100.0;
    Real strikeVal = 100.0;
    Volatility eqVol = 0.25;
    Rate rateVal = 0.04;
    Rate borrowVal = 0.005;

    // Rate model with near-zero vol (deterministic limit)
    Real meanReversion = 0.03;
    Real rateVol = 1e-6;
    Real borrowVol = 1e-6;

    auto spot = ext::make_shared<SimpleQuote>(spotVal);
    Handle<Quote> spotHandle(spot);

    Handle<BlackVolTermStructure> volSurface(
        ext::make_shared<BlackConstantVol>(today, NullCalendar(), eqVol, dc));

    Handle<YieldTermStructure> rateTS(
        ext::make_shared<FlatForward>(today, rateVal, dc));
    Handle<YieldTermStructure> borrowTS(
        ext::make_shared<FlatForward>(today, borrowVal, dc));

    auto rateModel = ext::make_shared<HullWhite>(rateTS, meanReversion, rateVol);
    auto borrowModel = ext::make_shared<HullWhite>(borrowTS, meanReversion, borrowVol);

    auto corrSR = ext::make_shared<SimpleQuote>(0.0);
    auto corrSB = ext::make_shared<SimpleQuote>(0.0);
    auto corrRB = ext::make_shared<SimpleQuote>(0.0);

    auto process = ext::make_shared<BlackScholesHWBorrowProcess>(
        spotHandle, volSurface,
        rateModel, borrowModel,
        Handle<Quote>(corrSR), Handle<Quote>(corrSB), Handle<Quote>(corrRB));

    // MC engine: 100k paths, 50 time steps
    auto mcEngine = ext::make_shared<MCPRAmericanHWBorrowEngine>(
        process, 50, 100000, 42);

    // FD benchmark: standard Black-Scholes with dividend yield = borrow rate
    auto bsProcess = ext::make_shared<BlackScholesMertonProcess>(
        spotHandle, borrowTS, rateTS, volSurface);
    auto fdEngine = ext::make_shared<FdBlackScholesVanillaEngine>(
        bsProcess, 200, 200, 0);

    // American put
    VanillaOption mcOption(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
        ext::make_shared<AmericanExercise>(today, maturity));
    mcOption.setPricingEngine(mcEngine);
    Real mcPrice = mcOption.NPV();

    VanillaOption fdOption(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
        ext::make_shared<AmericanExercise>(today, maturity));
    fdOption.setPricingEngine(fdEngine);
    Real fdPrice = fdOption.NPV();

    Real relError = std::abs(mcPrice - fdPrice) / fdPrice;

    BOOST_TEST_MESSAGE("  MC price:  " << mcPrice);
    BOOST_TEST_MESSAGE("  FD price:  " << fdPrice);
    BOOST_TEST_MESSAGE("  Rel error: " << relError);

    // 0.5% relative tolerance
    BOOST_CHECK_SMALL(relError, 0.005);
}


BOOST_AUTO_TEST_CASE(testStochasticRatesEffect) {
    BOOST_TEST_MESSAGE("Testing that stochastic rate vol changes the price...");

    DayCounter dc = Actual365Fixed();
    Date today(15, January, 2025);
    Settings::instance().evaluationDate() = today;
    Date maturity = today + Period(1, Years);

    Real spotVal = 100.0;
    Real strikeVal = 100.0;
    Volatility eqVol = 0.25;
    Rate rateVal = 0.04;
    Rate borrowVal = 0.005;
    Real meanReversion = 0.03;

    auto spot = ext::make_shared<SimpleQuote>(spotVal);
    Handle<Quote> spotHandle(spot);
    Handle<BlackVolTermStructure> volSurface(
        ext::make_shared<BlackConstantVol>(today, NullCalendar(), eqVol, dc));
    Handle<YieldTermStructure> rateTS(
        ext::make_shared<FlatForward>(today, rateVal, dc));
    Handle<YieldTermStructure> borrowTS(
        ext::make_shared<FlatForward>(today, borrowVal, dc));

    auto corrSR = ext::make_shared<SimpleQuote>(-0.20);
    auto corrSB = ext::make_shared<SimpleQuote>(0.10);
    auto corrRB = ext::make_shared<SimpleQuote>(0.15);

    // Near-zero vol baseline
    auto rateModelLow = ext::make_shared<HullWhite>(rateTS, meanReversion, 1e-6);
    auto borrowModelLow = ext::make_shared<HullWhite>(borrowTS, meanReversion, 1e-6);

    auto processLow = ext::make_shared<BlackScholesHWBorrowProcess>(
        spotHandle, volSurface,
        rateModelLow, borrowModelLow,
        Handle<Quote>(corrSR), Handle<Quote>(corrSB), Handle<Quote>(corrRB));

    auto engineLow = ext::make_shared<MCPRAmericanHWBorrowEngine>(
        processLow, 50, 50000, 42);

    // Meaningful rate/borrow vol with correlations
    auto rateModelHigh = ext::make_shared<HullWhite>(rateTS, meanReversion, 0.03);
    auto borrowModelHigh = ext::make_shared<HullWhite>(borrowTS, meanReversion, 0.05);

    auto processHigh = ext::make_shared<BlackScholesHWBorrowProcess>(
        spotHandle, volSurface,
        rateModelHigh, borrowModelHigh,
        Handle<Quote>(corrSR), Handle<Quote>(corrSB), Handle<Quote>(corrRB));

    auto engineHigh = ext::make_shared<MCPRAmericanHWBorrowEngine>(
        processHigh, 50, 50000, 42);

    VanillaOption optionLow(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
        ext::make_shared<AmericanExercise>(today, maturity));
    optionLow.setPricingEngine(engineLow);
    Real priceLow = optionLow.NPV();

    VanillaOption optionHigh(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
        ext::make_shared<AmericanExercise>(today, maturity));
    optionHigh.setPricingEngine(engineHigh);
    Real priceHigh = optionHigh.NPV();

    BOOST_TEST_MESSAGE("  Low rate vol price:  " << priceLow);
    BOOST_TEST_MESSAGE("  High rate vol price: " << priceHigh);

    // Prices should differ by a detectable amount
    BOOST_CHECK(std::abs(priceHigh - priceLow) > 0.01);
}


BOOST_AUTO_TEST_CASE(testDeepITMIntrinsicValue) {
    BOOST_TEST_MESSAGE("Testing deep ITM American put >= intrinsic value...");

    DayCounter dc = Actual365Fixed();
    Date today(15, January, 2025);
    Settings::instance().evaluationDate() = today;
    Date maturity = today + Period(1, Years);

    Real spotVal = 100.0;
    Real strikeVal = 130.0;
    Volatility eqVol = 0.25;
    Rate rateVal = 0.04;
    Rate borrowVal = 0.005;

    auto spot = ext::make_shared<SimpleQuote>(spotVal);
    Handle<Quote> spotHandle(spot);
    Handle<BlackVolTermStructure> volSurface(
        ext::make_shared<BlackConstantVol>(today, NullCalendar(), eqVol, dc));
    Handle<YieldTermStructure> rateTS(
        ext::make_shared<FlatForward>(today, rateVal, dc));
    Handle<YieldTermStructure> borrowTS(
        ext::make_shared<FlatForward>(today, borrowVal, dc));

    auto rateModel = ext::make_shared<HullWhite>(rateTS, 0.03, 0.01);
    auto borrowModel = ext::make_shared<HullWhite>(borrowTS, 0.50, 0.02);

    auto corrSR = ext::make_shared<SimpleQuote>(-0.15);
    auto corrSB = ext::make_shared<SimpleQuote>(0.05);
    auto corrRB = ext::make_shared<SimpleQuote>(0.10);

    auto process = ext::make_shared<BlackScholesHWBorrowProcess>(
        spotHandle, volSurface,
        rateModel, borrowModel,
        Handle<Quote>(corrSR), Handle<Quote>(corrSB), Handle<Quote>(corrRB));

    auto engine = ext::make_shared<MCPRAmericanHWBorrowEngine>(
        process, 50, 50000, 42);

    VanillaOption option(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
        ext::make_shared<AmericanExercise>(today, maturity));
    option.setPricingEngine(engine);
    Real price = option.NPV();

    Real intrinsic = strikeVal - spotVal;

    BOOST_TEST_MESSAGE("  Price:     " << price);
    BOOST_TEST_MESSAGE("  Intrinsic: " << intrinsic);

    BOOST_CHECK(price >= intrinsic * 0.99);
}


BOOST_AUTO_TEST_CASE(testObservableNotification) {
    BOOST_TEST_MESSAGE("Testing observable notification on correlation change...");

    DayCounter dc = Actual365Fixed();
    Date today(15, January, 2025);
    Settings::instance().evaluationDate() = today;
    Date maturity = today + Period(1, Years);

    auto spot = ext::make_shared<SimpleQuote>(100.0);
    Handle<Quote> spotHandle(spot);
    Handle<BlackVolTermStructure> volSurface(
        ext::make_shared<BlackConstantVol>(today, NullCalendar(), 0.25, dc));
    Handle<YieldTermStructure> rateTS(
        ext::make_shared<FlatForward>(today, 0.04, dc));
    Handle<YieldTermStructure> borrowTS(
        ext::make_shared<FlatForward>(today, 0.005, dc));

    auto rateModel = ext::make_shared<HullWhite>(rateTS, 0.03, 0.01);
    auto borrowModel = ext::make_shared<HullWhite>(borrowTS, 0.50, 0.02);

    auto corrSR = ext::make_shared<SimpleQuote>(-0.15);
    auto corrSB = ext::make_shared<SimpleQuote>(0.05);
    auto corrRB = ext::make_shared<SimpleQuote>(0.10);

    auto process = ext::make_shared<BlackScholesHWBorrowProcess>(
        spotHandle, volSurface,
        rateModel, borrowModel,
        Handle<Quote>(corrSR), Handle<Quote>(corrSB), Handle<Quote>(corrRB));

    // Use fewer paths for speed (we're testing notification, not convergence)
    auto engine = ext::make_shared<MCPRAmericanHWBorrowEngine>(
        process, 30, 10000, 42);

    VanillaOption option(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, 100.0),
        ext::make_shared<AmericanExercise>(today, maturity));
    option.setPricingEngine(engine);

    Real price1 = option.NPV();

    // Change correlation
    corrSR->setValue(-0.30);
    Real price2 = option.NPV();

    BOOST_TEST_MESSAGE("  Price (rho=-0.15): " << price1);
    BOOST_TEST_MESSAGE("  Price (rho=-0.30): " << price2);

    // The price should have changed
    BOOST_CHECK(price1 != price2);
}

BOOST_AUTO_TEST_CASE(testConvexityBiasMonotonicity) {
    BOOST_TEST_MESSAGE("Testing hidden optionality (convexity bias) is "
                       "monotone in rate vol [El Hassan, Maddah, Taleb 2026]...");

    // Replicate Table I setup from arXiv:2602.14350v1:
    // ATM equity put, S=100, K=100, T=1Y, sigma=40%, r~4.18%, no dividend.
    // The paper finds that the American put price increases monotonically
    // with the standard deviation of the stochastic interest rate.
    // At sigma_r=0: O_A = 14.3184, at sigma_r=4.78%: O_A = 14.5319 (pi_A=0.2136)

    DayCounter dc = Actual365Fixed();
    Date today(15, January, 2025);
    Settings::instance().evaluationDate() = today;
    Date maturity = today + Period(1, Years);

    Real spotVal = 100.0;
    Real strikeVal = 100.0;
    Volatility eqVol = 0.40;
    Rate rateVal = 0.0418;   // paper's expected rate at the fugit
    Rate borrowVal = 1e-6;   // equity put, no dividends

    Real meanReversion = 0.03;

    auto spot = ext::make_shared<SimpleQuote>(spotVal);
    Handle<Quote> spotHandle(spot);
    Handle<BlackVolTermStructure> volSurface(
        ext::make_shared<BlackConstantVol>(today, NullCalendar(), eqVol, dc));
    Handle<YieldTermStructure> rateTS(
        ext::make_shared<FlatForward>(today, rateVal, dc));
    Handle<YieldTermStructure> borrowTS(
        ext::make_shared<FlatForward>(today, borrowVal, dc));

    // Test monotonicity: increasing rate vol should increase the price
    // (convexity in the rate dimension favors American exercise)
    Real rateVols[] = { 1e-6, 0.01, 0.02, 0.03, 0.05 };
    Real prices[5];

    for (Size k = 0; k < 5; ++k) {
        auto rateModel = ext::make_shared<HullWhite>(
            rateTS, meanReversion, rateVols[k]);
        auto borrowModel = ext::make_shared<HullWhite>(
            borrowTS, meanReversion, 1e-6);

        auto process = ext::make_shared<BlackScholesHWBorrowProcess>(
            spotHandle, volSurface,
            rateModel, borrowModel,
            Handle<Quote>(ext::make_shared<SimpleQuote>(0.0)),
            Handle<Quote>(ext::make_shared<SimpleQuote>(0.0)),
            Handle<Quote>(ext::make_shared<SimpleQuote>(0.0)));

        auto engine = ext::make_shared<MCPRAmericanHWBorrowEngine>(
            process, 50, 80000, 42);

        VanillaOption option(
            ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
            ext::make_shared<AmericanExercise>(today, maturity));
        option.setPricingEngine(engine);
        prices[k] = option.NPV();

        BOOST_TEST_MESSAGE("  rate_vol=" << rateVols[k]
                           << "  price=" << prices[k]);
    }

    // Note: the paper uses r0=1% drifting to r_bar=4.18%, hence their
    // O_A(r*) = 14.3184 at an effective rate of ~3.35%.  Our flat curve
    // at 4.18% gives a lower put price (~13.98), as expected.  The key
    // qualitative finding — monotone convexity bias — is what we validate.

    // MC deterministic baseline should match FD
    auto bsProcess = ext::make_shared<BlackScholesMertonProcess>(
        spotHandle, borrowTS, rateTS, volSurface);
    auto fdEngine = ext::make_shared<FdBlackScholesVanillaEngine>(
        bsProcess, 200, 200, 0);
    VanillaOption fdOption(
        ext::make_shared<PlainVanillaPayoff>(Option::Put, strikeVal),
        ext::make_shared<AmericanExercise>(today, maturity));
    fdOption.setPricingEngine(fdEngine);
    Real fdPrice = fdOption.NPV();
    BOOST_TEST_MESSAGE("  FD benchmark:  " << fdPrice);
    BOOST_CHECK_SMALL(std::abs(prices[0] - fdPrice) / fdPrice, 0.005);

    // Monotonicity: each price should be >= previous (within MC noise)
    for (Size k = 1; k < 5; ++k) {
        BOOST_CHECK(prices[k] >= prices[k-1] - 0.05);
    }

    // The convexity bias at the highest rate vol should be positive
    Real piA = prices[4] - prices[0];
    BOOST_TEST_MESSAGE("  Convexity bias (sigma_r=5%): " << piA);
    BOOST_CHECK(piA > 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
