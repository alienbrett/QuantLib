/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 chloride contributors

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

/*! \file zabrswaptionvol.hpp
    \brief Swaption vol surface backed by ZABR smile sections with
           bilinear parameter interpolation across (expiry, tail) grid.

    Native ZABR evaluation at every strike — no tabulated interpolation.
    At gamma=1.0, this reduces to standard SABR (Hagan formula is used
    for alpha calibration; ZabrLocalVolatility for the smile section).

    Intended for use with LinearTsrPricer for CMS convexity-adjusted
    pricing where the full tail distribution matters.
*/

#ifndef quantlib_zabr_swaption_vol_hpp
#define quantlib_zabr_swaption_vol_hpp

#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/experimental/volatility/zabrsmilesection.hpp>
#include <ql/termstructures/volatility/sabr.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/time/schedule.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/instruments/makeois.hpp>
#include <ql/indexes/swapindex.hpp>
#include <utility>

namespace QuantLib {

    //! ZABR swaption volatility surface with bilinear parameter interpolation.
    /*!
        Stores SABR-compatible parameters (ATM normal vol, beta, nu, rho) on a
        2D grid of (option expiry, swap tail) and interpolates bilinearly.
        At each query point, calibrates alpha to match the interpolated ATM
        normal vol and constructs a ZabrSmileSection<ZabrLocalVolatility>.

        The gamma parameter controls tail weight:
        - gamma = 1.0: equivalent to standard SABR
        - gamma > 1.0: fatter tails (important for CMS convexity)

        Usage:
        \code
        auto vol = ext::make_shared<ZabrSwaptionVol>(
            referenceDate, calendar,
            expiryYears, tailYears,
            atmNormalVolMatrix, betaMatrix, nuMatrix, rhoMatrix,
            gamma, yieldCurve, swapIndex);
        Handle<SwaptionVolatilityStructure> h(vol);
        auto pricer = ext::make_shared<LinearTsrPricer>(h, meanReversion);
        \endcode
    */
    class ZabrSwaptionVol : public SwaptionVolatilityStructure {
      public:
        ZabrSwaptionVol(const Date& referenceDate,
                        const Calendar& calendar,
                        const std::vector<Real>& expiryYears,
                        const std::vector<Real>& tailYears,
                        const Matrix& atmNormalVol,
                        const Matrix& beta,
                        const Matrix& nu,
                        const Matrix& rho,
                        Real gamma,
                        Handle<YieldTermStructure> yieldCurve,
                        ext::shared_ptr<SwapIndex> swapIndex);

        //! \name TermStructure interface
        //@{
        Date maxDate() const override { return Date::maxDate(); }
        //@}

        //! \name VolatilityTermStructure interface
        //@{
        Rate minStrike() const override { return 0.0001; }
        Rate maxStrike() const override { return 10.0; }
        //@}

        //! \name SwaptionVolatilityStructure interface
        //@{
        const Period& maxSwapTenor() const override { return maxTenor_; }
        VolatilityType volatilityType() const override { return ShiftedLognormal; }
        //@}

      protected:
        Real shiftImpl(Time, Time) const override { return 0.0; }

        Volatility volatilityImpl(Time optionTime, Time swapLength,
                                  Rate strike) const override;

        ext::shared_ptr<SmileSection>
        smileSectionImpl(Time optionTime, Time swapLength) const override;

      private:
        Real forwardSwapRate(Time optionTime, Time swapLength) const;
        Real calibrateAlpha(Real forward, Real atmNormalVol, Time T,
                            Real beta, Real nu, Real rho) const;

        std::vector<Real> expiryYears_, tailYears_;
        Real gamma_;
        Handle<YieldTermStructure> yieldCurve_;
        ext::shared_ptr<SwapIndex> swapIndex_;
        Period maxTenor_;

        ext::shared_ptr<BilinearInterpolation> atmInterp_;
        ext::shared_ptr<BilinearInterpolation> betaInterp_;
        ext::shared_ptr<BilinearInterpolation> nuInterp_;
        ext::shared_ptr<BilinearInterpolation> rhoInterp_;
    };


    // -----------------------------------------------------------------------
    // inline / template implementations
    // -----------------------------------------------------------------------

    inline ZabrSwaptionVol::ZabrSwaptionVol(
            const Date& referenceDate,
            const Calendar& calendar,
            const std::vector<Real>& expiryYears,
            const std::vector<Real>& tailYears,
            const Matrix& atmNormalVol,
            const Matrix& beta,
            const Matrix& nu,
            const Matrix& rho,
            Real gamma,
            Handle<YieldTermStructure> yieldCurve,
            ext::shared_ptr<SwapIndex> swapIndex)
    : SwaptionVolatilityStructure(referenceDate, calendar, Following,
                                   Actual365Fixed()),
      expiryYears_(expiryYears), tailYears_(tailYears), gamma_(gamma),
      yieldCurve_(std::move(yieldCurve)),
      swapIndex_(std::move(swapIndex)),
      maxTenor_(static_cast<int>(tailYears.back()) * 12, Months) {

        // Bilinear interpolation on (tail, expiry) — QL convention: x=cols, y=rows
        atmInterp_ = ext::make_shared<BilinearInterpolation>(
            tailYears_.begin(), tailYears_.end(),
            expiryYears_.begin(), expiryYears_.end(), atmNormalVol);
        betaInterp_ = ext::make_shared<BilinearInterpolation>(
            tailYears_.begin(), tailYears_.end(),
            expiryYears_.begin(), expiryYears_.end(), beta);
        nuInterp_ = ext::make_shared<BilinearInterpolation>(
            tailYears_.begin(), tailYears_.end(),
            expiryYears_.begin(), expiryYears_.end(), nu);
        rhoInterp_ = ext::make_shared<BilinearInterpolation>(
            tailYears_.begin(), tailYears_.end(),
            expiryYears_.begin(), expiryYears_.end(), rho);

        for (auto* p : {&*atmInterp_, &*betaInterp_, &*nuInterp_, &*rhoInterp_})
            p->enableExtrapolation();
    }

    inline Real ZabrSwaptionVol::forwardSwapRate(Time optionTime,
                                                  Time swapLength) const {
        Date today = referenceDate();
        Calendar cal = calendar();
        Integer expDays = static_cast<Integer>(optionTime * 365.25 + 0.5);
        Date exDate = cal.advance(today, expDays, Days);
        auto fwdSwap = swapIndex_->clone(
            Period(static_cast<int>(swapLength * 12 + 0.5), Months))
            ->underlyingSwap(exDate);
        fwdSwap->setPricingEngine(
            ext::make_shared<DiscountingSwapEngine>(yieldCurve_));
        return fwdSwap->fairRate();
    }

    inline Real ZabrSwaptionVol::calibrateAlpha(
            Real forward, Real atmNV, Time T,
            Real beta, Real nu, Real rho) const {
        // Initial guess: α ≈ σ_N / F^(1-β)
        Real atmLN = atmNV / forward;
        Real alpha = atmLN * std::pow(forward, 1.0 - beta);
        // Newton iteration on Floch-Kennedy SABR
        for (int iter = 0; iter < 30; ++iter) {
            Real vol = sabrFlochKennedyVolatility(forward, forward,
                                                   std::max(T, 1e-6),
                                                   alpha, beta, nu, rho);
            Real da = alpha * 1e-4;
            Real vol2 = sabrFlochKennedyVolatility(forward, forward,
                                                    std::max(T, 1e-6),
                                                    alpha + da, beta, nu, rho);
            Real dvol = (vol2 - vol) / da;
            if (std::fabs(dvol) < 1e-14) break;
            alpha -= (vol - atmLN) / dvol;
            if (std::fabs(vol - atmLN) < 1e-12) break;
        }
        return std::max(alpha, 1e-10);
    }

    inline ext::shared_ptr<SmileSection>
    ZabrSwaptionVol::smileSectionImpl(Time tExp, Time tTail) const {
        Real atmNV = (*atmInterp_)(tTail, tExp);
        Real beta  = (*betaInterp_)(tTail, tExp);
        Real nu    = (*nuInterp_)(tTail, tExp);
        Real rho   = (*rhoInterp_)(tTail, tExp);

        Real fwd = forwardSwapRate(tExp, tTail);
        Real alpha = calibrateAlpha(fwd, atmNV, tExp, beta, nu, rho);

        std::vector<Real> params = {alpha, beta, nu, rho, gamma_};
        return ext::make_shared<ZabrSmileSection<ZabrLocalVolatility>>(
            tExp, fwd, params);
    }

    inline Volatility ZabrSwaptionVol::volatilityImpl(
            Time optionTime, Time swapLength, Rate strike) const {
        auto section = smileSectionImpl(optionTime, swapLength);
        return section->volatility(strike);
    }

}

#endif
