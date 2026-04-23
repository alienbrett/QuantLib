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

#include <ql/termstructures/volatility/equityfx/pwlpdfsurface.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    // =================================================================
    // Static helpers for wing vol computation
    // =================================================================

    // Compute call forward CF/F from a PdfSlice at moneyness xEval.
    // Same polynomial integration as the member function, but takes
    // PdfSlice directly so it can be called before construction.
    static Real callForwardOnSlice(
            const PwlPdfVolSurface::PdfSlice& sl, Real xEval) {
        Size m = sl.x.size();
        if (xEval >= sl.x[m - 1]) return 0.0;
        if (xEval <= sl.x[0]) xEval = sl.x[0];

        auto it = std::upper_bound(sl.x.begin(), sl.x.end(), xEval);
        Size seg = (it == sl.x.begin()) ? 0
                 : static_cast<Size>(it - sl.x.begin()) - 1;
        if (seg >= m - 1) seg = m - 2;

        bool onGrid = (std::abs(xEval - sl.x[seg]) < 1e-14);
        Real result = 0.0;

        if (!onGrid) {
            Real hi = sl.h[seg];
            Real r = sl.x[seg + 1] - xEval;
            Real tOff = xEval - sl.x[seg];
            result += sl.q[seg] * (r * r * r / (6.0 * hi));
            result += sl.q[seg + 1] * (r * r * (2.0 * r + 3.0 * tOff)
                                        / (6.0 * hi));
        }

        Size startSeg = onGrid ? seg : seg + 1;
        for (Size l = startSeg; l < m - 1; ++l) {
            Real hi = sl.h[l];
            Real d = sl.x[l] - xEval;
            result += sl.q[l] * (hi * hi / 6.0 + d * hi / 2.0);
            result += sl.q[l + 1] * (hi * hi / 3.0 + d * hi / 2.0);
        }

        return std::max(result, 0.0);
    }

    // Compute Black vol at a grid point on a single slice.
    // Returns negative if IV inversion fails (no valid vol).
    static Volatility sliceVol(
            const PwlPdfVolSurface::PdfSlice& sl, Real xEval) {
        Real F = sl.calForward;
        if (F <= 0.0) return -1.0;

        Real K = F * xEval;
        Real callFwd = callForwardOnSlice(sl, xEval);
        Real callPrice = F * callFwd;

        Real intrinsic = std::max(F - K, 0.0);
        if (callPrice <= intrinsic + 1e-12 * F) return -1.0;
        if (callPrice >= F * 0.9999) return -1.0;

        Option::Type optType;
        Real price;
        if (K >= F) {
            optType = Option::Call;
            price = callPrice;
        } else {
            optType = Option::Put;
            price = callPrice - (F - K);
            if (price <= 0.0) return -1.0;
        }

        try {
            Real stddev = blackFormulaImpliedStdDevLiRS(
                optType, K, F, price);
            Real vol = stddev / std::sqrt(sl.T);
            if (vol > 0.005 && vol < 5.0) return vol;
            return -1.0;
        } catch (...) {
            return -1.0;
        }
    }

    // For each slice, compute boundary vols for wing extrapolation.
    // Walks inward from each grid edge to find the outermost point
    // where Black IV inversion succeeds.  Strikes beyond this get
    // flat vol = boundary vol.  This keeps local vol finite and
    // well-defined (flat black vol → local vol = black vol).
    static void computeWingVols(
            std::vector<PwlPdfVolSurface::PdfSlice>& slices) {
        for (auto& sl : slices) {
            Size m = sl.x.size();
            if (m < 4) {
                sl.leftWingVol = 0.20;
                sl.rightWingVol = 0.20;
                continue;
            }

            // ATM vol as ultimate fallback
            Real atmVol = sliceVol(sl, 1.0);
            if (atmVol < 0.0) atmVol = 0.20;

            // Right wing: walk inward from right edge
            sl.rightWingVol = atmVol;
            for (Size i = m - 1; i > 0; --i) {
                Real vol = sliceVol(sl, sl.x[i]);
                if (vol > 0.0) {
                    sl.rightWingVol = vol;
                    break;
                }
            }

            // Left wing: walk inward from left edge
            sl.leftWingVol = atmVol;
            for (Size i = 0; i < m; ++i) {
                Real vol = sliceVol(sl, sl.x[i]);
                if (vol > 0.0) {
                    sl.leftWingVol = vol;
                    break;
                }
            }
        }
    }

    // =================================================================
    // Init: validate, build slices, sort, precompute moneyness grid
    // =================================================================

    static void initPdfSlices(
            std::vector<PwlPdfVolSurface::PdfSlice>& slices,
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& xGrids,
            const std::vector<std::vector<Real>>& qValues,
            const std::vector<Real>& calFwds,
            const DayCounter& dc) {

        Size N = dates.size();
        QL_REQUIRE(xGrids.size() == N && qValues.size() == N && N > 0,
                   "dates/xGrids/qValues size mismatch or empty");
        QL_REQUIRE(calFwds.size() == N,
                   "calibrationForwards size (" << calFwds.size()
                   << ") must match dates size (" << N << ")");

        slices.resize(N);
        for (Size i = 0; i < N; ++i) {
            Time T = dc.yearFraction(referenceDate, dates[i]);
            QL_REQUIRE(T > 0.0, "date[" << i << "] must be after reference date");

            Size m = xGrids[i].size();
            QL_REQUIRE(qValues[i].size() == m,
                       "slice " << i << ": xGrid/qValues size mismatch");
            QL_REQUIRE(m >= 2, "slice " << i << " needs >= 2 grid points");

            auto& sl = slices[i];
            sl.T = T;
            sl.calForward = calFwds[i];
            sl.x = xGrids[i];
            sl.q = qValues[i];

            // Precompute moneyness spacings
            sl.h.resize(m - 1);
            for (Size j = 0; j < m - 1; ++j)
                sl.h[j] = sl.x[j + 1] - sl.x[j];
        }

        std::sort(slices.begin(), slices.end(),
                  [](const PwlPdfVolSurface::PdfSlice& a,
                     const PwlPdfVolSurface::PdfSlice& b) {
                      return a.T < b.T;
                  });

        for (Size i = 1; i < N; ++i) {
            QL_REQUIRE(slices[i].T > slices[i-1].T,
                       "non-increasing maturity at slice " << i);
        }
    }

    // =================================================================
    // Constructors
    // =================================================================

    PwlPdfVolSurface::PwlPdfVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& xGrids,
            const std::vector<std::vector<Real>>& qValues,
            const std::vector<Real>& calibrationForwards,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_()
    {
        QL_REQUIRE(!spot_.empty() && !riskFreeRate_.empty()
                   && !dividendYield_.empty(),
                   "spot/rate/div handles must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
        initPdfSlices(slices_, referenceDate, dates, xGrids, qValues,
                      calibrationForwards, dc);
        computeWingVols(slices_);
    }

    PwlPdfVolSurface::PwlPdfVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& xGrids,
            const std::vector<std::vector<Real>>& qValues,
            const std::vector<Real>& calibrationForwards,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_(std::move(dividends))
    {
        QL_REQUIRE(!spot_.empty() && !riskFreeRate_.empty()
                   && !dividendYield_.empty(),
                   "spot/rate/div handles must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
        initPdfSlices(slices_, referenceDate, dates, xGrids, qValues,
                      calibrationForwards, dc);
        computeWingVols(slices_);
    }

    // =================================================================
    // TermStructure interface
    // =================================================================

    Date PwlPdfVolSurface::maxDate() const { return Date::maxDate(); }
    Real PwlPdfVolSurface::minStrike() const { return QL_MIN_REAL; }
    Real PwlPdfVolSurface::maxStrike() const { return QL_MAX_REAL; }

    // =================================================================
    // Forward
    // =================================================================

    Real PwlPdfVolSurface::forward(Time t) const {
        Real S = spot_->value();
        Real df = riskFreeRate_->discount(t);
        Real dq = dividendYield_->discount(t);
        Real pvDivs = 0.0;
        if (!dividends_.empty()) {
            Date refDate = referenceDate();
            DayCounter dc = dayCounter();
            for (const auto& div : dividends_) {
                Date exDate = div->date();
                if (exDate <= refDate) continue;
                Time tDiv = dc.yearFraction(refDate, exDate);
                if (tDiv < t)
                    pvDivs += div->amount() * riskFreeRate_->discount(tDiv);
            }
        }
        return (S - pvDivs) * dq / df;
    }

    Real PwlPdfVolSurface::forwardAtSlice(Size idx) const {
        Real F = slices_[idx].calForward;
        return (F > 0.0) ? F : forward(slices_[idx].T);
    }

    // =================================================================
    // callForwardSlice: delegates to static callForwardOnSlice
    // =================================================================

    Real PwlPdfVolSurface::callForwardSlice(Size idx, Real xEval) const {
        return callForwardOnSlice(slices_[idx], xEval);
    }

    // =================================================================
    // interpolatePdf: linear interpolation on moneyness grid
    // =================================================================

    Real PwlPdfVolSurface::interpolatePdf(Size idx, Real xEval) const {
        const PdfSlice& sl = slices_[idx];
        Size m = sl.x.size();
        if (xEval <= sl.x[0] || xEval >= sl.x[m - 1]) return 0.0;

        auto it = std::upper_bound(sl.x.begin(), sl.x.end(), xEval);
        Size j = static_cast<Size>(it - sl.x.begin()) - 1;
        if (j >= m - 1) j = m - 2;

        Real alpha = (xEval - sl.x[j]) / sl.h[j];
        return sl.q[j] * (1.0 - alpha) + sl.q[j + 1] * alpha;
    }

    // =================================================================
    // survivalSlice: S(x) = integral from x to x_max of q(u) du
    // =================================================================

    Real PwlPdfVolSurface::survivalSlice(Size idx, Real xEval) const {
        const PdfSlice& sl = slices_[idx];
        Size m = sl.x.size();
        if (xEval >= sl.x[m - 1]) return 0.0;
        if (xEval <= sl.x[0]) {
            Real total = 0.0;
            for (Size l = 0; l < m - 1; ++l)
                total += 0.5 * (sl.q[l] + sl.q[l + 1]) * sl.h[l];
            return total;
        }

        auto it = std::upper_bound(sl.x.begin(), sl.x.end(), xEval);
        Size seg = static_cast<Size>(it - sl.x.begin()) - 1;
        if (seg >= m - 1) seg = m - 2;

        Real result = 0.0;
        // Partial first segment
        Real alpha = (xEval - sl.x[seg]) / sl.h[seg];
        Real qAtX = sl.q[seg] * (1.0 - alpha) + sl.q[seg + 1] * alpha;
        Real r = sl.x[seg + 1] - xEval;
        result += 0.5 * (qAtX + sl.q[seg + 1]) * r;

        // Full segments
        for (Size l = seg + 1; l < m - 1; ++l)
            result += 0.5 * (sl.q[l] + sl.q[l + 1]) * sl.h[l];
        return result;
    }

    // =================================================================
    // callForward / pdfValue / survivalFunction: time-interpolated
    // =================================================================

    Real PwlPdfVolSurface::callForward(Real x, Time t) const {
        Size N = slices_.size();
        if (t < 1e-14) t = 1e-14;
        if (N == 1 || t <= slices_.front().T)
            return callForwardSlice(0, x) * t / slices_[0].T;
        if (t >= slices_.back().T)
            return callForwardSlice(N - 1, x) * t / slices_[N - 1].T;

        Size hi = 0;
        for (Size i = 1; i < N; ++i)
            if (slices_[i].T >= t) { hi = i; break; }
        Size lo = hi - 1;
        Real alpha = (t - slices_[lo].T) / (slices_[hi].T - slices_[lo].T);
        return (1.0 - alpha) * callForwardSlice(lo, x)
               + alpha * callForwardSlice(hi, x);
    }

    Real PwlPdfVolSurface::pdfValue(Real x, Time t) const {
        Size N = slices_.size();
        if (t < 1e-14) t = 1e-14;
        if (N == 1 || t <= slices_.front().T) return interpolatePdf(0, x);
        if (t >= slices_.back().T) return interpolatePdf(N - 1, x);

        Size hi = 0;
        for (Size i = 1; i < N; ++i)
            if (slices_[i].T >= t) { hi = i; break; }
        Size lo = hi - 1;
        Real alpha = (t - slices_[lo].T) / (slices_[hi].T - slices_[lo].T);
        return (1.0 - alpha) * interpolatePdf(lo, x)
               + alpha * interpolatePdf(hi, x);
    }

    Real PwlPdfVolSurface::survivalFunction(Real x, Time t) const {
        Size N = slices_.size();
        if (t < 1e-14) t = 1e-14;
        if (N == 1 || t <= slices_.front().T) return survivalSlice(0, x);
        if (t >= slices_.back().T) return survivalSlice(N - 1, x);

        Size hi = 0;
        for (Size i = 1; i < N; ++i)
            if (slices_[i].T >= t) { hi = i; break; }
        Size lo = hi - 1;
        Real alpha = (t - slices_[lo].T) / (slices_[hi].T - slices_[lo].T);
        return (1.0 - alpha) * survivalSlice(lo, x)
               + alpha * survivalSlice(hi, x);
    }

    // =================================================================
    // blackVolImpl — bilinear in total variance
    // =================================================================

    // Invert call forward price at a single slice to Black vol.
    static Volatility volAtSlice(
            const PwlPdfVolSurface& self,
            Size idx,
            Real strike) {
        const auto& sl = self.pdfSlices()[idx];
        Real F = (sl.calForward > 0.0) ? sl.calForward : 0.0;
        if (F <= 0.0) return 0.20;
        Real x = strike / F;
        Real callPrice = F * self.callForward(x, sl.T);

        // Wing extrapolation: use precomputed boundary vols
        // instead of hardcoded constants.  Flat vol in wings →
        // finite local vol (σ_loc = σ_BS for constant BS vol).
        Real intrinsic = std::max(F - strike, 0.0);
        if (callPrice <= intrinsic + 1e-12 * F)
            return sl.rightWingVol;   // deep OTM call / right wing
        if (callPrice >= F * 0.9999)
            return sl.leftWingVol;    // deep ITM call / left wing

        // OTM-side inversion for stability
        Option::Type optType;
        Real price;
        if (strike >= F) {
            optType = Option::Call;
            price = callPrice;
        } else {
            optType = Option::Put;
            price = callPrice - (F - strike);
            if (price <= 0.0) return sl.leftWingVol;
        }

        try {
            Real stddev = blackFormulaImpliedStdDevLiRS(
                optType, strike, F, price);
            return stddev / std::sqrt(sl.T);
        } catch (...) {
            return (x >= 1.0) ? sl.rightWingVol : sl.leftWingVol;
        }
    }

    Volatility PwlPdfVolSurface::blackVolImpl(Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;
        Size N = slices_.size();

        if (N == 1 || t <= slices_.front().T) {
            // Before/at first slice: flat vol extrapolation
            return volAtSlice(*this, 0, strike);
        }
        if (t >= slices_.back().T) {
            // After last slice: flat vol extrapolation
            return volAtSlice(*this, N - 1, strike);
        }

        // Between slices: bilinear in total variance
        Size hi = 0;
        for (Size i = 1; i < N; ++i)
            if (slices_[i].T >= t) { hi = i; break; }
        Size lo = hi - 1;
        Real alpha = (t - slices_[lo].T) /
                     (slices_[hi].T - slices_[lo].T);

        Real vol_lo = volAtSlice(*this, lo, strike);
        Real vol_hi = volAtSlice(*this, hi, strike);

        Real w_lo = vol_lo * vol_lo * slices_[lo].T;
        Real w_hi = vol_hi * vol_hi * slices_[hi].T;
        Real w = (1.0 - alpha) * w_lo + alpha * w_hi;

        return std::sqrt(std::max(w, 0.0) / t);
    }

} // namespace QuantLib
