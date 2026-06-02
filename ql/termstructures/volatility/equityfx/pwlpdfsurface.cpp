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

    // For each slice, compute boundary vols and wing slopes for
    // linear total-variance extrapolation.  Walks inward from each
    // grid edge to find the two outermost points where Black IV
    // inversion succeeds.  Beyond the boundary:
    //   w(k) = w_boundary + slope × (k - k_boundary)
    // where w = σ²T.  This is C¹ at the boundary → Dupire is
    // well-defined → no FD engine crash.
    static void computeWingVols(
            std::vector<PwlPdfVolSurface::PdfSlice>& slices) {
        for (auto& sl : slices) {
            Size m = sl.x.size();
            if (m < 4) {
                sl.leftWingVol = 0.20;
                sl.rightWingVol = 0.20;
                sl.leftBoundaryK = std::log(sl.x[0]);
                sl.rightBoundaryK = std::log(sl.x[m - 1]);
                sl.leftWingSlope = 0.0;
                sl.rightWingSlope = 0.0;
                continue;
            }

            Real atmVol = sliceVol(sl, 1.0);
            if (atmVol < 0.0) atmVol = 0.20;

            // Right wing: find two outermost valid vols for slope
            sl.rightWingVol = atmVol;
            sl.rightBoundaryK = 0.0;
            sl.rightWingSlope = 0.0;
            {
                Real k1 = 0.0, w1 = 0.0;
                bool found1 = false;
                for (Size i = m - 1; i > 0; --i) {
                    Real vol = sliceVol(sl, sl.x[i]);
                    if (vol > 0.0) {
                        Real k = std::log(sl.x[i]);
                        Real w = vol * vol * sl.T;
                        if (!found1) {
                            sl.rightWingVol = vol;
                            sl.rightBoundaryK = k;
                            k1 = k; w1 = w;
                            found1 = true;
                        } else {
                            // Second point → compute slope
                            if (std::abs(k1 - k) > 1e-10)
                                sl.rightWingSlope = std::max(
                                    (w1 - w) / (k1 - k), 0.0);
                            break;
                        }
                    }
                }
            }

            // Left wing: find two outermost valid vols for slope
            sl.leftWingVol = atmVol;
            sl.leftBoundaryK = 0.0;
            sl.leftWingSlope = 0.0;
            {
                Real k1 = 0.0, w1 = 0.0;
                bool found1 = false;
                for (Size i = 0; i < m; ++i) {
                    Real vol = sliceVol(sl, sl.x[i]);
                    if (vol > 0.0) {
                        Real k = std::log(sl.x[i]);
                        Real w = vol * vol * sl.T;
                        if (!found1) {
                            sl.leftWingVol = vol;
                            sl.leftBoundaryK = k;
                            k1 = k; w1 = w;
                            found1 = true;
                        } else {
                            // Second point → compute slope
                            if (std::abs(k - k1) > 1e-10)
                                sl.leftWingSlope = std::min(
                                    (w - w1) / (k - k1), 0.0);
                            break;
                        }
                    }
                }
            }
        }

        // Pass 2: enforce calendar monotonicity across slices.
        // At any fixed k in the wing, w(k, T) must be non-decreasing
        // in T.  Short-dated slices have narrower grids → steeper
        // extrapolated slopes.  Clamp each slice's wing w so it
        // never exceeds the next slice's w at the same k.
        //
        // Strategy: for each slice i (ascending T), compute w at the
        // previous slice's boundary k.  If w_i < w_{i-1}, flatten
        // the slope so w_i >= w_{i-1} everywhere beyond the boundary.
        if (slices.size() >= 2) {
            for (Size i = 1; i < slices.size(); ++i) {
                auto& prev = slices[i - 1];
                auto& curr = slices[i];

                // Left wing: check at prev's boundary k (most extreme)
                Real kL = prev.leftBoundaryK;
                Real w_prev_L = prev.leftWingVol * prev.leftWingVol * prev.T;
                // w at kL on prev = w_prev_L (it's the boundary)
                // w at kL on curr = curr boundary + curr slope * (kL - curr boundary)
                Real w_curr_L = curr.leftWingVol * curr.leftWingVol * curr.T
                    + curr.leftWingSlope * (kL - curr.leftBoundaryK);
                if (w_curr_L < w_prev_L) {
                    // Flatten curr's left slope so it matches prev at kL
                    Real w_curr_bdy = curr.leftWingVol * curr.leftWingVol * curr.T;
                    if (std::abs(kL - curr.leftBoundaryK) > 1e-10) {
                        curr.leftWingSlope = std::min(
                            (w_prev_L - w_curr_bdy) / (kL - curr.leftBoundaryK),
                            0.0);
                    } else {
                        curr.leftWingSlope = 0.0;
                    }
                }

                // Right wing: check at prev's boundary k
                Real kR = prev.rightBoundaryK;
                Real w_prev_R = prev.rightWingVol * prev.rightWingVol * prev.T;
                Real w_curr_R = curr.rightWingVol * curr.rightWingVol * curr.T
                    + curr.rightWingSlope * (kR - curr.rightBoundaryK);
                if (w_curr_R < w_prev_R) {
                    Real w_curr_bdy = curr.rightWingVol * curr.rightWingVol * curr.T;
                    if (std::abs(kR - curr.rightBoundaryK) > 1e-10) {
                        curr.rightWingSlope = std::max(
                            (w_prev_R - w_curr_bdy) / (kR - curr.rightBoundaryK),
                            0.0);
                    } else {
                        curr.rightWingSlope = 0.0;
                    }
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
    // blackVarianceImpl — strict calendar monotonicity
    // =================================================================
    //
    // blackForwardVariance (non-virtual in BlackVolTermStructure) asserts
    // var(t2) >= var(t1).  Our blackVolImpl enforces calendar monotonicity
    // via a floor scan, but the piecewise nature of the wing extrapolation
    // and IV inversion can produce sub-epsilon variance decreases across
    // slice boundaries.
    //
    // Fix: add a tiny strictly-increasing perturbation ε·t (ε = 1e-14)
    // so var(t+dt) > var(t) numerically, without affecting pricing.
    // ε·1year = 1e-14 → unmeasurable in vol space.

    Real PwlPdfVolSurface::blackVarianceImpl(Time t, Real strike) const {
        Volatility vol = blackVolImpl(t, strike);
        return vol * vol * t + 1e-14 * t;
    }

    // =================================================================
    // blackVolImpl — bilinear in total variance
    // =================================================================

    // Linear total variance wing extrapolation.
    // w(k) = w_boundary + slope × (k - k_boundary), then σ = √(w/T).
    // C¹ continuous → Dupire well-defined → no FD crash.
    static Volatility wingVol(
            const PwlPdfVolSurface::PdfSlice& sl,
            Real k, bool isRight) {
        Real bvol = isRight ? sl.rightWingVol : sl.leftWingVol;
        Real bk   = isRight ? sl.rightBoundaryK : sl.leftBoundaryK;
        Real slope = isRight ? sl.rightWingSlope : sl.leftWingSlope;
        Real w = bvol * bvol * sl.T + slope * (k - bk);
        return std::sqrt(std::max(w, 1e-8) / sl.T);
    }

    // Invert call forward price at a single slice to Black vol.
    static Volatility volAtSlice(
            const PwlPdfVolSurface& self,
            Size idx,
            Real strike) {
        const auto& sl = self.pdfSlices()[idx];
        Real F = (sl.calForward > 0.0) ? sl.calForward : 0.0;
        if (F <= 0.0) return 0.20;
        Real x = strike / F;
        Real k = std::log(x);
        Real callPrice = F * self.callForward(x, sl.T);

        // Wing extrapolation: linear total variance beyond boundary.
        Real intrinsic = std::max(F - strike, 0.0);
        if (callPrice <= intrinsic + 1e-12 * F)
            return wingVol(sl, k, true);    // right wing
        if (callPrice >= F * 0.9999)
            return wingVol(sl, k, false);   // left wing

        // OTM-side inversion for stability
        Option::Type optType;
        Real price;
        if (strike >= F) {
            optType = Option::Call;
            price = callPrice;
        } else {
            optType = Option::Put;
            price = callPrice - (F - strike);
            if (price <= 0.0) return wingVol(sl, k, false);
        }

        try {
            Real stddev = blackFormulaImpliedStdDevLiRS(
                optType, strike, F, price);
            return stddev / std::sqrt(sl.T);
        } catch (...) {
            return wingVol(sl, k, x >= 1.0);
        }
    }

    Volatility PwlPdfVolSurface::blackVolImpl(Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;
        Size N = slices_.size();

        if (N == 1 || t <= slices_.front().T) {
            return volAtSlice(*this, 0, strike);
        }
        if (t >= slices_.back().T) {
            // After last slice: ensure w(t) >= w(last slice)
            Real vol_last = volAtSlice(*this, N - 1, strike);
            Real w_last = vol_last * vol_last * slices_[N - 1].T;
            // Check all earlier slices for calendar monotonicity
            for (Size i = 0; i < N - 1; ++i) {
                Real vi = volAtSlice(*this, i, strike);
                Real wi = vi * vi * slices_[i].T;
                if (wi > w_last) w_last = wi;
            }
            return std::sqrt(std::max(w_last, 0.0) / t);
        }

        // Between slices: bilinear in total variance with calendar floor.
        // Scan ALL slices up to and including hi to find the maximum
        // total variance at this strike.  w(t) must never decrease.
        Size hi = 0;
        for (Size i = 1; i < N; ++i)
            if (slices_[i].T >= t) { hi = i; break; }
        Size lo = hi - 1;

        // Compute w at all slices up to hi, enforce monotonicity
        Real w_lo = 0.0, w_hi = 0.0;
        Real w_floor = 0.0;
        for (Size i = 0; i <= hi; ++i) {
            Real vi = volAtSlice(*this, i, strike);
            Real wi = vi * vi * slices_[i].T;
            if (wi < w_floor) wi = w_floor;
            w_floor = wi;
            if (i == lo) w_lo = wi;
            if (i == hi) w_hi = wi;
        }

        Real alpha = (t - slices_[lo].T) /
                     (slices_[hi].T - slices_[lo].T);
        Real w = (1.0 - alpha) * w_lo + alpha * w_hi;

        return std::sqrt(std::max(w, 0.0) / t);
    }

    // =================================================================
    // Variance helpers and Dupire local vol
    //
    // Dupire in (K, T) space for undiscounted call C:
    //   sigma_loc^2 = [dC/dT + (r-q)*K*dC/dK] / [0.5 * K^2 * p(K)]
    // with dC/dK = -S(K/F), p(K) = q_x(x)/F, x = K/F.
    // =================================================================

    Real PwlPdfVolSurface::totalVariance(Real k, Time t) const {
        if (t < 1e-14) t = 1e-14;
        Real F = forward(t);
        Real K = F * std::exp(k);
        Volatility v = blackVolImpl(t, K);
        return v * v * t;
    }

    Real PwlPdfVolSurface::totalVarianceTimeDerivative(Real k, Time t) const {
        const Real dt = std::max(t * 1e-3, 1e-5);
        Time tLo = std::max(t - dt, 1e-6);
        Time tHi = t + dt;
        Real wLo = totalVariance(k, tLo);
        Real wHi = totalVariance(k, tHi);
        return (wHi - wLo) / (tHi - tLo);
    }

    Real PwlPdfVolSurface::localVol(Real k, Time t) const {
        return std::sqrt(localVariance(k, t));
    }

    Real PwlPdfVolSurface::localVariance(Real k, Time t) const {
        if (t < 1e-14) t = 1e-14;

        const Real volFloor = 0.005;
        const Real volCap = 5.0;
        const Real pdfFloor = 1e-8;

        Real F = forward(t);
        Real K = F * std::exp(k);
        Real x = K / F;

        // Black vol used both as floor and as wing fallback.
        Real bvol = std::max(blackVolImpl(t, K), volFloor);
        Real floorVar = bvol * bvol;

        // Denominator: 0.5 * K^2 * p(K) = 0.5 * K^2 * q_x / F.
        Real q_val = pdfValue(x, t);
        if (q_val < pdfFloor)
            return floorVar;
        Real denominator = 0.5 * K * K * q_val / F;

        // dC/dT at fixed K via two-slice FD on undiscounted call forwards.
        Size N = slices_.size();
        if (N == 1)
            return floorVar;

        Size lo, hi;
        if (t <= slices_.front().T) {
            lo = 0; hi = 1;
        } else if (t >= slices_.back().T) {
            lo = N - 2; hi = N - 1;
        } else {
            hi = 0;
            for (Size i = 1; i < N; ++i)
                if (slices_[i].T >= t) { hi = i; break; }
            lo = hi - 1;
        }

        Real T_lo = slices_[lo].T;
        Real T_hi = slices_[hi].T;
        Real S = spot_->value();
        Real F_lo = slices_[lo].calForward > 0.0
            ? slices_[lo].calForward
            : S * dividendYield_->discount(T_lo, true)
                / riskFreeRate_->discount(T_lo, true);
        Real F_hi = slices_[hi].calForward > 0.0
            ? slices_[hi].calForward
            : S * dividendYield_->discount(T_hi, true)
                / riskFreeRate_->discount(T_hi, true);

        Real x_lo = K / F_lo;
        Real x_hi = K / F_hi;
        Real C_lo = F_lo * callForward(x_lo, T_lo);
        Real C_hi = F_hi * callForward(x_hi, T_hi);
        Real dCdT = (C_hi - C_lo) / (T_hi - T_lo);
        if (dCdT < 0.0) dCdT = 0.0;

        // (r-q) at t — instantaneous forward rates from yield curves.
        Real r_inst = riskFreeRate_->forwardRate(
            t, t + 1e-4, Continuous, NoFrequency).rate();
        Real q_inst = dividendYield_->forwardRate(
            t, t + 1e-4, Continuous, NoFrequency).rate();
        Real drift = r_inst - q_inst;

        Real survK = survivalFunction(x, t);
        Real dCdK = -survK;

        Real numerator = dCdT + drift * K * dCdK;
        if (numerator < 0.0) numerator = 0.0;

        Real lv = numerator / denominator;
        if (lv <= 0.0 || !std::isfinite(lv))
            return floorVar;

        // Relative cap: local vol ≤ 2× black vol → local var ≤ 4 × black var.
        Real capVar = std::min(4.0 * bvol * bvol, volCap * volCap);
        return std::min(std::max(lv, volFloor * volFloor), capVar);
    }

} // namespace QuantLib
