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
    // callForwardSlice: moneyness-space polynomial integration
    //
    // CF/F at moneyness xEval:
    //   c(xEval) = integral from xEval to x_max of (u - xEval) q(u) du
    //
    // For full segment [x_l, x_{l+1}] with linear q, h = x_{l+1} - x_l:
    //   coeff of q_l:     h^2/6 + (x_l - xEval) * h/2
    //   coeff of q_{l+1}: h^2/3 + (x_l - xEval) * h/2
    //
    // For partial segment [xEval, x_{l+1}]:
    //   r = x_{l+1} - xEval, t_off = xEval - x_l
    //   coeff of q_l:     r^3 / (6*h)
    //   coeff of q_{l+1}: r^2 * (2*r + 3*t_off) / (6*h)
    // =================================================================

    Real PwlPdfVolSurface::callForwardSlice(Size idx, Real xEval) const {
        const PdfSlice& sl = slices_[idx];
        Size m = sl.x.size();

        if (xEval >= sl.x[m - 1]) return 0.0;
        if (xEval <= sl.x[0]) xEval = sl.x[0];

        // Find segment: x[seg] <= xEval < x[seg+1]
        auto it = std::upper_bound(sl.x.begin(), sl.x.end(), xEval);
        Size seg = (it == sl.x.begin()) ? 0
                 : static_cast<Size>(it - sl.x.begin()) - 1;
        if (seg >= m - 1) seg = m - 2;

        bool onGrid = (std::abs(xEval - sl.x[seg]) < 1e-14);
        Real result = 0.0;

        // Partial first segment [xEval, x[seg+1]]
        if (!onGrid) {
            Real hi = sl.h[seg];
            Real r = sl.x[seg + 1] - xEval;
            Real tOff = xEval - sl.x[seg];
            result += sl.q[seg] * (r * r * r / (6.0 * hi));
            result += sl.q[seg + 1] * (r * r * (2.0 * r + 3.0 * tOff)
                                        / (6.0 * hi));
        }

        // Full segments
        Size startSeg = onGrid ? seg : seg + 1;
        for (Size l = startSeg; l < m - 1; ++l) {
            Real hi = sl.h[l];
            Real d = sl.x[l] - xEval;  // >= 0
            result += sl.q[l] * (hi * hi / 6.0 + d * hi / 2.0);
            result += sl.q[l + 1] * (hi * hi / 3.0 + d * hi / 2.0);
        }

        return std::max(result, 0.0);
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
    // blackVolImpl
    // =================================================================

    Volatility PwlPdfVolSurface::blackVolImpl(Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;
        Size N = slices_.size();
        Real F, callPrice;

        if (N == 1 || t <= slices_.front().T) {
            F = forwardAtSlice(0);
            Real x = strike / F;
            Real c = callForwardSlice(0, x) * t / slices_[0].T;
            callPrice = F * c;
        } else if (t >= slices_.back().T) {
            F = forwardAtSlice(N - 1);
            Real x = strike / F;
            Real c = callForwardSlice(N - 1, x) * t / slices_[N - 1].T;
            callPrice = F * c;
        } else {
            Size hi = 0;
            for (Size i = 1; i < N; ++i)
                if (slices_[i].T >= t) { hi = i; break; }
            Size lo = hi - 1;
            Real alpha = (t - slices_[lo].T) /
                         (slices_[hi].T - slices_[lo].T);

            Real F_lo = forwardAtSlice(lo);
            Real F_hi = forwardAtSlice(hi);
            Real x_lo = strike / F_lo;
            Real x_hi = strike / F_hi;

            Real C_lo = F_lo * callForwardSlice(lo, x_lo);
            Real C_hi = F_hi * callForwardSlice(hi, x_hi);

            callPrice = (1.0 - alpha) * C_lo + alpha * C_hi;
            F = (1.0 - alpha) * F_lo + alpha * F_hi;
        }

        // Edge: deep OTM
        Real intrinsic = std::max(F - strike, 0.0);
        if (callPrice <= intrinsic + 1e-12 * F) {
            try {
                Real cATM = callForwardSlice(0, 1.0);  // x=1 is ATM
                Real F0 = forwardAtSlice(0);
                if (cATM > 1e-12) {
                    Real sd = blackFormulaImpliedStdDevLiRS(
                        Option::Call, F0, F0, F0 * cATM);
                    return std::max(sd / std::sqrt(slices_[0].T), 0.01);
                }
            } catch (...) {}
            return 0.01;
        }
        // Edge: deep ITM
        if (callPrice >= F * 0.9999) {
            try {
                Real cATM = callForwardSlice(0, 1.0);
                Real F0 = forwardAtSlice(0);
                if (cATM > 1e-12 && t > 1e-8) {
                    Real sd = blackFormulaImpliedStdDevLiRS(
                        Option::Call, F0, F0, F0 * cATM);
                    return sd / std::sqrt(slices_[0].T);
                }
            } catch (...) {}
            return 0.20;
        }

        // OTM-side inversion for stability
        Option::Type optType;
        Real price;
        if (strike >= F) {
            optType = Option::Call;
            price = callPrice;
        } else {
            optType = Option::Put;
            price = callPrice - (F - strike);
            if (price <= 0.0) return 0.01;
        }

        try {
            Real stddev = blackFormulaImpliedStdDevLiRS(
                optType, strike, F, price);
            return stddev / std::sqrt(t);
        } catch (...) {
            return 0.20;
        }
    }

} // namespace QuantLib
