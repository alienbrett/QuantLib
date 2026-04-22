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

#include <ql/termstructures/volatility/equityfx/pwlpdfblackvolsurface.hpp>
#include <ql/math/interpolations/linearinterpolation.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    // Shared init: validate, build slices, sort, set interpolation.
    // Called from both constructors after members are initialized.
    static void initSlices(
            std::vector<PwlPdfBlackVolSurface::Slice>& slices,
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& kGrids,
            const std::vector<std::vector<Real>>& impliedVols,
            const std::vector<Real>& calFwds,
            const DayCounter& dc) {
        Size N = dates.size();
        bool hasFwds = !calFwds.empty();
        if (hasFwds) {
            QL_REQUIRE(calFwds.size() == N,
                       "calibrationForwards size (" << calFwds.size()
                       << ") must match dates size (" << N << ")");
        }

        slices.resize(N);
        for (Size i = 0; i < N; ++i) {
            Time T = dc.yearFraction(referenceDate, dates[i]);
            QL_REQUIRE(T > 0.0,
                       "date[" << i << "] must be after reference date");
            Size m = kGrids[i].size();
            QL_REQUIRE(impliedVols[i].size() == m,
                       "slice " << i << ": kGrid/IV size mismatch");
            QL_REQUIRE(m >= 2,
                       "slice " << i << " needs >= 2 grid points");

            slices[i].T = T;
            slices[i].calForward = hasFwds ? calFwds[i] : 0.0;
            slices[i].k = kGrids[i];
            slices[i].totalVar.resize(m);
            for (Size j = 0; j < m; ++j) {
                Real v = impliedVols[i][j];
                slices[i].totalVar[j] = v * v * T;
            }
        }

        std::sort(slices.begin(), slices.end(),
                  [](const PwlPdfBlackVolSurface::Slice& a,
                     const PwlPdfBlackVolSurface::Slice& b) {
                      return a.T < b.T;
                  });
        for (Size i = 1; i < N; ++i) {
            QL_REQUIRE(slices[i].T > slices[i-1].T,
                       "non-increasing maturity at slice " << i);
        }
    }

    // =================================================================
    // Constructor: continuous dividends only
    // =================================================================

    PwlPdfBlackVolSurface::PwlPdfBlackVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& kGrids,
            const std::vector<std::vector<Real>>& impliedVols,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc,
            const std::vector<Real>& calibrationForwards)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_()
    {
        Size N = dates.size();
        QL_REQUIRE(kGrids.size() == N && impliedVols.size() == N && N > 0,
                   "dates/kGrids/impliedVols size mismatch or empty");
        QL_REQUIRE(!spot_.empty() && !riskFreeRate_.empty()
                   && !dividendYield_.empty(),
                   "spot/rate/div handles must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);

        initSlices(slices_, referenceDate, dates, kGrids, impliedVols,
                   calibrationForwards, dc);
        setInterpolation<Linear>();
    }

    // =================================================================
    // Constructor: with discrete dividends (matches eSSVI)
    // =================================================================

    PwlPdfBlackVolSurface::PwlPdfBlackVolSurface(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<std::vector<Real>>& kGrids,
            const std::vector<std::vector<Real>>& impliedVols,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
            const DayCounter& dc,
            const std::vector<Real>& calibrationForwards)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_(std::move(dividends))
    {
        Size N = dates.size();
        QL_REQUIRE(kGrids.size() == N && impliedVols.size() == N && N > 0,
                   "dates/kGrids/impliedVols size mismatch or empty");
        QL_REQUIRE(!spot_.empty() && !riskFreeRate_.empty()
                   && !dividendYield_.empty(),
                   "spot/rate/div handles must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);

        initSlices(slices_, referenceDate, dates, kGrids, impliedVols,
                   calibrationForwards, dc);
        setInterpolation<Linear>();
    }

    // =================================================================
    // TermStructure interface
    // =================================================================

    Date PwlPdfBlackVolSurface::maxDate() const {
        return Date::maxDate();
    }

    Real PwlPdfBlackVolSurface::minStrike() const {
        return QL_MIN_REAL;
    }

    Real PwlPdfBlackVolSurface::maxStrike() const {
        return QL_MAX_REAL;
    }

    // =================================================================
    // Forward — matches EssviVolatilityTermStructure::forward exactly
    // =================================================================

    Real PwlPdfBlackVolSurface::forward(Time t) const {
        Real S = spot_->value();
        Real df = riskFreeRate_->discount(t);
        Real dq = dividendYield_->discount(t);

        // Subtract PV of discrete dividends before t
        Real pvDivs = 0.0;
        if (!dividends_.empty()) {
            Date refDate = referenceDate();
            DayCounter dc = dayCounter();
            for (const auto& div : dividends_) {
                Date exDate = div->date();
                if (exDate <= refDate)
                    continue;
                Time tDiv = dc.yearFraction(refDate, exDate);
                if (tDiv < t) {
                    pvDivs += div->amount() * riskFreeRate_->discount(tDiv);
                }
            }
        }

        return (S - pvDivs) * dq / df;
    }

    Real PwlPdfBlackVolSurface::forwardAtSlice(Size idx) const {
        Real F = slices_[idx].calForward;
        if (F > 0.0)
            return F;
        return forward(slices_[idx].T);
    }

    // =================================================================
    // Per-slice interpolation
    // =================================================================

    Real PwlPdfBlackVolSurface::interpolateSlice(Size idx, Real k) const {
        const Slice& sl = slices_[idx];
        Real kc = std::max(k, sl.k.front());
        kc = std::min(kc, sl.k.back());
        Real w = sl.interp(kc, true);
        return std::max(w, 0.0);
    }

    // =================================================================
    // Core: blackVolImpl
    // =================================================================

    Volatility PwlPdfBlackVolSurface::blackVolImpl(Time t, Real strike) const {
        if (t < 1e-14) t = 1e-14;

        Real w;
        Size N = slices_.size();

        if (N == 1 || t <= slices_.front().T) {
            Real fwd = forwardAtSlice(0);
            Real k = std::log(strike / fwd);
            w = interpolateSlice(0, k);
            w = w * t / slices_[0].T;
        } else if (t >= slices_.back().T) {
            Real fwd = forwardAtSlice(N - 1);
            Real k = std::log(strike / fwd);
            w = interpolateSlice(N - 1, k);
            w = w * t / slices_[N - 1].T;
        } else {
            // Between slices — each uses its own calibration forward
            Size hi = 0;
            for (Size i = 1; i < N; ++i) {
                if (slices_[i].T >= t) {
                    hi = i;
                    break;
                }
            }
            Size lo = hi - 1;
            Real T_lo = slices_[lo].T;
            Real T_hi = slices_[hi].T;
            Real alpha = (t - T_lo) / (T_hi - T_lo);

            Real fwd_lo = forwardAtSlice(lo);
            Real fwd_hi = forwardAtSlice(hi);
            Real k_lo = std::log(strike / fwd_lo);
            Real k_hi = std::log(strike / fwd_hi);

            Real w_lo = interpolateSlice(lo, k_lo);
            Real w_hi = interpolateSlice(hi, k_hi);
            w = (1.0 - alpha) * w_lo + alpha * w_hi;
        }

        w = std::max(w, 0.0);
        return std::sqrt(w / t);
    }

} // namespace QuantLib
