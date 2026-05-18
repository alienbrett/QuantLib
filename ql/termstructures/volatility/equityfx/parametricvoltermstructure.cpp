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

#include <ql/termstructures/volatility/equityfx/parametricvoltermstructure.hpp>
#include <ql/cashflows/dividend.hpp>
#include <ql/errors.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace QuantLib {

    // ===== ParametricVolShape default d2fdz2 / dfdParams (central FD) ========

    Real ParametricVolShape::d2fdz2(Real z,
                                    const std::vector<Real>& params,
                                    Real h) const {
        // Central FD on dfdz.  Goes through the vtable so Python
        // proxies (which override dfdz) work transparently.
        return (dfdz(z + h, params) - dfdz(z - h, params)) / (2.0 * h);
    }

    std::vector<Real> ParametricVolShape::dfdParams(
            Real z, const std::vector<Real>& params) const {
        std::vector<Real> out(params.size());
        std::vector<Real> p = params;
        for (std::size_t i = 0; i < params.size(); ++i) {
            Real h = std::max(std::abs(params[i]) * 1e-5, 1e-7);
            p[i] = params[i] + h;
            Real f_plus = f(z, p);
            p[i] = params[i] - h;
            Real f_minus = f(z, p);
            p[i] = params[i];
            out[i] = (f_plus - f_minus) / (2.0 * h);
        }
        return out;
    }

    // ===== S3Shape ============================================================

    Real S3Shape::f(Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 2,
                   "S3Shape: expected 2 params (s2, c2), got " << p.size());
        Real s2 = p[0], c2 = p[1];
        QL_REQUIRE(c2 >= 0.0,
                   "S3Shape: c2 (=" << c2 << ") must be >= 0 for positivity");
        Real a = 1.0 + s2 * z;
        Real radicand = 0.25 * a * a + 0.5 * c2 * z * z;
        return 0.5 * a + std::sqrt(radicand);
    }

    Real S3Shape::dfdz(Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 2,
                   "S3Shape: expected 2 params (s2, c2), got " << p.size());
        Real s2 = p[0], c2 = p[1];
        Real a = 1.0 + s2 * z;
        Real radicand = 0.25 * a * a + 0.5 * c2 * z * z;
        Real root = std::sqrt(radicand);
        // d/dz [0.5 a]       = 0.5 s2
        // d/dz [sqrt(rad)]   = (0.5 a s2 + c2 z) / (2 root)
        return 0.5 * s2 + (0.5 * a * s2 + c2 * z) / (2.0 * root);
    }

    Real S3Shape::d2fdz2(Real z, const std::vector<Real>& p, Real /*h*/) const {
        QL_REQUIRE(p.size() == 2,
                   "S3Shape: expected 2 params (s2, c2), got " << p.size());
        Real s2 = p[0], c2 = p[1];
        Real a = 1.0 + s2 * z;
        Real R  = 0.25 * a * a + 0.5 * c2 * z * z;
        Real r  = std::sqrt(R);
        Real Rp  = 0.5 * a * s2 + c2 * z;       // R'
        Real Rpp = 0.5 * s2 * s2 + c2;          // R''
        // 0.5·a is linear in z so contributes nothing to f''.
        // f''(z) = R''/(2 r) − R'² / (4 r³)
        return Rpp / (2.0 * r) - (Rp * Rp) / (4.0 * r * R);
    }

    std::vector<Real> S3Shape::dfdParams(
            Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 2,
                   "S3Shape: expected 2 params (s2, c2), got " << p.size());
        Real s2 = p[0], c2 = p[1];
        Real a = 1.0 + s2 * z;
        Real R = 0.25 * a * a + 0.5 * c2 * z * z;
        Real sqrtR = std::sqrt(R);
        // ∂f/∂s2 = 0.5·z + (0.5·a·z) / (2·sqrtR) = 0.5·z + 0.25·a·z/sqrtR
        // ∂f/∂c2 = (0.5·z²) / (2·sqrtR) = z²/(4·sqrtR)
        Real ds2 = 0.5 * z + 0.25 * a * z / sqrtR;
        Real dc2 = z * z / (4.0 * sqrtR);
        return {ds2, dc2};
    }


    // ===== JWShape ============================================================

    namespace {
        // S3-with-given-c kernel, used by both halves of JW.
        inline Real s3KernelF(Real z, Real s2, Real c) {
            Real a = 1.0 + s2 * z;
            Real radicand = 0.25 * a * a + 0.5 * c * z * z;
            return 0.5 * a + std::sqrt(radicand);
        }
        inline Real s3KernelDfdz(Real z, Real s2, Real c) {
            Real a = 1.0 + s2 * z;
            Real radicand = 0.25 * a * a + 0.5 * c * z * z;
            Real root = std::sqrt(radicand);
            return 0.5 * s2 + (0.5 * a * s2 + c * z) / (2.0 * root);
        }
        inline Real s3KernelD2fdz2(Real z, Real s2, Real c) {
            Real a = 1.0 + s2 * z;
            Real R  = 0.25 * a * a + 0.5 * c * z * z;
            Real r  = std::sqrt(R);
            Real Rp  = 0.5 * a * s2 + c * z;
            Real Rpp = 0.5 * s2 * s2 + c;
            return Rpp / (2.0 * r) - (Rp * Rp) / (4.0 * r * R);
        }
    }

    Real JWShape::f(Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 3,
                   "JWShape: expected 3 params (s2, c_minus, c_plus), got "
                   << p.size());
        Real s2 = p[0], cm = p[1], cp = p[2];
        QL_REQUIRE(cm >= 0.0 && cp >= 0.0,
                   "JWShape: c_minus (" << cm << ") and c_plus (" << cp
                   << ") must be >= 0");
        return s3KernelF(z, s2, (z < 0.0 ? cm : cp));
    }

    Real JWShape::dfdz(Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 3,
                   "JWShape: expected 3 params (s2, c_minus, c_plus), got "
                   << p.size());
        Real s2 = p[0], cm = p[1], cp = p[2];
        return s3KernelDfdz(z, s2, (z < 0.0 ? cm : cp));
    }

    Real JWShape::d2fdz2(Real z, const std::vector<Real>& p, Real /*h*/) const {
        QL_REQUIRE(p.size() == 3,
                   "JWShape: expected 3 params (s2, c_minus, c_plus), got "
                   << p.size());
        Real s2 = p[0], cm = p[1], cp = p[2];
        // Piecewise-S3.  f''(0) is discontinuous when cm != cp; right-
        // continuous convention matches f / dfdz at z = 0 (uses c_plus).
        return s3KernelD2fdz2(z, s2, (z < 0.0 ? cm : cp));
    }

    std::vector<Real> JWShape::dfdParams(
            Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 3,
                   "JWShape: expected 3 params (s2, c_minus, c_plus), got "
                   << p.size());
        Real s2 = p[0], cm = p[1], cp = p[2];
        Real c_eff = (z < 0.0) ? cm : cp;
        Real a = 1.0 + s2 * z;
        Real R = 0.25 * a * a + 0.5 * c_eff * z * z;
        Real sqrtR = std::sqrt(R);
        Real ds2 = 0.5 * z + 0.25 * a * z / sqrtR;
        Real dc = z * z / (4.0 * sqrtR);
        if (z < 0.0)
            return {ds2, dc, 0.0};
        else
            return {ds2, 0.0, dc};
    }


    // ===== K7Shape ============================================================

    namespace {
        // Fixed blend widths.  Could be made per-slice params in a future K9.
        constexpr Real K7_MU_1 = 0.5;   // inner-wing knot
        constexpr Real K7_MU_2 = 2.0;   // outer-wing knot

        struct K7Blend {
            Real T1, T2;             // tanh^2(z/mu_k)
            Real T1p, T2p;           // d/dz tanh^2(z/mu_k)
            Real T1pp, T2pp;         // d^2/dz^2 tanh^2(z/mu_k)
        };

        inline K7Blend k7BlendAt(Real z) {
            Real u1 = z / K7_MU_1;
            Real u2 = z / K7_MU_2;
            Real th1 = std::tanh(u1);
            Real th2 = std::tanh(u2);
            Real sech2_1 = 1.0 - th1 * th1;
            Real sech2_2 = 1.0 - th2 * th2;
            K7Blend b;
            b.T1 = th1 * th1;
            b.T2 = th2 * th2;
            // d/dz tanh^2(z/mu) = (2/mu) * tanh(z/mu) * sech^2(z/mu)
            b.T1p = (2.0 / K7_MU_1) * th1 * sech2_1;
            b.T2p = (2.0 / K7_MU_2) * th2 * sech2_2;
            // d^2/dz^2 tanh^2(z/mu) = (2/mu^2) * sech^2(z/mu) * (1 - 3*tanh^2)
            b.T1pp = (2.0 / (K7_MU_1 * K7_MU_1)) * sech2_1
                     * (1.0 - 3.0 * b.T1);
            b.T2pp = (2.0 / (K7_MU_2 * K7_MU_2)) * sech2_2
                     * (1.0 - 3.0 * b.T2);
            return b;
        }

        inline void k7CEff(Real z,
                           Real c, Real cmm, Real cmp, Real cfm, Real cfp,
                           Real& c_eff, Real& dc_dz, Real& d2c_dz2) {
            K7Blend b = k7BlendAt(z);
            Real c_mid = (z < 0.0) ? cmm : cmp;
            Real c_far = (z < 0.0) ? cfm : cfp;
            c_eff   = c + (c_mid - c) * b.T1 + (c_far - c_mid) * b.T2;
            dc_dz   = (c_mid - c) * b.T1p + (c_far - c_mid) * b.T2p;
            d2c_dz2 = (c_mid - c) * b.T1pp + (c_far - c_mid) * b.T2pp;
        }
    }

    Real K7Shape::f(Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 6,
                   "K7Shape: expected 6 params "
                   "(s, c, c_mid_minus, c_mid_plus, c_far_minus, c_far_plus), got "
                   << p.size());
        Real s = p[0], c = p[1];
        Real cmm = p[2], cmp = p[3];
        Real cfm = p[4], cfp = p[5];
        QL_REQUIRE(c >= 0.0 && cmm >= 0.0 && cmp >= 0.0
                       && cfm >= 0.0 && cfp >= 0.0,
                   "K7Shape: c, c_mid_*, c_far_* must all be >= 0");
        Real a = 1.0 + s * z;
        Real c_eff, dc_dz, d2c_dz2;
        k7CEff(z, c, cmm, cmp, cfm, cfp, c_eff, dc_dz, d2c_dz2);
        Real R = 0.25 * a * a + 0.5 * c_eff * z * z;
        return 0.5 * a + std::sqrt(std::max(R, 0.0));
    }

    Real K7Shape::dfdz(Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 6, "K7Shape: expected 6 params");
        Real s = p[0], c = p[1];
        Real cmm = p[2], cmp = p[3];
        Real cfm = p[4], cfp = p[5];
        Real a = 1.0 + s * z;
        Real c_eff, dc_dz, d2c_dz2;
        k7CEff(z, c, cmm, cmp, cfm, cfp, c_eff, dc_dz, d2c_dz2);
        Real R = 0.25 * a * a + 0.5 * c_eff * z * z;
        if (R <= 0.0)
            return 0.5 * s;
        // R' = (1/2) s (1 + s z) + (1/2) c_eff' z^2 + c_eff z
        Real Rp = 0.5 * s * a + 0.5 * dc_dz * z * z + c_eff * z;
        return 0.5 * s + Rp / (2.0 * std::sqrt(R));
    }

    Real K7Shape::d2fdz2(Real z, const std::vector<Real>& p,
                         Real /*h*/) const {
        QL_REQUIRE(p.size() == 6, "K7Shape: expected 6 params");
        Real s = p[0], c = p[1];
        Real cmm = p[2], cmp = p[3];
        Real cfm = p[4], cfp = p[5];
        Real a = 1.0 + s * z;
        Real c_eff, dc_dz, d2c_dz2;
        k7CEff(z, c, cmm, cmp, cfm, cfp, c_eff, dc_dz, d2c_dz2);
        Real R = 0.25 * a * a + 0.5 * c_eff * z * z;
        if (R <= 0.0)
            return 0.0;
        Real r = std::sqrt(R);
        Real Rp = 0.5 * s * a + 0.5 * dc_dz * z * z + c_eff * z;
        // R'' = (1/2) s^2 + (1/2) c_eff'' z^2 + 2 c_eff' z + c_eff
        Real Rpp = 0.5 * s * s + 0.5 * d2c_dz2 * z * z
                   + 2.0 * dc_dz * z + c_eff;
        // f' = (1/2) s + R'/(2 sqrt(R))   (linear-in-z half drops out of f'')
        // f'' = R''/(2 r) - R'^2 / (4 r^3)
        return Rpp / (2.0 * r) - (Rp * Rp) / (4.0 * r * R);
    }

    std::vector<Real> K7Shape::dfdParams(
            Real z, const std::vector<Real>& p) const {
        QL_REQUIRE(p.size() == 6, "K7Shape: expected 6 params");
        Real s = p[0], c = p[1];
        Real cmm = p[2], cmp = p[3];
        Real cfm = p[4], cfp = p[5];
        Real a = 1.0 + s * z;
        Real c_eff, dc_dz, d2c_dz2;
        k7CEff(z, c, cmm, cmp, cfm, cfp, c_eff, dc_dz, d2c_dz2);
        Real R = 0.25 * a * a + 0.5 * c_eff * z * z;
        Real sqrtR = (R > 0.0) ? std::sqrt(R) : 1.0;  // guarded; R>0 in fit
        K7Blend b = k7BlendAt(z);

        // f = (1/2) a + sqrt(R), so df/dq = (dR/dq) / (2 sqrt(R)) for any
        // shape param q != s (where the (1/2) a term also contributes).
        // c_eff = c (1 - T1) + c_mid (T1 - T2) + c_far T2
        // dR/dq = (1/2) (dc_eff/dq) z^2

        // d/ds : also contributes via (1/2) a
        Real ds = 0.5 * z + 0.25 * a * z / sqrtR;

        // d/dc
        Real dc = 0.5 * z * z * (1.0 - b.T1) / (2.0 * sqrtR);

        // d/dc_mid_minus, d/dc_mid_plus : (T1 - T2) on respective side
        Real dc_mid_factor = 0.5 * z * z * (b.T1 - b.T2) / (2.0 * sqrtR);
        Real dcmm = (z < 0.0) ? dc_mid_factor : 0.0;
        Real dcmp = (z >= 0.0) ? dc_mid_factor : 0.0;

        // d/dc_far_minus, d/dc_far_plus : T2 on respective side
        Real dc_far_factor = 0.5 * z * z * b.T2 / (2.0 * sqrtR);
        Real dcfm = (z < 0.0) ? dc_far_factor : 0.0;
        Real dcfp = (z >= 0.0) ? dc_far_factor : 0.0;

        return {ds, dc, dcmm, dcmp, dcfm, dcfp};
    }


    // ===== ParametricVolTermStructure =========================================

    namespace {

        std::vector<Time> datesToTimes(const Date& refDate,
                                       const std::vector<Date>& dates,
                                       const DayCounter& dc) {
            std::vector<Time> times(dates.size());
            for (Size i = 0; i < dates.size(); ++i) {
                times[i] = dc.yearFraction(refDate, dates[i]);
                QL_REQUIRE(times[i] > 0.0,
                           "date[" << i << "] (" << dates[i]
                           << ") must be after reference date ("
                           << refDate << ")");
                if (i > 0)
                    QL_REQUIRE(times[i] > times[i-1],
                               "dates must be strictly increasing");
            }
            return times;
        }

        void validateSlices(const std::vector<Time>& T,
                            const std::vector<ParametricVolSlice>& slices) {
            QL_REQUIRE(T.size() == slices.size(),
                       "dates size (" << T.size()
                       << ") must match slices size (" << slices.size() << ")");
            QL_REQUIRE(!slices.empty(), "at least one slice required");
            for (Size i = 0; i < slices.size(); ++i) {
                QL_REQUIRE(slices[i].atmIv > 0.0,
                           "slice[" << i << "] atmIv must be > 0");
            }
        }

    }  // namespace

    ParametricVolTermStructure::ParametricVolTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<ParametricVolSlice>& slices,
            ext::shared_ptr<ParametricVolShape> shape,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          T_(datesToTimes(referenceDate, dates, dc)),
          slices_(slices),
          shape_(std::move(shape)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_()
    {
        validateSlices(T_, slices_);
        QL_REQUIRE(shape_, "shape must not be null");
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    ParametricVolTermStructure::ParametricVolTermStructure(
            const Date& referenceDate,
            const std::vector<Date>& dates,
            const std::vector<ParametricVolSlice>& slices,
            ext::shared_ptr<ParametricVolShape> shape,
            Handle<Quote> spot,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            DividendSchedule dividends,
            const DayCounter& dc)
        : BlackVolatilityTermStructure(referenceDate, Calendar(), Following, dc),
          T_(datesToTimes(referenceDate, dates, dc)),
          slices_(slices),
          shape_(std::move(shape)),
          spot_(std::move(spot)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          dividends_(std::move(dividends))
    {
        validateSlices(T_, slices_);
        QL_REQUIRE(shape_, "shape must not be null");
        QL_REQUIRE(!spot_.empty(), "spot handle must not be empty");
        QL_REQUIRE(!riskFreeRate_.empty(), "risk-free rate handle must not be empty");
        QL_REQUIRE(!dividendYield_.empty(), "dividend yield handle must not be empty");
        registerWith(spot_);
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
    }

    Real ParametricVolTermStructure::forward(Time t) const {
        Real S = spot_->value();
        Real df = riskFreeRate_->discount(t);
        Real dq = dividendYield_->discount(t);

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

    Real ParametricVolTermStructure::z(Size i, Real k) const {
        const auto& s = slices_.at(i);
        Real sigmaHat = s.atmIv * std::sqrt(T_[i]);
        return k / sigmaHat;
    }

    Real ParametricVolTermStructure::totalVariance(Size i, Real k) const {
        const auto& s = slices_.at(i);
        Real sigmaHat = s.atmIv * std::sqrt(T_[i]);
        Real zv = k / sigmaHat;
        Real fv = shape_->f(zv, s.params);
        // w = T * atm_iv^2 * f(z)  =  sigmaHat^2 * f(z)
        return sigmaHat * sigmaHat * fv;
    }

    Real ParametricVolTermStructure::totalVarianceStrikeDerivative(
            Size i, Real k) const {
        const auto& s = slices_.at(i);
        Real sigmaHat = s.atmIv * std::sqrt(T_[i]);
        Real zv = k / sigmaHat;
        Real dfdz = shape_->dfdz(zv, s.params);
        // dw/dk = sigmaHat^2 * df/dk = sigmaHat^2 * df/dz * dz/dk
        //       = sigmaHat^2 * df/dz / sigmaHat = sigmaHat * df/dz
        return sigmaHat * dfdz;
    }

    ParametricVolTermStructure::TimeWeights
    ParametricVolTermStructure::timeWeights(Time t) const {
        const Size N = T_.size();
        TimeWeights tw{};

        if (t <= T_.front()) {
            tw.sLo = 0;
            tw.sHi = 0;
            tw.wLo = t / T_.front();
            tw.wHi = 0.0;
            tw.leftExtrap = true;
            return tw;
        }
        tw.leftExtrap = false;

        if (N == 1 || t >= T_.back()) {
            if (N == 1) {
                tw.sLo = 0;
                tw.sHi = 0;
                tw.wLo = 1.0;
                tw.wHi = 0.0;
                return tw;
            }
            tw.sLo = N - 2;
            tw.sHi = N - 1;
            Time dT = T_[N-1] - T_[N-2];
            Real a = (t - T_[N-2]) / dT;
            tw.wLo = 1.0 - a;
            tw.wHi = a;
            return tw;
        }

        auto it = std::upper_bound(T_.begin(), T_.end(), t);
        Size hi = static_cast<Size>(it - T_.begin());
        Size lo = hi - 1;
        Time dT = T_[hi] - T_[lo];
        Real a = (t - T_[lo]) / dT;
        tw.sLo = lo;
        tw.sHi = hi;
        tw.wLo = 1.0 - a;
        tw.wHi = a;
        return tw;
    }

    Real ParametricVolTermStructure::totalVariance(Real k, Time t) const {
        QL_REQUIRE(t > 0.0, "time must be > 0");
        TimeWeights tw = timeWeights(t);
        Real w_lo = totalVariance(tw.sLo, k);
        if (tw.leftExtrap || tw.sHi == tw.sLo)
            return tw.wLo * w_lo;
        Real w_hi = totalVariance(tw.sHi, k);
        return tw.wLo * w_lo + tw.wHi * w_hi;
    }

    Real ParametricVolTermStructure::totalVarianceStrikeDerivative(
            Real k, Time t) const {
        QL_REQUIRE(t > 0.0, "time must be > 0");
        TimeWeights tw = timeWeights(t);
        Real d_lo = totalVarianceStrikeDerivative(tw.sLo, k);
        if (tw.leftExtrap || tw.sHi == tw.sLo)
            return tw.wLo * d_lo;
        Real d_hi = totalVarianceStrikeDerivative(tw.sHi, k);
        return tw.wLo * d_lo + tw.wHi * d_hi;
    }

    Real ParametricVolTermStructure::totalVarianceStrikeSecondDerivative(
            Size i, Real k) const {
        // w(T,k) = sigmaHat² · f(z),  z = k/sigmaHat,  dz/dk = 1/sigmaHat
        // ⇒ d²w/dk² = sigmaHat² · f''(z) · (1/sigmaHat)² = f''(z)
        const auto& s = slices_.at(i);
        Real sigmaHat = s.atmIv * std::sqrt(T_[i]);
        Real zv = k / sigmaHat;
        return shape_->d2fdz2(zv, s.params);
    }

    Real ParametricVolTermStructure::totalVarianceStrikeSecondDerivative(
            Real k, Time t) const {
        QL_REQUIRE(t > 0.0, "time must be > 0");
        TimeWeights tw = timeWeights(t);
        Real d_lo = totalVarianceStrikeSecondDerivative(tw.sLo, k);
        if (tw.leftExtrap || tw.sHi == tw.sLo)
            return tw.wLo * d_lo;
        Real d_hi = totalVarianceStrikeSecondDerivative(tw.sHi, k);
        return tw.wLo * d_lo + tw.wHi * d_hi;
    }

    Real ParametricVolTermStructure::totalVarianceTimeDerivative(
            Real k, Time t) const {
        QL_REQUIRE(t > 0.0, "time must be > 0");
        const Size N = T_.size();
        QL_REQUIRE(N >= 2,
                   "totalVarianceTimeDerivative: needs >=2 pillars "
                   "(local vol is degenerate on a single-pillar surface)");

        // Linear-in-w time interpolation gives a piecewise-constant slope
        // per pillar interval.  We use the right-continuous convention:
        // for t exactly on a pillar, return the slope of the interval
        // starting at that pillar.
        Size lo, hi;
        if (t <= T_.front()) {
            // Left extrap: w(t,k) = (t/T_0)·w_0(k)  ⇒  dw/dT = w_0(k)/T_0
            return totalVariance(Size(0), k) / T_.front();
        } else if (t >= T_.back()) {
            // Right extrap: hold slope of last interior interval.
            lo = N - 2;
            hi = N - 1;
        } else {
            auto it = std::upper_bound(T_.begin(), T_.end(), t);
            hi = static_cast<Size>(it - T_.begin());
            lo = hi - 1;
        }
        Real w_lo = totalVariance(lo, k);
        Real w_hi = totalVariance(hi, k);
        return (w_hi - w_lo) / (T_[hi] - T_[lo]);
    }

    Real ParametricVolTermStructure::localVariance(Real k, Time t) const {
        QL_REQUIRE(t > 0.0, "localVariance: time must be > 0");

        Real w     = totalVariance(k, t);
        Real dwdk  = totalVarianceStrikeDerivative(k, t);
        Real d2wdk = totalVarianceStrikeSecondDerivative(k, t);
        Real dwdt  = totalVarianceTimeDerivative(k, t);

        // Flat-vol degenerate case: only the time derivative survives.
        if (dwdk == 0.0 && d2wdk == 0.0)
            return dwdt;

        QL_REQUIRE(w > 0.0,
                   "localVariance: non-positive total variance " << w
                   << " at k=" << k << ", t=" << t);

        // Gatheral / Dupire denominator in y = log(K/F):
        //   D = 1 − k/w·dw/dk + ¼(−¼ − 1/w + k²/w²)·(dw/dk)² + ½·d²w/dk²
        Real den1 = 1.0 - k / w * dwdk;
        Real den2 = 0.25 * (-0.25 - 1.0 / w + k * k / (w * w)) * dwdk * dwdk;
        Real den3 = 0.5 * d2wdk;
        Real den  = den1 + den2 + den3;

        QL_REQUIRE(den > 0.0,
                   "ParametricVolTermStructure::localVariance: non-positive "
                   "Dupire denominator " << den << " at k=" << k
                   << ", t=" << t << " (calendar/butterfly arb in surface)");
        return dwdt / den;
    }

    Real ParametricVolTermStructure::localVol(Real k, Time t) const {
        Real lv2 = localVariance(k, t);
        QL_REQUIRE(lv2 >= 0.0,
                   "ParametricVolTermStructure::localVol: negative local "
                   "variance " << lv2 << " at k=" << k << ", t=" << t);
        return std::sqrt(lv2);
    }

    Volatility ParametricVolTermStructure::blackVolImpl(Time t,
                                                        Real strike) const {
        if (t < 1e-14) t = 1e-14;
        Real F = forward(t);
        Real k = std::log(strike / F);
        Real w = totalVariance(k, t);
        QL_REQUIRE(w >= 0.0, "negative total variance at (k=" << k
                              << ", t=" << t << ")");
        return std::sqrt(w / t);
    }

    void ParametricVolTermStructure::setSlice(Size i,
                                              Real atmIv,
                                              const std::vector<Real>& params) {
        QL_REQUIRE(i < slices_.size(), "slice index out of range");
        QL_REQUIRE(atmIv > 0.0, "atmIv must be > 0");
        slices_[i].atmIv = atmIv;
        slices_[i].params = params;
        notifyObservers();
    }

    void ParametricVolTermStructure::setSlices(
            const std::vector<ParametricVolSlice>& slices) {
        QL_REQUIRE(slices.size() == T_.size(),
                   "slices size must match pillar count");
        for (Size i = 0; i < slices.size(); ++i)
            QL_REQUIRE(slices[i].atmIv > 0.0,
                       "slice[" << i << "] atmIv must be > 0");
        slices_ = slices;
        notifyObservers();
    }

    // ----- Batch evaluation ---------------------------------------------------

    std::vector<Real> ParametricVolTermStructure::batchBlackVol(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes) const {
        QL_REQUIRE(sliceIndices.size() == strikes.size(),
                   "sliceIndices size (" << sliceIndices.size()
                   << ") must match strikes size (" << strikes.size() << ")");
        std::vector<Real> out(strikes.size());
        for (Size j = 0; j < strikes.size(); ++j) {
            Size i = sliceIndices[j];
            QL_REQUIRE(i < slices_.size(),
                       "sliceIndices[" << j << "] = " << i
                       << " out of range (N=" << slices_.size() << ")");
            const auto& s = slices_[i];
            Real T = T_[i];
            Real F = forward(T);
            Real k = std::log(strikes[j] / F);
            Real sigmaHat = s.atmIv * std::sqrt(T);
            Real zv = k / sigmaHat;
            Real w = sigmaHat * sigmaHat * shape_->f(zv, s.params);
            QL_REQUIRE(w >= 0.0, "negative total variance at slice " << i);
            out[j] = std::sqrt(w / T);
        }
        return out;
    }

    std::vector<Real> ParametricVolTermStructure::batchBlackVolAtTimes(
            const std::vector<Real>& times,
            const std::vector<Real>& strikes) const {
        QL_REQUIRE(times.size() == strikes.size(),
                   "times size (" << times.size()
                   << ") must match strikes size (" << strikes.size() << ")");
        std::vector<Real> out(times.size());
        for (Size j = 0; j < times.size(); ++j) {
            Time t = times[j];
            if (t < 1e-14) t = 1e-14;
            Real F = forward(t);
            Real k = std::log(strikes[j] / F);
            Real w = totalVariance(k, t);
            QL_REQUIRE(w >= 0.0, "negative total variance at t=" << t);
            out[j] = std::sqrt(w / t);
        }
        return out;
    }

    std::vector<Real> ParametricVolTermStructure::batchImpliedVolGradient(
            const std::vector<Size>& sliceIndices,
            const std::vector<Real>& strikes) const {
        QL_REQUIRE(sliceIndices.size() == strikes.size(),
                   "sliceIndices size (" << sliceIndices.size()
                   << ") must match strikes size (" << strikes.size() << ")");
        QL_REQUIRE(!slices_.empty(), "no slices");
        const Size n_shape = slices_.front().params.size();
        const Size stride = 1 + n_shape;
        std::vector<Real> out(strikes.size() * stride);
        for (Size j = 0; j < strikes.size(); ++j) {
            Size i = sliceIndices[j];
            QL_REQUIRE(i < slices_.size(),
                       "sliceIndices[" << j << "] = " << i
                       << " out of range");
            const auto& s = slices_[i];
            QL_REQUIRE(s.params.size() == n_shape,
                       "slice " << i << " has " << s.params.size()
                       << " params, expected " << n_shape);
            Real T = T_[i];
            Real F = forward(T);
            Real k = std::log(strikes[j] / F);
            Real atm = s.atmIv;
            Real sigmaHat = atm * std::sqrt(T);
            Real zv = k / sigmaHat;
            Real fv = shape_->f(zv, s.params);
            QL_REQUIRE(fv > 0.0, "f(z) <= 0 at slice " << i << ", z=" << zv);
            Real sqrtF = std::sqrt(fv);
            Real dfdz1 = shape_->dfdz(zv, s.params);
            // ∂σ/∂atm = (2f - z·f') / (2·√f)
            out[j * stride + 0] = (2.0 * fv - zv * dfdz1) / (2.0 * sqrtF);
            // ∂σ/∂params_p = atm · ∂f/∂params_p / (2·√f)
            auto dfdp = shape_->dfdParams(zv, s.params);
            QL_REQUIRE(dfdp.size() == n_shape,
                       "shape returned " << dfdp.size()
                       << " grads, expected " << n_shape);
            Real scale = atm / (2.0 * sqrtF);
            for (Size p = 0; p < n_shape; ++p)
                out[j * stride + 1 + p] = scale * dfdp[p];
        }
        return out;
    }

    // ----- Arbitrage checks (Klassen 2017 §3) ---------------------------------

    Real ParametricVolTermStructure::butterflyDensity(Size i,
                                                      Real z,
                                                      Real h) const {
        const auto& s = slices_.at(i);
        Real sigmaHat = s.atmIv * std::sqrt(T_[i]);
        Real fv  = shape_->f(z, s.params);
        QL_REQUIRE(fv > 0.0, "f(z) <= 0 at z=" << z << " (pillar " << i << ")");
        Real dfdz1   = shape_->dfdz(z, s.params);
        Real d2fdz2v = shape_->d2fdz2(z, s.params, h);
        Real term1 = 1.0 - z * dfdz1 / (2.0 * fv);
        Real g = term1 * term1
                 - 0.25 * dfdz1 * dfdz1 / fv
                 - (sigmaHat * sigmaHat / 16.0) * dfdz1 * dfdz1
                 + 0.5 * d2fdz2v;
        return g;
    }

    Real ParametricVolTermStructure::butterflyArbViolation(
            Size i,
            const std::vector<Real>& zGrid,
            Real h) const {
        Real minG = std::numeric_limits<Real>::infinity();
        for (Real z : zGrid) {
            Real g = butterflyDensity(i, z, h);
            if (g < minG) minG = g;
        }
        return std::max(0.0, -minG);
    }

    Real ParametricVolTermStructure::calendarArbViolation(
            const std::vector<Real>& kGrid) const {
        Real worst = 0.0;
        for (Size i = 0; i + 1 < slices_.size(); ++i) {
            for (Real k : kGrid) {
                Real wLo = totalVariance(i, k);
                Real wHi = totalVariance(i + 1, k);
                Real def = wLo - wHi;  // > 0 means calendar arb
                if (def > worst) worst = def;
            }
        }
        return worst;
    }

    // ----- Constraint-grid evaluators (for SLSQP / trust-constr) --------------

    std::vector<Real> ParametricVolTermStructure::butterflyDensityGrid(
            const std::vector<Real>& zGrid, Real h) const {
        const Size N = slices_.size();
        const Size n_z = zGrid.size();
        std::vector<Real> out(N * n_z);
        for (Size i = 0; i < N; ++i) {
            for (Size j = 0; j < n_z; ++j) {
                out[i * n_z + j] = butterflyDensity(i, zGrid[j], h);
            }
        }
        return out;
    }

    namespace {
        // Stateless g(z) on a free-standing slice spec — used to FD-bump
        // params/atm_iv for the gradient without mutating the stored
        // surface slices.
        inline Real gValue(const ParametricVolShape& shape,
                           Real atmIv, Real T,
                           const std::vector<Real>& params,
                           Real z, Real h) {
            Real sigmaHat = atmIv * std::sqrt(T);
            Real fv = shape.f(z, params);
            QL_REQUIRE(fv > 0.0, "f(z) <= 0 in butterfly gradient");
            Real dfdz1 = shape.dfdz(z, params);
            Real d2fdz2v = shape.d2fdz2(z, params, h);
            Real term1 = 1.0 - z * dfdz1 / (2.0 * fv);
            return term1 * term1
                   - 0.25 * dfdz1 * dfdz1 / fv
                   - (sigmaHat * sigmaHat / 16.0) * dfdz1 * dfdz1
                   + 0.5 * d2fdz2v;
        }
    }

    std::vector<Real>
    ParametricVolTermStructure::butterflyDensityGridGradient(
            const std::vector<Real>& zGrid, Real h) const {
        QL_REQUIRE(!slices_.empty(), "no slices");
        const Size N = slices_.size();
        const Size n_z = zGrid.size();
        const Size n_shape = slices_.front().params.size();
        const Size stride = 1 + n_shape;
        std::vector<Real> out(N * n_z * stride);

        for (Size i = 0; i < N; ++i) {
            const auto& s = slices_[i];
            Real T = T_[i];
            std::vector<Real> params_buf = s.params;  // local mutable copy

            for (Size j = 0; j < n_z; ++j) {
                Real z = zGrid[j];

                // ∂g/∂atm_iv via central FD.
                Real h_a = std::max(s.atmIv * 1e-5, 1e-7);
                Real gp_a = gValue(*shape_, s.atmIv + h_a, T, params_buf, z, h);
                Real gm_a = gValue(*shape_, s.atmIv - h_a, T, params_buf, z, h);
                out[(i * n_z + j) * stride + 0] = (gp_a - gm_a) / (2.0 * h_a);

                // ∂g/∂params[p] via central FD on the copied params buffer.
                for (Size p = 0; p < n_shape; ++p) {
                    Real saved = params_buf[p];
                    Real h_p = std::max(std::abs(saved) * 1e-5, 1e-7);
                    params_buf[p] = saved + h_p;
                    Real gp_p = gValue(*shape_, s.atmIv, T, params_buf, z, h);
                    params_buf[p] = saved - h_p;
                    Real gm_p = gValue(*shape_, s.atmIv, T, params_buf, z, h);
                    params_buf[p] = saved;
                    out[(i * n_z + j) * stride + 1 + p] =
                        (gp_p - gm_p) / (2.0 * h_p);
                }
            }
        }
        return out;
    }

    std::vector<Real> ParametricVolTermStructure::calendarDeficitGrid(
            const std::vector<Real>& kGrid) const {
        const Size N = slices_.size();
        if (N < 2) return {};
        const Size n_k = kGrid.size();
        std::vector<Real> out((N - 1) * n_k);
        for (Size i = 1; i < N; ++i) {
            for (Size j = 0; j < n_k; ++j) {
                Real w_lo = totalVariance(i - 1, kGrid[j]);
                Real w_hi = totalVariance(i, kGrid[j]);
                // Deficit ≥ 0  iff w_i ≥ w_{i-1}  iff calendar-arb-free.
                out[(i - 1) * n_k + j] = w_hi - w_lo;
            }
        }
        return out;
    }

    std::vector<Real>
    ParametricVolTermStructure::calendarDeficitGridGradient(
            const std::vector<Real>& kGrid) const {
        QL_REQUIRE(!slices_.empty(), "no slices");
        const Size N = slices_.size();
        if (N < 2) return {};
        const Size n_k = kGrid.size();
        const Size n_shape = slices_.front().params.size();
        const Size grad_per_slice = 1 + n_shape;
        // Row layout for each (pair_idx, k):
        //   [grad_w_lower (atm + params), grad_w_upper (atm + params)]
        const Size row = 2 * grad_per_slice;
        std::vector<Real> out((N - 1) * n_k * row);

        // ∂w/∂atm_iv = 2·atm_iv·T·f(z)  +  atm_iv²·T·f'(z)·∂z/∂atm_iv
        // where z = k/(atm_iv·√T), ∂z/∂atm_iv = -z/atm_iv.
        // Combined: ∂w/∂atm_iv = 2·atm_iv·T·f − atm_iv·T·z·f'·1
        //                       = atm_iv·T·(2·f − z·f')
        // ∂w/∂params_p = atm_iv²·T·∂f/∂params_p
        for (Size i = 1; i < N; ++i) {
            for (Size j = 0; j < n_k; ++j) {
                Real k = kGrid[j];
                Size base = ((i - 1) * n_k + j) * row;

                // Lower pillar contribution: gradient w.r.t. slice (i-1)
                {
                    const auto& s = slices_[i - 1];
                    Real T = T_[i - 1];
                    Real sigmaHat = s.atmIv * std::sqrt(T);
                    Real z = k / sigmaHat;
                    Real fv = shape_->f(z, s.params);
                    Real dfdz1 = shape_->dfdz(z, s.params);
                    auto dfdp = shape_->dfdParams(z, s.params);
                    Real datm = s.atmIv * T * (2.0 * fv - z * dfdz1);
                    // Sign: deficit = w_i − w_{i-1}, so wrt w_{i-1} → negate
                    out[base + 0] = -datm;
                    for (Size p = 0; p < n_shape; ++p) {
                        out[base + 1 + p] = -(s.atmIv * s.atmIv) * T * dfdp[p];
                    }
                }
                // Upper pillar contribution: gradient w.r.t. slice i
                {
                    const auto& s = slices_[i];
                    Real T = T_[i];
                    Real sigmaHat = s.atmIv * std::sqrt(T);
                    Real z = k / sigmaHat;
                    Real fv = shape_->f(z, s.params);
                    Real dfdz1 = shape_->dfdz(z, s.params);
                    auto dfdp = shape_->dfdParams(z, s.params);
                    Real datm = s.atmIv * T * (2.0 * fv - z * dfdz1);
                    out[base + grad_per_slice + 0] = datm;
                    for (Size p = 0; p < n_shape; ++p) {
                        out[base + grad_per_slice + 1 + p] =
                            (s.atmIv * s.atmIv) * T * dfdp[p];
                    }
                }
            }
        }
        return out;
    }

}  // namespace QuantLib
