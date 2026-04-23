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

#include <ql/termstructures/volatility/equityfx/pwlpdflocalvolsurface.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    PwlPdfLocalVolSurface::PwlPdfLocalVolSurface(
            ext::shared_ptr<PwlPdfVolSurface> pdfSurface,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Handle<Quote> spot)
        : LocalVolTermStructure(pdfSurface->businessDayConvention(),
                                pdfSurface->dayCounter()),
          pdfSurface_(std::move(pdfSurface)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          spot_(std::move(spot))
    {
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
        registerWith(spot_);
    }

    const Date& PwlPdfLocalVolSurface::referenceDate() const {
        return pdfSurface_->referenceDate();
    }

    DayCounter PwlPdfLocalVolSurface::dayCounter() const {
        return pdfSurface_->dayCounter();
    }

    Date PwlPdfLocalVolSurface::maxDate() const {
        return pdfSurface_->maxDate();
    }

    Real PwlPdfLocalVolSurface::minStrike() const {
        return pdfSurface_->minStrike();
    }

    Real PwlPdfLocalVolSurface::maxStrike() const {
        return pdfSurface_->maxStrike();
    }

    void PwlPdfLocalVolSurface::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<PwlPdfLocalVolSurface>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            LocalVolTermStructure::accept(v);
    }

    // =================================================================
    // Dupire local vol from PWL PDF
    //
    // Full Dupire in (K, T) space:
    //   sigma_loc^2 = [dC/dT + (r-q)*K*dC/dK + q*C] / [0.5 * K^2 * d^2C/dK^2]
    //
    // Dupire in (K, T) space for undiscounted call C:
    //   sigma_loc^2 = [dC/dT + (r-q)*K*dC/dK] / [0.5 * K^2 * p(K)]
    //
    // Denominator: p(K) = q_x(x) / F => 0.5 * K^2 * q_x / F
    // Numerator: dC/dT via FD of analytically-computed call prices
    //            + drift correction via survival function
    // dC/dK = -S_x(x) where S_x = integral of q_x from x to x_max
    // =================================================================

    Volatility PwlPdfLocalVolSurface::localVolImpl(Time t,
                                                    Real underlyingLevel)
                                                                      const {
        if (t < 1e-14) t = 1e-14;

        const Real volFloor = 0.005;
        const Real volCap = 5.0;
        const Real pdfFloor = 1e-12;

        // Forward and moneyness
        DiscountFactor dr = riskFreeRate_->discount(t, true);
        DiscountFactor dq = dividendYield_->discount(t, true);
        Real S = spot_->value();
        Real F = S * dq / dr;
        Real K = underlyingLevel;
        Real x = K / F;  // moneyness

        // Denominator: 0.5 * K^2 * p(K)
        // where p(K) = q_x(x) / F (strike-space PDF from moneyness PDF)
        // => 0.5 * K^2 * q_x / F
        Real q_val = pdfSurface_->pdfValue(x, t);
        if (q_val < pdfFloor) {
            // Wing: PDF too small for Dupire.  Fall back to black vol
            // which uses flat wing extrapolation at boundary vols.
            // Flat black vol → local vol = black vol, so this is
            // self-consistent and keeps local vol finite.
            return std::max(pdfSurface_->blackVol(t, K, true), volFloor);
        }

        Real denominator = 0.5 * K * K * q_val / F;

        // Numerator: dC/dT at fixed K
        // C(K,T) = DF(T) * F(T) * c(k(K,T), T)
        //
        // Use a small time bump for the numerator.
        // Since the call-forward interpolation is piecewise-linear in T,
        // the derivative is piecewise-constant. We evaluate the left
        // and right limits to get the appropriate slope at T.
        const auto& slices = pdfSurface_->pdfSlices();
        Size N = slices.size();

        Real dCdT;

        if (N == 1) {
            // Single slice: no term structure info — use black vol
            return std::max(pdfSurface_->blackVol(t, K, true), volFloor);
        } else {
            // Find bracketing slices
            Size lo, hi;
            if (t <= slices.front().T) {
                // Before first slice: use first two slices for FD
                lo = 0; hi = 1;
            } else if (t >= slices.back().T) {
                // After last slice: use last two slices
                lo = N - 2; hi = N - 1;
            } else {
                hi = 0;
                for (Size i = 1; i < N; ++i) {
                    if (slices[i].T >= t) { hi = i; break; }
                }
                lo = hi - 1;
            }

            // Use calibration forwards (same as blackVolImpl)
            Real T_lo = slices[lo].T;
            Real T_hi = slices[hi].T;
            Real F_lo = slices[lo].calForward > 0.0
                ? slices[lo].calForward
                : S * dividendYield_->discount(T_lo, true)
                    / riskFreeRate_->discount(T_lo, true);
            Real F_hi = slices[hi].calForward > 0.0
                ? slices[hi].calForward
                : S * dividendYield_->discount(T_hi, true)
                    / riskFreeRate_->discount(T_hi, true);

            Real x_lo = K / F_lo;
            Real x_hi = K / F_hi;

            Real C_lo = F_lo * pdfSurface_->callForward(x_lo, T_lo);
            Real C_hi = F_hi * pdfSurface_->callForward(x_hi, T_hi);

            // Calendar spread = dC/dT at fixed K
            dCdT = (C_hi - C_lo) / (T_hi - T_lo);
            if (dCdT < 0.0) dCdT = 0.0;
        }

        // Full Dupire for undiscounted call:
        //   sigma^2 = [dC/dT + (r-q)*K*dC/dK] / [0.5 * K^2 * p(K)]
        //
        // dC/dK = -S(k) where S(k) = survival function = integral q_k from k
        // The (r-q)*K*(-S(k)) drift term is significant for equities.
        //
        // Instantaneous rates from yield curves:
        Real r_inst = riskFreeRate_->forwardRate(t, t + 1e-4,
                                                 Continuous, NoFrequency).rate();
        Real q_inst = dividendYield_->forwardRate(t, t + 1e-4,
                                                  Continuous, NoFrequency).rate();
        Real drift = r_inst - q_inst;

        // Survival function: dC/dK = -S(x)
        Real survK = pdfSurface_->survivalFunction(x, t);
        Real dCdK = -survK;

        Real numerator = dCdT + drift * K * dCdK;
        if (numerator < 0.0) numerator = 0.0;

        Real localVar = numerator / denominator;

        if (localVar <= 0.0 || !std::isfinite(localVar))
            return std::max(pdfSurface_->blackVol(t, K, true), volFloor);

        Real localVol = std::sqrt(localVar);

        // Relative cap: local vol ≤ 3× black vol at same (t, K).
        // Prevents transition-zone spikes where PDF is tiny but
        // above pdfFloor.  blackVol uses flat wing extrapolation
        // at boundary vols, so the cap is always finite and smooth.
        Real bvol = std::max(pdfSurface_->blackVol(t, K, true), volFloor);
        Real relCap = std::min(3.0 * bvol, volCap);
        return std::min(std::max(localVol, volFloor), relCap);
    }

} // namespace QuantLib
