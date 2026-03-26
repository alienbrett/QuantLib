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

#include <ql/termstructures/volatility/equityfx/essvilocalvolsurface.hpp>
#include <utility>

namespace QuantLib {

    EssviLocalVolSurface::EssviLocalVolSurface(
            ext::shared_ptr<EssviVolatilityTermStructure> essviSurface,
            Handle<YieldTermStructure> riskFreeRate,
            Handle<YieldTermStructure> dividendYield,
            Handle<Quote> spot)
        : LocalVolTermStructure(essviSurface->businessDayConvention(),
                                essviSurface->dayCounter()),
          essviSurface_(std::move(essviSurface)),
          riskFreeRate_(std::move(riskFreeRate)),
          dividendYield_(std::move(dividendYield)),
          spot_(std::move(spot))
    {
        registerWith(riskFreeRate_);
        registerWith(dividendYield_);
        registerWith(spot_);
    }

    const Date& EssviLocalVolSurface::referenceDate() const {
        return essviSurface_->referenceDate();
    }

    DayCounter EssviLocalVolSurface::dayCounter() const {
        return essviSurface_->dayCounter();
    }

    Date EssviLocalVolSurface::maxDate() const {
        return essviSurface_->maxDate();
    }

    Real EssviLocalVolSurface::minStrike() const {
        return essviSurface_->minStrike();
    }

    Real EssviLocalVolSurface::maxStrike() const {
        return essviSurface_->maxStrike();
    }

    void EssviLocalVolSurface::accept(AcyclicVisitor& v) {
        auto* v1 = dynamic_cast<Visitor<EssviLocalVolSurface>*>(&v);
        if (v1 != nullptr)
            v1->visit(*this);
        else
            LocalVolTermStructure::accept(v);
    }

    Volatility EssviLocalVolSurface::localVolImpl(Time t, Real underlyingLevel)
                                                                          const {
        if (t < 1e-14) t = 1e-14;

        // Forward price
        DiscountFactor dr = riskFreeRate_->discount(t, true);
        DiscountFactor dq = dividendYield_->discount(t, true);
        Real fwd = spot_->value() * dq / dr;

        // Log-moneyness
        Real k = std::log(underlyingLevel / fwd);

        // Delegate to the inner EssviSurface's analytic Dupire formula.
        // We need to use the forward from the eSSVI surface's own curves
        // for consistency. But the eSSVI surface computes its own forward
        // internally, and we need log-moneyness relative to that forward.
        // Since we constructed with the same yield curves, our forward
        // matches the eSSVI surface's forward, so k is consistent.

        // Access the inner surface's analytic local vol
        const auto& surface = essviSurface_->surface();
        return surface.localVol(k, t);
    }

} // namespace QuantLib
