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
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/math/interpolations/cubicinterpolation.hpp>
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
    // Gatheral-Dupire local vol via PwlPdfVolSurface::localVol(k, t)
    // =================================================================

    Volatility PwlPdfLocalVolSurface::localVolImpl(Time t,
                                                    Real underlyingLevel)
                                                                      const {
        if (t < 1e-14) t = 1e-14;
        Real fwd = pdfSurface_->forward(t);
        Real k = std::log(underlyingLevel / fwd);
        return pdfSurface_->localVol(k, t);
    }

    std::vector<Volatility> PwlPdfLocalVolSurface::localVolGrid(
            const std::vector<Time>& times,
            const std::vector<Real>& underlyingLevels) const {
        const Size nT = times.size();
        const Size nS = underlyingLevels.size();
        std::vector<Volatility> out(nT * nS);
        std::vector<Real> logS(nS);
        for (Size j = 0; j < nS; ++j)
            logS[j] = std::log(underlyingLevels[j]);
        for (Size i = 0; i < nT; ++i) {
            Time t = times[i];
            if (t < 1e-14) t = 1e-14;
            Real fwd = pdfSurface_->forward(t);
            Real logFwd = std::log(fwd);
            Volatility* row = out.data() + i * nS;
            for (Size j = 0; j < nS; ++j) {
                Real k = logS[j] - logFwd;
                row[j] = pdfSurface_->localVol(k, t);
            }
        }
        return out;
    }

} // namespace QuantLib
