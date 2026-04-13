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

/*! \file g2copulabondoptionengine.hpp
    \brief European bond option engine using G2 copula + SABR marginals

    Prices European call/put options on fixed-rate bonds by combining:
    - G2 two-factor model for swap rate correlation structure
    - SABR (or flat Bachelier) marginal distributions per co-initial tenor
    - Gaussian copula to join the marginals
    - Gauss-Hermite quadrature for deterministic evaluation

    At each quadrature point, the engine bootstraps discount factors from
    the realized swap rates and computes the bond price as sum(cashflow * DF).
    Option payoff = max(putCall * (bondPrice - strike), 0).

    Complexity: O(n_quad^2 * n_tenors) with n_quad=16 → 256 evals.
    Runtime: sub-millisecond in C++ (vs 30ms in Python).
*/

#ifndef quantlib_g2_copula_bond_option_engine_hpp
#define quantlib_g2_copula_bond_option_engine_hpp

#include <ql/models/shortrate/twofactormodels/g2.hpp>
#include <ql/instruments/bond.hpp>
#include <ql/option.hpp>
#include <ql/math/array.hpp>
#include <ql/math/matrix.hpp>
#include <vector>

namespace QuantLib {

    //! SABR parameters for a single co-initial swap rate tenor.
    struct SabrTenorParams {
        Real forward;       //!< ATM forward swap rate
        Real atmNormalVol;   //!< ATM normal (Bachelier) vol
        // Optional SABR smile params (if alpha > 0, use SABR; else flat Bachelier)
        Real alpha;
        Real beta;
        Real nu;
        Real rho;
        Real shift;

        SabrTenorParams()
        : forward(0.0), atmNormalVol(0.0),
          alpha(0.0), beta(0.0), nu(0.0), rho(0.0), shift(0.0) {}

        bool hasSabr() const { return alpha > 0.0; }
    };

    //! Single cashflow entry: (amount, yearFractionFromExpiry)
    struct CashflowEntry {
        Real amount;
        Real timeFromExpiry;

        CashflowEntry() : amount(0.0), timeFromExpiry(0.0) {}
    };

    //! European bond option via G2 copula + SABR marginals.
    /*! Standalone pricing function — not tied to the QL instrument framework.
        Takes bond cashflows, option terms, G2 model, and SABR params directly.

        \ingroup bondengines
    */

    //! Compute European bond option PV via G2 copula quadrature.
    /*!
        \param model         Calibrated G2 model
        \param cashflows     Bond cashflows after option expiry
        \param tenors        Sparse co-initial tenor grid (e.g., {0.5, 1, 2, 5, 10, 20})
        \param sabrParams    SABR params per tenor (must match tenors.size())
        \param strike        Option strike (dirty price per 100 par)
        \param expiryTime    Time to option expiry in years
        \param type          Option::Call or Option::Put
        \param quadOrder     Gauss-Hermite quadrature order per dimension (default 16)
        \param nCdfStrikes   Number of strikes for CDF grid (default 200)
        \param nCdfStd       Standard deviations for CDF strike range (default 4.0)
        \return              Option PV (NOT discounted to today — caller applies DF)
    */
    Real g2CopulaBondOptionPV(
        const ext::shared_ptr<G2>& model,
        const std::vector<CashflowEntry>& cashflows,
        const std::vector<Real>& tenors,
        const std::vector<SabrTenorParams>& sabrParams,
        Real strike,
        Time expiryTime,
        Option::Type type,
        Size quadOrder = 16,
        Size nCdfStrikes = 200,
        Real nCdfStd = 4.0);

}

#endif
