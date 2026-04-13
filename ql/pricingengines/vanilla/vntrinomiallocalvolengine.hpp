/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/*! \file vntrinomiallocalvolengine.hpp
    \brief Trinomial local-vol tree with VN discrete-dividend interpolation
*/

#ifndef quantlib_vn_trinomial_localvol_engine_hpp
#define quantlib_vn_trinomial_localvol_engine_hpp

#include <ql/instruments/vanillaoption.hpp>
#include <ql/instruments/dividendschedule.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    //! Trinomial tree with Dupire local vol and VN dividend interpolation
    /*! Prices European and American vanilla options using a recombining
        trinomial tree where each node's transition probabilities are
        calibrated to the Dupire local volatility surface.

        Three probabilities (pu, pm, pd) at each node satisfy three
        constraints: sum-to-one, forward match, and local variance match.
        This exactly captures the smile dynamics at every node.

        Discrete cash dividends are handled via Vellekoop-Nieuwenhuis
        interpolation, identical to VNBinomialVanillaEngine.

        Grid spacing: dx = sigma_ref * sqrt(3 * dt), where sigma_ref is
        the ATM implied vol.  This ensures pu, pm, pd stay well within
        (0, 1) for typical local vol values.

        \ingroup vanillaengines
    */
    class VNTrinomialLocalVolEngine : public VanillaOption::engine {
      public:
        VNTrinomialLocalVolEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size timeSteps);
        void calculate() const override;
      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        Size timeSteps_;
    };

}

#endif
