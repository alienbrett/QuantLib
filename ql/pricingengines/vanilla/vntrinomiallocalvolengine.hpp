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

        Adaptive branching: when local vol exceeds the standard trinomial
        threshold (σ_loc > √3·σ_ref), the tree branches to j±k instead
        of j±1, keeping probabilities non-negative without clamping.

        Discrete cash dividends are handled via Vellekoop-Nieuwenhuis
        interpolation, identical to VNBinomialVanillaEngine.

        Optional coarse-grid local vol: when lvGridStride > 0, the engine
        precomputes local vol on a coarser spatial grid (every stride-th
        node) and linearly interpolates during tree traversal.  This is
        accurate for smooth analytic surfaces (eSSVI) and gives 3-6×
        speedup.

        \ingroup vanillaengines
    */
    class VNTrinomialLocalVolEngine : public VanillaOption::engine {
      public:
        //! Standard: extracts local vol from process (generic Dupire)
        VNTrinomialLocalVolEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size timeSteps);

        //! Fast path: uses an explicit local vol surface (e.g. EssviLocalVolSurface)
        VNTrinomialLocalVolEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size timeSteps,
            ext::shared_ptr<LocalVolTermStructure> localVol);

        //! Fast path with coarse-grid interpolation (lvGridStride > 0)
        VNTrinomialLocalVolEngine(
            ext::shared_ptr<GeneralizedBlackScholesProcess> process,
            DividendSchedule dividends,
            Size timeSteps,
            ext::shared_ptr<LocalVolTermStructure> localVol,
            Size lvGridStride);

        void calculate() const override;
      private:
        ext::shared_ptr<GeneralizedBlackScholesProcess> process_;
        DividendSchedule dividends_;
        Size timeSteps_;
        ext::shared_ptr<LocalVolTermStructure> explicitLocalVol_;
        Size lvGridStride_ = 0;
    };

}

#endif
