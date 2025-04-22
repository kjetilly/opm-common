// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::EclTwoPhaseMaterialParams
 */
#ifndef OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <memory>

#include <opm/common/utility/gpuDecorators.hpp>
#if OPM_IS_COMPILING_WITH_GPU_COMPILER
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#endif
#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {

enum class EclTwoPhaseApproach {
    GasOil,
    OilWater,
    GasWater
};

#if OPM_IS_COMPILING_WITH_GPU_COMPILER
template<class T>
struct NoPointer {
    using type = T;
};
template<class OriginalContainer>
struct TransformContainer {
    using type = NoPointer<OriginalContainer>; //Opm::gpuistl::GpuView<typename OriginalContainer::value_type>;
};
#else
template<class OriginalContainer>
struct TransformContainer {
    using type = OriginalContainer;
};
#endif

/*!
 * \brief Implementation for the parameters required by the material law for two-phase
 *        simulations.
 *
 * Essentially, this class just stores the two parameter objects for
 * the twophase capillary pressure laws.
 */
template<class Traits, class GasOilParamsT, class OilWaterParamsT, class GasWaterParamsT>
class EclTwoPhaseMaterialParams : public EnsureFinalized
{
    using Scalar = typename Traits::Scalar;
    enum { numPhases = 3 };
public:
    using EnsureFinalized :: finalize;

    using GasOilParams = typename TransformContainer<GasOilParamsT>::type;
    using OilWaterParams = typename TransformContainer<OilWaterParamsT>::type;
    using GasWaterParams = typename TransformContainer<GasWaterParamsT>::type;

    #if OPM_IS_COMPILING_WITH_GPU_COMPILER
    template<class T>
    using SmartPointer = Opm::gpuistl::PointerView<T>;
    #else
    template<class T>
    using SmartPointer = std::shared_ptr<T>;
    #endif

    /*!
     * \brief The default constructor.
     */
    OPM_HOST_DEVICE EclTwoPhaseMaterialParams()
    {
    }

    OPM_HOST_DEVICE void setApproach(EclTwoPhaseApproach newApproach)
    { approach_ = newApproach; }

    OPM_HOST_DEVICE EclTwoPhaseApproach approach() const
    { return approach_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
     OPM_HOST_DEVICE const GasOilParams& gasOilParams() const
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
     OPM_HOST_DEVICE GasOilParams& gasOilParams()
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief Set the parameter object for the gas-oil twophase law.
     */
     OPM_HOST_DEVICE void setGasOilParams(SmartPointer<GasOilParams> val)
    { gasOilParams_ = val; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
     OPM_HOST_DEVICE const OilWaterParams& oilWaterParams() const
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
     OPM_HOST_DEVICE OilWaterParams& oilWaterParams()
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief Set the parameter object for the oil-water twophase law.
     */
    void setOilWaterParams(SmartPointer<OilWaterParams> val)
    { oilWaterParams_ = val; }

  /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    const GasWaterParams& gasWaterParams() const
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    GasWaterParams& gasWaterParams()
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief Set the parameter object for the gas-water twophase law.
     */
    void setGasWaterParams(SmartPointer<GasWaterParams> val)
    { gasWaterParams_ = val; }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        // This is for restart serialization.
        // Only dynamic state in the parameters need to be stored.
        serializer(*gasOilParams_);
        serializer(*oilWaterParams_);
        serializer(*gasWaterParams_);
    }

    void setSwl(Scalar) {}

private:
    EclTwoPhaseApproach approach_{EclTwoPhaseApproach::GasOil};

    SmartPointer<GasOilParams> gasOilParams_{nullptr};
    SmartPointer<OilWaterParams> oilWaterParams_{nullptr};
    SmartPointer<GasWaterParams> gasWaterParams_{nullptr};
};

} // namespace Opm

#endif
