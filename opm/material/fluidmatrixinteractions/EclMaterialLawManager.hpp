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
 * \copydoc Opm::EclMaterialLawManager
 */

#ifndef OPM_ECL_MATERIAL_LAW_MANAGER_HPP
#define OPM_ECL_MATERIAL_LAW_MANAGER_HPP

#include "opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp"
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/WagHysteresisConfig.hpp>

#include <opm/common/utility/VectorWithDefaultAllocator.hpp>
#include <opm/common/utility/gpuistl_if_available.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawTwoPhaseTypes.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/SatCurveMultiplexer.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/DirectionalMaterialLawParams.hpp>

#include <cassert>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace Opm {

class EclipseState;
class EclEpsGridProperties;
template<class Scalar> class EclEpsScalingPoints;
template<class Scalar> struct EclEpsScalingPointsInfo;
class EclHysteresisConfig;
enum class EclTwoPhaseSystemType;
class FieldPropsManager;
class Runspec;
class SgfnTable;
class SgofTable;
class SlgofTable;
class TableColumn;

}

namespace Opm::EclMaterialLaw {

template<class Traits> class InitParams;

template<class, class, class, class>
struct DefaultParams {};

/// Helper to detect if two template-template parameters are the same.
template <template <class, class, class, class> class A,
          template <class, class, class, class> class B>
struct IsSameTemplateTemplate : std::false_type {};
template <template <class, class, class, class> class A>
struct IsSameTemplateTemplate<A, A> : std::true_type {};

/*!
 * \ingroup fluidmatrixinteractions
 *
 * \brief Provides an simple way to create and manage the material law objects
 *        for a complete ECL deck.
 *
 * \tparam TraitsT  The material traits type (phase indices, scalar type, etc.)
 * \tparam ParamsT  A template-template parameter for the material law parameters type.
 *                  When set to DefaultParams (the default), the multiplexer params
 *                  (EclMultiplexerMaterialParams) are used. Otherwise, ParamsT is
 *                  instantiated with <Traits, GasOilLaw, OilWaterLaw, GasWaterLaw>
 *                  and used directly as the MaterialLawParams type.
 * \tparam Storage  A template for container types (default: VectorWithDefaultAllocator,
 *                  i.e. std::vector). Can be replaced with e.g. GpuBuffer for GPU storage.
 * \tparam SharedPointer  A template for shared pointer types (default: std::shared_ptr).
 */
template <class TraitsT,
          template <class, class, class, class> class ParamsT = DefaultParams,
          template<class> class Storage = VectorWithDefaultAllocator,
          template<class> class SharedPointer = std::shared_ptr>
class Manager
{
    using Traits = TraitsT;
    using Scalar = typename Traits::Scalar;
    static constexpr int gasPhaseIdx = Traits::gasPhaseIdx;
    static constexpr int oilPhaseIdx = Traits::nonWettingPhaseIdx;
    static constexpr int waterPhaseIdx = Traits::wettingPhaseIdx;
    static constexpr int numPhases = Traits::numPhases;
    using GasOilEffectiveParamVector = typename EclMaterialLaw::TwoPhaseTypes<Traits>::GasOilEffectiveParamVector;
    using GasWaterEffectiveParamVector = typename EclMaterialLaw::TwoPhaseTypes<Traits>::GasWaterEffectiveParamVector;
    using OilWaterEffectiveParamVector = typename EclMaterialLaw::TwoPhaseTypes<Traits>::OilWaterEffectiveParamVector;

    using GasOilLaw = typename EclMaterialLaw::TwoPhaseTypes<Traits>::GasOilLaw;
    using OilWaterLaw = typename EclMaterialLaw::TwoPhaseTypes<Traits>::OilWaterLaw;
    using GasWaterLaw = typename EclMaterialLaw::TwoPhaseTypes<Traits>::GasWaterLaw;

public:
    // the three-phase material law used by the simulation
    using MaterialLaw = EclMultiplexerMaterial<Traits, GasOilLaw, OilWaterLaw, GasWaterLaw>;

    /*!
     * \brief The material law parameters type.
     *
     * When ParamsT is not DefaultParams, ParamsT<Traits, GasOilLaw, OilWaterLaw, GasWaterLaw>
     * is used directly. Otherwise, the MaterialLaw's default Params type
     * (EclMultiplexerMaterialParams) is used.
     */
    using MaterialLawParams = std::conditional_t<
        IsSameTemplateTemplate<ParamsT, DefaultParams>::value,
        typename MaterialLaw::Params,
        ParamsT<Traits, GasOilLaw, OilWaterLaw, GasWaterLaw>>;

    using DirectionalMaterialLawParamsPtr = std::unique_ptr<DirectionalMaterialLawParams<MaterialLawParams>>;

private:
    using GasOilScalingPointsVector = Storage<SharedPointer<EclEpsScalingPoints<Scalar>>>;
    using OilWaterScalingPointsVector = Storage<SharedPointer<EclEpsScalingPoints<Scalar>>>;
    using GasWaterScalingPointsVector = Storage<SharedPointer<EclEpsScalingPoints<Scalar>>>;
    using OilWaterScalingInfoVector = Storage<EclEpsScalingPointsInfo<Scalar>>;
    using MaterialLawParamsVector = Storage<SharedPointer<MaterialLawParams>>;

public:
    struct Params
    {
        OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage{};
        GasOilEffectiveParamVector gasOilEffectiveParamVector{};
        OilWaterEffectiveParamVector oilWaterEffectiveParamVector{};
        GasWaterEffectiveParamVector gasWaterEffectiveParamVector{};
        GasOilScalingPointsVector gasOilUnscaledPointsVector{};
        OilWaterScalingPointsVector oilWaterUnscaledPointsVector{};
        GasWaterScalingPointsVector gasWaterUnscaledPointsVector{};
        Storage<int> krnumXArray{};
        Storage<int> krnumYArray{};
        Storage<int> krnumZArray{};
        Storage<int> imbnumXArray{};
        Storage<int> imbnumYArray{};
        Storage<int> imbnumZArray{};
        Storage<int> satnumRegionArray{};
        Storage<int> imbnumRegionArray{};
        Storage<MaterialLawParams> materialLawParams{};
        DirectionalMaterialLawParamsPtr dirMaterialLawParams{};
        bool onlyPiecewiseLinear = true;

        bool hasDirectionalRelperms() const
        {
            return !krnumXArray.empty() ||
                   !krnumYArray.empty() ||
                   !krnumZArray.empty();
        }

        bool hasDirectionalImbnum() const
        {
            return !imbnumXArray.empty() ||
                   !imbnumYArray.empty() ||
                   !imbnumZArray.empty();
        }
    };

    void initFromState(const EclipseState& eclState);

    // \brief Function argument 'fieldPropIntOnLeadAssigner' needed to lookup
    //        field properties of cells on the leaf grid view for CpGrid with local grid refinement.
    //        Function argument 'lookupIdxOnLevelZeroAssigner' is added to lookup, for each
    //        leaf gridview cell with index 'elemIdx', its 'lookupIdx' (index of the parent/equivalent cell on level zero).
    void initParamsForElements(const EclipseState& eclState, size_t numCompressedElems,
                               const std::function<std::vector<int>(const FieldPropsManager&, const std::string&, bool)>&
                               fieldPropIntOnLeafAssigner,
                               const std::function<unsigned(unsigned)>& lookupIdxOnLevelZeroAssigner);

    /*!
     * \brief Modify the initial condition according to the SWATINIT keyword.
     *
     * The method returns the water saturation which yields a givenn capillary
     * pressure. The reason this method is not folded directly into initFromState() is
     * that the capillary pressure given depends on the particuars of how the simulator
     * calculates its initial condition.
     */
    std::pair<Scalar, bool>
    applySwatinit(unsigned elemIdx,
                  Scalar pcow,
                  Scalar Sw);

    /// Apply SWATINIT-like scaling of oil/water capillary pressure curve at
    /// simulation restart.
    ///
    /// \param[in] elemIdx Active cell index
    ///
    /// \param[in] maxPcow Scaled maximum oil/water capillary pressure.
    ///   Typically the PPCW restart file array's entry for the
    ///   corresponding cell.
    void applyRestartSwatInit(const unsigned elemIdx, const Scalar maxPcow);

    bool enableEndPointScaling() const
    { return enableEndPointScaling_; }

    bool enablePpcwmax() const
    { return enablePpcwmax_; }

    const EclHysteresisConfig& hysteresisConfig() const
    { return hysteresisConfig_; }

    bool enableHysteresis() const
    { return hysteresisConfig_.enableHysteresis(); }

    bool enablePCHysteresis() const
    { return hysteresisConfig_.enablePCHysteresis(); }

    bool enableWettingHysteresis() const
    { return hysteresisConfig_.enableWettingHysteresis(); }

    bool enableNonWettingHysteresis() const
    { return hysteresisConfig_.enableNonWettingHysteresis(); }

    bool hasGas() const
    { return hasGas_; }

    bool hasOil() const
    { return hasOil_; }

    bool hasWater() const
    { return hasWater_; }

    const EclEpsScalingPointsInfo<Scalar>& unscaledEpsInfo(unsigned satRegionIdx) const
    { return unscaledEpsInfo_[satRegionIdx]; }

    SharedPointer<WagHysteresisConfig::WagHysteresisConfigRecord>
    wagHystersisConfig(unsigned satRegionIdx) const
    { return wagHystersisConfig_[satRegionIdx]; }

    const EclEpsConfig& gasOilConfig() const
    { return gasOilConfig_; }

    const EclEpsConfig& gasWaterConfig() const
    { return gasWaterConfig_; }

    const EclEpsConfig& oilWaterConfig() const
    { return oilWaterConfig_; }

    MaterialLawParams& materialLawParams(unsigned elemIdx)
    {
        assert(elemIdx <  params_.materialLawParams.size());
        return params_.materialLawParams[elemIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
        assert(elemIdx <  params_.materialLawParams.size());
        return params_.materialLawParams[elemIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx, FaceDir::DirEnum facedir) const
    { return materialLawParamsFunc_(elemIdx, facedir); }

    MaterialLawParams& materialLawParams(unsigned elemIdx, FaceDir::DirEnum facedir)
    { return const_cast<MaterialLawParams&>(materialLawParamsFunc_(elemIdx, facedir)); }

    /*!
     * \brief Returns a material parameter object for a given element and saturation region.
     *
     * This method changes the saturation table idx in the original material law parameter object.
     * In the context of ECL reservoir simulators, this is required to properly handle
     * wells with its own saturation table idx. In order to reset the saturation table idx
     * in the materialLawparams_ call the method with the cells satRegionIdx
     */
    const MaterialLawParams& connectionMaterialLawParams(unsigned satRegionIdx, unsigned elemIdx) const;

    int satnumRegionIdx(unsigned elemIdx) const
    { return params_.satnumRegionArray[elemIdx]; }

    int getKrnumSatIdx(unsigned elemIdx, FaceDir::DirEnum facedir) const;

    bool hasDirectionalRelperms() const
    { return params_.hasDirectionalRelperms(); }

    bool hasDirectionalImbnum() const
    { return params_.hasDirectionalImbnum(); }

    int imbnumRegionIdx(unsigned elemIdx) const
    { return params_.imbnumRegionArray[elemIdx]; }

    EclMultiplexerApproach threePhaseApproach() const
    { return threePhaseApproach_; }

    EclTwoPhaseApproach twoPhaseApproach() const
    { return twoPhaseApproach_; }

    const Storage<Scalar>& stoneEtas() const
    { return stoneEtas_; }

    template <class FluidState>
    bool updateHysteresis(const FluidState& fluidState, unsigned elemIdx)
    {
        OPM_TIMEFUNCTION_LOCAL(Subsystem::SatProps);
        if (!enableHysteresis())
            return false;
        bool changed = MaterialLaw::updateHysteresis(materialLawParams(elemIdx), fluidState);
        if (hasDirectionalRelperms() || hasDirectionalImbnum()) {
            using Dir = FaceDir::DirEnum;
            constexpr int ndim = 3;
            const Dir facedirs[] = {Dir::XPlus, Dir::YPlus, Dir::ZPlus};
            for (int i = 0; i<ndim; i++) {
                const bool ischanged =
                    MaterialLaw::updateHysteresis(materialLawParams(elemIdx, facedirs[i]), fluidState);
                changed = changed || ischanged;
            }
        }
        return changed;
    }

    void oilWaterHysteresisParams(Scalar& soMax,
                                  Scalar& swMax,
                                  Scalar& swMin,
                                  unsigned elemIdx) const;

    void setOilWaterHysteresisParams(const Scalar& soMax,
                                     const Scalar& swMax,
                                     const Scalar& swMin,
                                     unsigned elemIdx);

    void gasOilHysteresisParams(Scalar& sgmax,
                                Scalar& shmax,
                                Scalar& somin,
                                unsigned elemIdx) const;

    void setGasOilHysteresisParams(const Scalar& sgmax,
                                   const Scalar& shmax,
                                   const Scalar& somin,
                                   unsigned elemIdx);

    EclEpsScalingPoints<Scalar>& oilWaterScaledEpsPointsDrainage(unsigned elemIdx);

    const EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(size_t elemIdx) const
    { return params_.oilWaterScaledEpsInfoDrainage[elemIdx]; }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        // This is for restart serialization.
        // Only dynamic state in the parameters need to be stored.
        // For that reason we do not serialize the vector
        // as that would recreate the objects inside.
        for (auto& mat : params_.materialLawParams) {
            serializer(mat);
        }
    }

    bool satCurveIsAllPiecewiseLinear() const
    {
        return this->params_.onlyPiecewiseLinear;
    }

    // Getters for GPU copy support
    const Storage<MaterialLawParams>& getMaterialLawParams() const
    { return params_.materialLawParams; }

    EclMultiplexerApproach getThreePhaseApproach() const
    { return threePhaseApproach_; }

    EclTwoPhaseApproach getTwoPhaseApproach() const
    { return twoPhaseApproach_; }

    /*!
     * \brief Construct a Manager from pre-built GPU data.
     *
     * This constructor is intended for creating GPU (GpuBuffer/GpuView) variants
     * of the Manager. Init-phase data (scaling points, effective params, etc.)
     * is not copied -- only the runtime-essential materialLawParams are stored.
     */
    Manager(Storage<MaterialLawParams>&& materialLawParams,
            EclMultiplexerApproach threePhaseApproach,
            EclTwoPhaseApproach twoPhaseApproach,
            bool hasGas, bool hasOil, bool hasWater)
        : threePhaseApproach_(threePhaseApproach)
        , twoPhaseApproach_(twoPhaseApproach)
        , hasGas_(hasGas)
        , hasOil_(hasOil)
        , hasWater_(hasWater)
    {
        params_.materialLawParams = std::move(materialLawParams);
    }

    Manager() = default;

#if HAVE_CUDA
    // Forward declare friend functions for GPU support
    template <class T,
              template<class, class, class, class> class P,
              template<class> class S,
              template<class> class SP>
    friend Manager<T, P, gpuistl::GpuBuffer, SP>
    gpuistl::copy_to_gpu(const Manager<T, P, S, SP>&);

    template <class T,
              template<class, class, class, class> class P,
              template<class> class SP>
    friend Manager<T, P, gpuistl::GpuView, SP>
    gpuistl::make_view(Manager<T, P, gpuistl::GpuBuffer, SP>&);
#endif // HAVE_CUDA

private:
    const MaterialLawParams& materialLawParamsFunc_(unsigned elemIdx, FaceDir::DirEnum facedir) const;

    void readGlobalEpsOptions_(const EclipseState& eclState);

    void readGlobalHysteresisOptions_(const EclipseState& state);

    void readGlobalThreePhaseOptions_(const Runspec& runspec);

    bool enableEndPointScaling_{false};
    EclHysteresisConfig hysteresisConfig_;
    Storage<SharedPointer<WagHysteresisConfig::WagHysteresisConfigRecord>> wagHystersisConfig_;

    Storage<EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;

    Params params_;

    EclMultiplexerApproach threePhaseApproach_ = EclMultiplexerApproach::Default;
    // this attribute only makes sense for twophase simulations!
    EclTwoPhaseApproach twoPhaseApproach_ = EclTwoPhaseApproach::GasOil;

    Storage<Scalar> stoneEtas_;

    bool enablePpcwmax_{false};
    Storage<Scalar> maxAllowPc_;
    Storage<bool> modifySwl_;

    bool hasGas_{true};
    bool hasOil_{true};
    bool hasWater_{true};

    EclEpsConfig gasOilConfig_;
    EclEpsConfig oilWaterConfig_;
    EclEpsConfig gasWaterConfig_;
};

} // namespace Opm::EclMaterialLaw

#if HAVE_CUDA
#include <opm/common/utility/gpuistl_if_available.hpp>

namespace Opm::gpuistl {

/*!
 * \brief Copy a CPU-based EclMaterialLaw::Manager to a GpuBuffer-based variant.
 *
 * Only runtime-essential data (materialLawParams, phase approach, phase flags)
 * is copied to GPU memory. Init-phase data (scaling points, effective params,
 * hysteresis config) is not transferred.
 */
template <class TraitsT,
          template<class, class, class, class> class ParamsT,
          template<class> class Storage,
          template<class> class SharedPointer>
EclMaterialLaw::Manager<TraitsT, ParamsT, GpuBuffer, SharedPointer>
copy_to_gpu(const EclMaterialLaw::Manager<TraitsT, ParamsT, Storage, SharedPointer>& cpuManager)
{
    using MaterialLawParams = typename EclMaterialLaw::Manager<TraitsT, ParamsT, Storage, SharedPointer>::MaterialLawParams;

    return EclMaterialLaw::Manager<TraitsT, ParamsT, GpuBuffer, SharedPointer>(
        GpuBuffer<MaterialLawParams>(cpuManager.getMaterialLawParams()),
        cpuManager.getThreePhaseApproach(),
        cpuManager.getTwoPhaseApproach(),
        cpuManager.hasGas(),
        cpuManager.hasOil(),
        cpuManager.hasWater());
}

/*!
 * \brief Create a GpuView-based Manager from a GpuBuffer-based Manager.
 */
template <class TraitsT,
          template<class, class, class, class> class ParamsT,
          template<class> class SharedPointer>
EclMaterialLaw::Manager<TraitsT, ParamsT, GpuView, SharedPointer>
make_view(EclMaterialLaw::Manager<TraitsT, ParamsT, GpuBuffer, SharedPointer>& gpuManager)
{
    using MaterialLawParams = typename EclMaterialLaw::Manager<TraitsT, ParamsT, GpuBuffer, SharedPointer>::MaterialLawParams;

    auto mlpView = make_view<MaterialLawParams>(gpuManager.params_.materialLawParams);

    return EclMaterialLaw::Manager<TraitsT, ParamsT, GpuView, SharedPointer>(
        std::move(mlpView),
        gpuManager.getThreePhaseApproach(),
        gpuManager.getTwoPhaseApproach(),
        gpuManager.hasGas(),
        gpuManager.hasOil(),
        gpuManager.hasWater());
}

} // namespace Opm::gpuistl
#endif // HAVE_CUDA

#endif
