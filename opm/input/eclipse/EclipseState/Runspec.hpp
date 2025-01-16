/*
  Copyright 2016  Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify it under the terms
  of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  OPM is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_RUNSPEC_HPP
#define OPM_RUNSPEC_HPP

#include <opm/common/OpmLog/KeywordLocation.hpp>

#include <opm/input/eclipse/EclipseState/EndpointScaling.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Regdims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Tabdims.hpp>

#include <opm/input/eclipse/Schedule/Action/Actdims.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQParams.hpp>

#include <bitset>
#include <cstddef>
#include <ctime>
#include <optional>

namespace Opm {

    class Deck;

} // namespace Opm

namespace Opm {

class Phases
{
public:
    Phases() noexcept = default;
    Phases(bool oil, bool gas, bool wat,
           bool solvent = false,
           bool polymer = false,
           bool energy = false,
           bool polymw = false,
           bool foam = false,
           bool brine = false,
           bool zfraction = false) noexcept;

    static Phases serializationTestObject();

    bool active( Phase ) const noexcept;
    size_t size() const noexcept;

    bool operator==(const Phases& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(bits);
    }

private:
    std::bitset<NUM_PHASES_IN_ENUM> bits;
};

class Welldims {
public:
    Welldims() = default;
    explicit Welldims(const Deck& deck);

    static Welldims serializationTestObject();

    long long maxConnPerWell() const
    {
        return this->nCWMax;
    }

    long long maxWellsPerGroup() const
    {
        return this->nWGMax;
    }

    long long maxGroupsInField() const
    {
        return this->nGMax;
    }

    long long maxWellsInField() const
    {
        return this->nWMax;
    }

    long long maxWellListsPrWell() const
    {
        return this->nWlistPrWellMax;
    }

    long long maxDynamicWellLists() const
    {
        return this->nDynWlistMax;
    }

    const std::optional<KeywordLocation>& location() const
    {
        return this->m_location;
    }

    static bool rst_cmp(const Welldims& full_dims, const Welldims& rst_dims) {
        return full_dims.maxConnPerWell() == rst_dims.maxConnPerWell() &&
            full_dims.maxWellsPerGroup() == rst_dims.maxWellsPerGroup() &&
            full_dims.maxGroupsInField() == rst_dims.maxGroupsInField() &&
            full_dims.maxWellsInField() == rst_dims.maxWellsInField() &&
            full_dims.maxWellListsPrWell() == rst_dims.maxWellListsPrWell() &&
            full_dims.maxDynamicWellLists() == rst_dims.maxDynamicWellLists();
    }

    bool operator==(const Welldims& data) const {
        return this->location() == data.location() &&
            rst_cmp(*this, data);
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(nWMax);
        serializer(nCWMax);
        serializer(nWGMax);
        serializer(nGMax);
        serializer(nWlistPrWellMax);
        serializer(nDynWlistMax);
        serializer(m_location);
    }

private:
    long long nWMax  { 0 };
    long long nCWMax { 0 };
    long long nWGMax { 0 };
    long long nGMax  { 0 };
    long long nWlistPrWellMax  { 1 };
    long long nDynWlistMax  { 1 };
    std::optional<KeywordLocation> m_location;
};

class WellSegmentDims {
public:
    WellSegmentDims();
    explicit WellSegmentDims(const Deck& deck);

    static WellSegmentDims serializationTestObject();

    long long maxSegmentedWells() const
    {
        return this->nSegWellMax;
    }

    long long maxSegmentsPerWell() const
    {
        return this->nSegmentMax;
    }

    long long maxLateralBranchesPerWell() const
    {
        return this->nLatBranchMax;
    }

    const std::optional<KeywordLocation>& location() const
    {
        return this->location_;
    }

    bool operator==(const WellSegmentDims& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(nSegWellMax);
        serializer(nSegmentMax);
        serializer(nLatBranchMax);
        serializer(location_);
    }

private:
    long long nSegWellMax;
    long long nSegmentMax;
    long long nLatBranchMax;
    std::optional<KeywordLocation> location_;
};

class NetworkDims {
public:
    NetworkDims();
    explicit NetworkDims(const Deck& deck);

    static NetworkDims serializationTestObject();

    long long maxNONodes() const
    {
        return this->nMaxNoNodes;
    }

    long long maxNoBranches() const
    {
        return this->nMaxNoBranches;
    }

    long long maxNoBranchesConToNode() const
    {
        return this->nMaxNoBranchesConToNode;
    }

    bool extendedNetwork() const
    {
        return this->type_ == Type::Extended;
    }

    bool standardNetwork() const
    {
        return this->type_ == Type::Standard;
    }

    bool active() const
    {
        return this->extendedNetwork()
            || this->standardNetwork();
    }

    bool operator==(const NetworkDims& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(nMaxNoNodes);
        serializer(nMaxNoBranches);
        serializer(nMaxNoBranchesConToNode);
    }

private:
    enum class Type { None, Extended, Standard, };

    long long nMaxNoNodes;
    long long nMaxNoBranches;
    long long nMaxNoBranchesConToNode;
    Type type_{ Type::None };
};

class AquiferDimensions {
public:
    AquiferDimensions();
    explicit AquiferDimensions(const Deck& deck);

    static AquiferDimensions serializationTestObject();

    long long maxAnalyticAquifers() const
    {
        return this->maxNumAnalyticAquifers;
    }

    long long maxAnalyticAquiferConnections() const
    {
        return this->maxNumAnalyticAquiferConn;
    }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(this->maxNumAnalyticAquifers);
        serializer(this->maxNumAnalyticAquiferConn);
    }

private:
    long long maxNumAnalyticAquifers;
    long long maxNumAnalyticAquiferConn;
};

bool operator==(const AquiferDimensions& lhs, const AquiferDimensions& rhs);

class EclHysterConfig
{
public:
    EclHysterConfig() = default;
    explicit EclHysterConfig(const Deck& deck);

    static EclHysterConfig serializationTestObject();

    /*!
     * \brief Specify whether hysteresis is enabled or not.
     */
    //void setActive(bool yesno);

    /*!
     * \brief Returns whether hysteresis is enabled (active).
     */
    bool active() const;

    /*!
     * \brief Return the type of the hysteresis model which is used for capillary pressure.
     *
     * -1: capillary pressure hysteresis is disabled
     * 0: use the Killough model for capillary pressure hysteresis
     */
    long long pcHysteresisModel() const;

    /*!
     * \brief Return the type of the hysteresis model which is used for relative permeability.
     *
     * -1: relperm hysteresis is disabled
     * 0: use the Carlson model for relative permeability hysteresis
     */
    long long krHysteresisModel() const;

    /*!
     * \brief Regularisation parameter used for Killough model.
     *
     * default: 0.1
     */
    double modParamTrapped() const;

    /*!
     * \brief Curvature parameter used for capillary pressure hysteresis.
     *
     * default: 0.1
     */
    double curvatureCapPrs() const;

    /*!
     * \brief Wag hysteresis.
     */
    bool activeWag() const;

    bool operator==(const EclHysterConfig& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(activeHyst);
        serializer(pcHystMod);
        serializer(krHystMod);
        serializer(modParamTrappedValue);
        serializer(curvatureCapPrsValue);
        serializer(activeWagHyst);
    }

private:
    // enable hysteresis at all
    bool activeHyst  { false };

    // the capillary pressure and the relperm hysteresis models to be used
    long long pcHystMod { -1 };
    long long krHystMod { -1 };
    // regularisation parameter used for Killough model
    double modParamTrappedValue { 0.1 };
    // curvature parameter for capillary pressure
    double curvatureCapPrsValue { 0.1 };

    // enable WAG hysteresis
    bool activeWagHyst  { false };
};

class SatFuncControls {
public:
    enum class ThreePhaseOilKrModel {
        Default,
        Stone1,
        Stone2
    };

    enum class KeywordFamily {
        Family_I,               // SGOF, SWOF, SLGOF
        Family_II,              // SGFN, SOF{2,3}, SWFN, SGWFN
        Family_III,             // GSF, WSF

        Undefined,
    };

    SatFuncControls();
    explicit SatFuncControls(const Deck& deck);
    explicit SatFuncControls(const double tolcritArg,
                             const ThreePhaseOilKrModel model,
                             const KeywordFamily family);

    static SatFuncControls serializationTestObject();

    double minimumRelpermMobilityThreshold() const
    {
        return this->tolcrit;
    }

    ThreePhaseOilKrModel krModel() const
    {
        return this->krmodel;
    }

    KeywordFamily family() const
    {
        return this->satfunc_family;
    }

    bool operator==(const SatFuncControls& rhs) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(tolcrit);
        serializer(krmodel);
        serializer(satfunc_family);
    }

private:
    double tolcrit;
    ThreePhaseOilKrModel krmodel = ThreePhaseOilKrModel::Default;
    KeywordFamily satfunc_family = KeywordFamily::Undefined;
};


class Nupcol {
public:
    Nupcol();
    explicit Nupcol(long long min_value);
    void update(long long value);
    long long value() const;

    static Nupcol serializationTestObject();
    bool operator==(const Nupcol& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer) {
        serializer(this->nupcol_value);
        serializer(this->min_nupcol);
    }

private:
    long long min_nupcol;
    long long nupcol_value;
};


class Tracers {
public:
    Tracers() = default;

    explicit Tracers(const Deck& );
    long long water_tracers() const;

    template<class Serializer>
    void serializeOp(Serializer& serializer) {
        serializer(this->m_oil_tracers);
        serializer(this->m_water_tracers);
        serializer(this->m_gas_tracers);
        serializer(this->m_env_tracers);
        serializer(this->diffusion_control);
        serializer(this->max_iter);
        serializer(this->min_iter);
    }

    static Tracers serializationTestObject();
    bool operator==(const Tracers& data) const;

private:
    long long m_oil_tracers{};
    long long m_water_tracers{};
    long long m_gas_tracers{};
    long long m_env_tracers{};
    bool diffusion_control{false};
    long long max_iter{};
    long long min_iter{};
    // The TRACERS keyword has some additional options which seem quite arcane,
    // for now not included here.
};


class Runspec {
public:
    Runspec() = default;
    explicit Runspec( const Deck& );

    static Runspec serializationTestObject();

    std::time_t start_time() const noexcept;
    const UDQParams& udqParams() const noexcept;
    const Phases& phases() const noexcept;
    const Tabdims&  tabdims() const noexcept;
    const Regdims&  regdims() const noexcept;
    const EndpointScaling& endpointScaling() const noexcept;
    const Welldims& wellDimensions() const noexcept;
    const WellSegmentDims& wellSegmentDimensions() const noexcept;
    const NetworkDims& networkDimensions() const noexcept;
    const AquiferDimensions& aquiferDimensions() const noexcept;
    long long eclPhaseMask( ) const noexcept;
    const EclHysterConfig& hysterPar() const noexcept;
    const Actdims& actdims() const noexcept;
    const SatFuncControls& saturationFunctionControls() const noexcept;
    const Nupcol& nupcol() const noexcept;
    const Tracers& tracers() const;
    bool compositionalMode() const;
    size_t numComps() const;
    bool co2Storage() const noexcept;
    bool co2Sol() const noexcept;
    bool h2Sol() const noexcept;
    bool h2Storage() const noexcept;
    bool micp() const noexcept;
    bool mech() const noexcept;
    bool temp() const noexcept;
    bool compositional() const noexcept;

    bool operator==(const Runspec& data) const;
    static bool rst_cmp(const Runspec& full_state, const Runspec& rst_state);

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(this->m_start_time);
        serializer(active_phases);
        serializer(m_tabdims);
        serializer(m_regdims);
        serializer(endscale);
        serializer(welldims);
        serializer(wsegdims);
        serializer(netwrkdims);
        serializer(aquiferdims);
        serializer(udq_params);
        serializer(hystpar);
        serializer(m_actdims);
        serializer(m_sfuncctrl);
        serializer(m_nupcol);
        serializer(m_tracers);
        serializer(m_comps);
        serializer(m_co2storage);
        serializer(m_co2sol);
        serializer(m_h2sol);
        serializer(m_h2storage);
        serializer(m_micp);
        serializer(m_mech);
        serializer(m_temp);
    }

private:
    std::time_t m_start_time{};
    Phases active_phases{};
    Tabdims m_tabdims{};
    Regdims m_regdims{};
    EndpointScaling endscale{};
    Welldims welldims{};
    WellSegmentDims wsegdims{};
    NetworkDims netwrkdims{};
    AquiferDimensions aquiferdims{};
    UDQParams udq_params{};
    EclHysterConfig hystpar{};
    Actdims m_actdims{};
    SatFuncControls m_sfuncctrl{};
    Nupcol m_nupcol{};
    Tracers m_tracers{};
    size_t m_comps = 0;
    bool m_co2storage{false};
    bool m_co2sol{false};
    bool m_h2sol{false};
    bool m_h2storage{false};
    bool m_micp{false};
    bool m_mech{false};
    bool m_temp{false};
};

std::size_t declaredMaxRegionID(const Runspec& rspec);

} // namespace Opm

#endif // OPM_RUNSPEC_HPP
