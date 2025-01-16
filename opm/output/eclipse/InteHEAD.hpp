/*
  Copyright 2021-2024 Equinor ASA.
  Copyright 2016, 2017, 2018 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_INTEHEAD_HEADER_INCLUDED
#define OPM_INTEHEAD_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>

#include <array>
#include <cstddef>
#include <ctime>
#include <memory>
#include <vector>

namespace Opm {

class EclipseGrid;
class EclipseState;
class Phases;
class Schedule;
class ScheduleState;
class UDQInput;
class UnitSystem;

} // namespace Opm

namespace Opm { namespace RestartIO {

    class InteHEAD
    {
    public:
        struct WellTableDim {
            long long numWells{};
            long long maxPerf{};
            long long maxWellInGroup{};
            long long maxGroupInField{};
            long long maxWellsInField{};
            long long mxwlstprwel{};
            long long mxdynwlst{};
        };

        struct WellSegDims {
            long long nsegwl{};
            long long nswlmx{};
            long long nsegmx{};
            long long nlbrmx{};
            long long nisegz{};
            long long nrsegz{};
            long long nilbrz{};
        };

        struct RegDims {
            long long ntfip{};
            long long nmfipr{};
            long long nrfreg{};
            long long ntfreg{};
            long long nplmix{};
        };

        struct RockOpts {
            long long ttyp{};
        };

        struct TimePoint {
            long long year{};
            long long month{};          // 1..12
            long long day{};            // 1..31

            long long hour{};           // 0..23
            long long minute{};         // 0..59
            long long second{};         // 0..59

            long long microseconds{};   // 0..999999
        };

        struct Phases {
            long long oil{};
            long long water{};
            long long gas{};
        };

        struct TuningPar {
            long long newtmx{};
            long long newtmn{};
            long long litmax{};
            long long litmin{};
            long long mxwsit{};
            long long mxwpit{};
            long long wseg_mx_rst{};
        };

        struct Group {
            long long ngroups{};
        };

        struct UdqParam {
            long long udqParam_1{};
            long long num_iuads{};
            long long num_iuaps{};

            std::array<long long, static_cast<std::size_t>(UDQVarType::NumTypes)> numUDQs{};
        };

        struct ActionParam {
            long long no_actions{};
            long long max_no_sched_lines_per_action{};
            long long max_no_conditions_per_action{};
            long long max_no_characters_per_line{};
        };

        struct GuideRateNominatedPhase {
            long long nominated_phase;
        };

        struct ActiveNetwork {
            long long actnetwrk;
        };

        struct NetworkDims {
            long long noactnod{};
            long long noactbr{};
            long long nodmax{};
            long long nbrmax{};
            long long nibran{};
            long long nrbran{};
            long long ninode{};
            long long nrnode{};
            long long nznode{};
            long long ninobr{};
        };

        struct NetBalanceDims {
            long long maxNoIterationsNBC{};
            long long maxNoIterationsTHP{};
        };

        struct AquiferDims {
            // Number of active analytic aquifers (# unique aquifer IDs)
            long long numAquifers {0};

            // Declared maximum number of analytic aquifers in model
            // (AQUDIMS(5))
            long long maxNumAquifers {0};

            // Declared maximum number of connections in any analytic
            // aquifer (AQUDIMS(6))
            long long maxNumAquiferConn {0};

            // Maximum number of *active* connections in any analytic aquifer
            long long maxNumActiveAquiferConn {0};

            // Maximum aquifer ID across all of the model's analytic aquifers.
            long long maxAquiferID {0};

            // Number of numeric aquifer records (lines of AQUNUM data, AQUDIMS(1))
            long long numNumericAquiferRecords {0};

            // Number of data elements per aquifer in IAAQ array.
            long long numIntAquiferElem {18};

            // Number of data elements per aquifer in SAAQ array.
            long long numRealAquiferElem {24};

            // Number of data elements per aquifer in XAAQ array.
            long long numDoubAquiferElem {10};

            // Number of data elements in IAQN array per numeric aquifer record.
            long long numNumericAquiferIntElem {10};

            // Number of data elements in RAQN array per numeric aquifer record.
            long long numNumericAquiferDoubleElem {13};

            // Number of data elements per coonnection in ICAQ array.
            long long numIntConnElem {7};

            // Number of data elements per connecetion in SCAQ array.
            long long numRealConnElem {2};

            // Number of data elements per connection in ACAQ array.
            long long numDoubConnElem {4};
        };

        InteHEAD();
        ~InteHEAD() = default;

        InteHEAD(const InteHEAD& rhs) = default;
        InteHEAD(InteHEAD&& rhs) = default;

        InteHEAD& operator=(const InteHEAD& rhs) = default;
        InteHEAD& operator=(InteHEAD&& rhs) = default;

        InteHEAD& dimensions(const long long nx, const long long ny, const long long nz);
        InteHEAD& dimensions(const std::array<long long,3>& cartDims);
        InteHEAD& numActive(const long long nactive);

        InteHEAD& unitConventions(const UnitSystem& usys);
        InteHEAD& wellTableDimensions(const WellTableDim& wtdim);
        InteHEAD& aquiferDimensions(const AquiferDims& aqudims);

        InteHEAD& calendarDate(const TimePoint& date);
        InteHEAD& activePhases(const Phases& phases);

        InteHEAD& drsdt(const Schedule&   sched,
                        const std::size_t lookup_step);

        InteHEAD& params_NWELZ(const long long niwelz, const long long nswelz, const long long nxwelz, const long long nzwelz);
        InteHEAD& params_NCON(const long long niconz, const long long nsconz, const long long nxconz);
        InteHEAD& params_GRPZ(const std::array<long long, 4>& grpz);
        InteHEAD& params_NGCTRL(const long long gct);

        InteHEAD& stepParam(const long long tstep, const long long report_step);
        InteHEAD& tuningParam(const TuningPar& tunpar);
        InteHEAD& variousParam(const long long version, const long long iprog);
        InteHEAD& wellSegDimensions(const WellSegDims& wsdim);
        InteHEAD& activeNetwork(const ActiveNetwork& actntwrk);
        InteHEAD& networkDimensions(const NetworkDims& nwdim);
        InteHEAD& netBalanceData(const NetBalanceDims& nwbaldim);
        InteHEAD& regionDimensions(const RegDims& rdim);
        InteHEAD& rockOpts(const RockOpts& rckop);
        InteHEAD& ngroups(const Group& gr);
        InteHEAD& udqParam_1(const UdqParam& udqpar);
        InteHEAD& actionParam(const ActionParam& act_par);
        InteHEAD& variousUDQ_ACTIONXParam();
        InteHEAD& nominatedPhaseGuideRate(GuideRateNominatedPhase nphase);
        InteHEAD& whistControlMode(long long mode);
        InteHEAD& liftOptParam(long long in_enc);

        static long long numRsegElem(const Opm::Phases& phase);

        const std::vector<long long>& data() const
        {
            return this->data_;
        }

    private:
        std::vector<long long> data_;
    };

    InteHEAD::TimePoint
    getSimulationTimePoint(const std::time_t start,
                           const double      elapsed);

    InteHEAD::AquiferDims
    inferAquiferDimensions(const EclipseState& es);

    InteHEAD::AquiferDims
    inferAquiferDimensions(const EclipseState&  es,
                           const ScheduleState& sched);
}} // Opm::RestartIO

#endif // OPM_INTEHEAD_HEADER_INCLUDED
