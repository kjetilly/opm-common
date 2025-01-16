/*
  Copyright 2020 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  OPM is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along
  with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RST_HEADER
#define RST_HEADER

#include <opm/input/eclipse/EclipseState/Runspec.hpp>

#include <cstddef>
#include <ctime>
#include <utility>
#include <vector>

namespace Opm {

class UnitSystem;

} // namespace Opm

namespace Opm::RestartIO {

struct RstHeader
{
    RstHeader(const Runspec& runspec,
              const UnitSystem& unit_system,
              const std::vector<long long>& intehead,
              const std::vector<bool>& logihead,
              const std::vector<double>& doubhead);

    Runspec runspec;
    long long nx;
    long long ny;
    long long nz;
    long long nactive;
    long long num_wells;
    long long ncwmax;
    long long max_wells_in_group;
    long long max_groups_in_field;
    long long max_wells_in_field;
    long long year;
    long long month;
    long long mday;
    long long hour;
    long long minute;
    long long microsecond;
    long long phase_sum;
    long long niwelz;
    long long nswelz;
    long long nxwelz;
    long long nzwelz;
    long long niconz;
    long long nsconz;
    long long nxconz;
    long long nigrpz;
    long long nsgrpz;
    long long nxgrpz;
    long long nzgrpz;
    long long ncamax;
    long long niaaqz;
    long long nsaaqz;
    long long nxaaqz;
    long long nicaqz;
    long long nscaqz;
    long long nacaqz;
    long long tstep;
    long long report_step;
    long long histctl_override;
    long long newtmx;
    long long newtmn;
    long long litmax;
    long long litmin;
    long long mxwsit;
    long long mxwpit;
    long long version;
    long long iprog;
    long long nsegwl;
    long long nswlmx;
    long long nsegmx;
    long long nlbrmx;
    long long nisegz;
    long long nrsegz;
    long long nilbrz;
    long long ntfip ;
    long long nmfipr;
    long long ngroup;
    long long nwgmax;
    long long nfield_udq;
    long long ngroup_udq;
    long long nsegment_udq;
    long long nwell_udq;
    long long num_action;
    long long guide_rate_nominated_phase;
    long long max_wlist;

    bool e300_radial;
    bool e100_radial;
    bool enable_hysteris;
    bool enable_msw;
    bool is_live_oil;
    bool is_wet_gas;
    bool const_comp_oil;
    bool dir_relperm;
    bool reversible_relperm;
    bool endscale;
    bool dir_eps;
    bool reversible_eps;
    bool alt_eps;
    bool group_control_active;
    bool glift_all_nupcol;

    double next_timestep1;
    double next_timestep2;
    double max_timestep;
    double guide_rate_a;
    double guide_rate_b;
    double guide_rate_c;
    double guide_rate_d;
    double guide_rate_e;
    double guide_rate_f;
    double guide_rate_delay;
    double guide_rate_damping;
    double udq_range;
    double udq_undefined;
    double udq_eps;
    double glift_min_wait;
    double glift_rate_delta;
    double glift_min_eco_grad;

    std::time_t sim_time() const;
    std::pair<std::time_t, std::size_t> restart_info() const;
    long long num_udq() const;
};

} // namespace Opm::RestartIO

#endif  // RST_HEADER
