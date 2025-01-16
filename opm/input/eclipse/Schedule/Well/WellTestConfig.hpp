/*
  Copyright 2018 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

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
#ifndef WELLTEST_CONFIG_H
#define WELLTEST_CONFIG_H

#include <cstddef>
#include <string>
#include <unordered_map>


namespace Opm {

namespace RestartIO {
struct RstState;
}

namespace WTest {
/*
  Different numerical values are used in the restart file to enumarate the
  possible WTEST modes and the actual reason a well has been closed.
*/

namespace EclConfigReason {
constexpr long long NONE      =  1;
constexpr long long PHYSICAL   = 2;
constexpr long long ECONOMIC   = 3;
constexpr long long GCON       = 5;
constexpr long long THPLimit   = 7;
constexpr long long CONNECTION = 11;
}

namespace EclCloseReason {
constexpr long long NONE     = 1; // May be written to UNRST during history
constexpr long long PHYSICAL = 3;
constexpr long long ECONOMIC = 5;
constexpr long long GCON     = 6;
constexpr long long THPLimit = 9;
}

enum class Reason {
    NONE     = 0,
    PHYSICAL = 1,
    ECONOMIC = 2,
    GROUP = 4,
    THP_DESIGN=8,
    COMPLETION=16,
};

}

class WellTestConfig {

public:
    using Reason = WTest::Reason;
    struct WTESTWell {
        std::string name{};
        long long reasons{};
        double test_interval{};
        long long num_test{};
        double startup_time{};
        // the related WTEST keywords is entered and will begin
        // taking effects since this report step
        long long begin_report_step{};

        bool operator==(const WTESTWell& data) const {
            return name == data.name &&
                   reasons == data.reasons &&
                   test_interval == data.test_interval &&
                   num_test == data.num_test &&
                   startup_time == data.startup_time &&
                   begin_report_step == data.begin_report_step;
        }

        WTESTWell() = default;
        WTESTWell(const std::string& name, long long reasons, double test_interval, long long num_test, double startup_time, long long begin_report_step);
        bool test_well(long long num_attempt, double elapsed) const;

        static long long inverse_ecl_reasons(long long ecl_reasons);
        static WTESTWell serializationTestObject();
        long long ecl_reasons() const;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(name);
            serializer(reasons);
            serializer(test_interval);
            serializer(num_test);
            serializer(startup_time);
            serializer(begin_report_step);
        }
    };

    static WellTestConfig serializationTestObject();

    WellTestConfig() = default;
    WellTestConfig(const RestartIO::RstState& rst_state, long long report_step);
    void add_well(const std::string& well, long long reasons, double test_interval,
                  long long num_test, double startup_time, long long current_step);
    void add_well(const std::string& well, const std::string& reasons, double test_interval,
                  long long num_test, double startup_time, long long current_step);
    void drop_well(const std::string& well);
    bool has(const std::string& well) const;
    bool has(const std::string& well, Reason reason) const;
    const WTESTWell& get(const std::string& well) const;

    static std::string reasonToString(const Reason reason);
    bool empty() const;

    bool operator==(const WellTestConfig& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(wells);
    }

private:
    std::unordered_map<std::string, WTESTWell> wells;
};
}

#endif

