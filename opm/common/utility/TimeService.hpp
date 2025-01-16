/*
  Copyright 2019 Equinor ASA.

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

#ifndef OPM_TIMESERVICE_HEADER_INCLUDED
#define OPM_TIMESERVICE_HEADER_INCLUDED

#include <chrono>
#include <ctime>
#include <string>
#include <unordered_map>

namespace Opm {

    class DeckRecord;

    using time_point = std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<int64_t, std::ratio<1,1000>>>;

    namespace TimeService {
    std::time_t to_time_t(const time_point& tp);
    time_point from_time_t(std::time_t t);
    time_point now();

    std::time_t advance(const std::time_t tp, const double sec);
    std::time_t makeUTCTime(std::tm timePoint);
    const std::unordered_map<std::string , long long>& eclipseMonthIndices();
    const std::unordered_map<long long, std::string>& eclipseMonthNames();
    long long eclipseMonth(const std::string& name);
    bool valid_month(const std::string& month_name);

    std::time_t mkdatetime(long long in_year, long long in_month, long long in_day, long long hour, long long minute, long long second);
    std::time_t mkdate(long long in_year, long long in_month, long long in_day);
    std::time_t portable_timegm(const std::tm* t);
    std::time_t timeFromEclipse(const DeckRecord &dateRecord);
    }

    class TimeStampUTC
    {
    public:
        struct YMD {
            long long year{0};
            long long month{0};
            long long day{0};

            bool operator==(const YMD& data) const
            {
                return year == data.year &&
                       month == data.month &&
                       day == data.day;
            }

            template<class Serializer>
            void serializeOp(Serializer& serializer)
            {
                serializer(year);
                serializer(month);
                serializer(day);
            }
        };

        TimeStampUTC() = default;

        explicit TimeStampUTC(const std::time_t tp);
        explicit TimeStampUTC(const YMD& ymd);
        TimeStampUTC(long long year, long long month, long long day);
        TimeStampUTC(const YMD& ymd,
                     long long hour,
                     long long minutes,
                     long long seconds,
                     long long usec);

        TimeStampUTC& operator=(const std::time_t tp);
        bool operator==(const TimeStampUTC& data) const;

        TimeStampUTC& hour(const long long h);
        TimeStampUTC& minutes(const long long m);
        TimeStampUTC& seconds(const long long s);
        TimeStampUTC& microseconds(const long long us);

        const YMD& ymd() const { return ymd_; }
        long long year()         const { return this->ymd_.year;  }
        long long month()        const { return this->ymd_.month; }
        long long day()          const { return this->ymd_.day;   }
        long long hour()         const { return this->hour_;      }
        long long minutes()      const { return this->minutes_;   }
        long long seconds()      const { return this->seconds_;   }
        long long microseconds() const { return this->usec_;      }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(ymd_);
            serializer(hour_);
            serializer(minutes_);
            serializer(seconds_);
            serializer(usec_);
        }

    private:

        YMD ymd_{};
        long long hour_{0};
        long long minutes_{0};
        long long seconds_{0};
        long long usec_{0};
    };

    TimeStampUTC operator+(const TimeStampUTC& lhs, std::chrono::duration<double> delta);
    std::time_t asTimeT(const TimeStampUTC& tp);
    std::time_t asLocalTimeT(const TimeStampUTC& tp);
    time_point asTimePoint(const TimeStampUTC& tp);


} // namespace Opm

#endif // OPM_TIMESERVICE_HEADER_INCLUDED
