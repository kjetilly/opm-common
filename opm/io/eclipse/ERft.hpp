/*
   Copyright 2019 Equinor ASA.

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM.  If not, see <http://www.gnu.org/licenses/>.
   */

#ifndef OPM_IO_ERFT_HPP
#define OPM_IO_ERFT_HPP

#include <opm/io/eclipse/EclFile.hpp>

#include <ctime>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Opm { namespace EclIO {

class ERft : public EclFile
{
public:
    explicit ERft(const std::string &filename);

    using RftDate = std::tuple<long long,long long,long long>;
    template <typename T>
    const std::vector<T>& getRft(const std::string& name, const std::string& wellName,
                                 const RftDate& date) const;

    template <typename T>
    const std::vector<T>& getRft(const std::string& name, const std::string& wellName,
                                 long long year, long long month, long long day) const;
    template <typename T>
    const std::vector<T>& getRft(const std::string& name, long long reportIndex) const;

    std::vector<std::string> listOfWells() const;
    std::vector<RftDate> listOfdates() const;

    using RftReportList = std::vector<std::tuple<std::string, RftDate, float>>;
    const RftReportList& listOfRftReports() const { return rftReportList; }

    bool hasRft(const std::string& wellName, const RftDate& date) const;
    bool hasRft(const std::string& wellName, long long year, long long month, long long day) const;

    std::vector<EclEntry> listOfRftArrays(long long reportIndex ) const;

    std::vector<EclEntry> listOfRftArrays(const std::string& wellName,
                                          const RftDate& date) const;

    std::vector<EclEntry> listOfRftArrays(const std::string& wellName,
                                          long long year, long long month, long long day) const;

    bool hasArray(const std::string& arrayName, const std::string& wellName,
                  const RftDate& date) const;

    bool hasArray(const std::string& arrayName, long long reportInd) const;

    long long numberOfReports() { return numReports; }

private:
    std::map<long long, std::tuple<long long,long long>> arrIndexRange;
    long long numReports;
    std::vector<float> timeList;

    std::set<std::string> wellList;
    std::set<RftDate> dateList;
    RftReportList rftReportList;

    std::map<std::tuple<std::string,RftDate>,long long> reportIndices;  //  mapping report index to wellName and date (tupe)

    long long getReportIndex(const std::string& wellName, const RftDate& date) const;

    long long getArrayIndex(const std::string& name, long long reportIndex) const;
    long long getArrayIndex(const std::string& name, const std::string& wellName,
                      const RftDate& date) const;
};

}} // namespace Opm::EclIO

#endif // OPM_IO_ERFT_HPP
