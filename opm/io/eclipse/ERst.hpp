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

#ifndef OPM_IO_ERST_HPP
#define OPM_IO_ERST_HPP

#include <opm/io/eclipse/EclFile.hpp>

#include <ios>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>


namespace Opm { namespace EclIO { namespace OutputStream {
    class Restart;
}}}

namespace Opm { namespace EclIO {

class ERst : public EclFile
{
public:
    explicit ERst(const std::string& filename);

    bool hasReportStepNumber(long long number) const;
    bool hasArray(const std::string& name, long long number) const;
    bool hasLGR(const std::string& gridname, long long reportStepNumber) const;

    void loadReportStepNumber(long long number);

    template <typename T>
    const std::vector<T>& getRestartData(const std::string& name, long long reportStepNumber)
    {
        return getRestartData<T>(name,reportStepNumber, 0);
    }

    template <typename T>
    const std::vector<T>& getRestartData(const std::string& name, long long reportStepNumber, long long occurrence);

    template <typename T>
    const std::vector<T>& getRestartData(long long index, long long reportStepNumber)
    {
        auto indRange = this->getIndexRange(reportStepNumber);
        return  this->get<T>(index + std::get<0>(indRange));
    }

    template <typename T>
    const std::vector<T>& getRestartData(const std::string& name, long long reportStepNumber, const std::string& lgr_name);

    template <typename T>
    const std::vector<T>& getRestartData(long long index, long long reportStepNumber, const std::string& lgr_name);

    long long occurrence_count(const std::string& name, long long reportStepNumber) const;
    size_t numberOfReportSteps() const { return seqnum.size(); };

    const std::vector<long long>& listOfReportStepNumbers() const { return seqnum; }

    std::vector<EclEntry> listOfRstArrays(long long reportStepNumber);
    std::vector<EclEntry> listOfRstArrays(long long reportStepNumber, const std::string& lgr_name);

    friend class OutputStream::Restart;

private:
    long long nReports;
    std::vector<long long> seqnum;                           // report step numbers, from SEQNUM array in restart file
    mutable std::unordered_map<long long,bool> reportLoaded;
    std::map<long long, std::pair<long long,long long>> arrIndexRange;   // mapping report step number to array indeces (start and end)
    std::vector<std::vector<std::string>> lgr_names;                           // report step numbers, from SEQNUM array in restart file

    void initUnified();
    void initSeparate(const long long number);

    long long get_start_index_lgrname(long long number, const std::string& lgr_name);

    long long getArrayIndex(const std::string& name, long long seqnum, long long occurrence);
    long long getArrayIndex(const std::string& name, long long number, const std::string& lgr_name);

    std::tuple<long long,long long> getIndexRange(long long reportStepNumber) const;

    std::streampos
    restartStepWritePosition(const long long seqnumValue) const;

};

}} // namespace Opm::EclIO

#endif // OPM_IO_ERST_HPP
