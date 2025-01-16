/*
  Copyright 2019 Equinor ASA.

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

#ifndef OPM_IO_ECLUTIL_HPP
#define OPM_IO_ECLUTIL_HPP

#include <opm/io/eclipse/EclIOdata.hpp>

#include <cstdint>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

namespace Opm { namespace EclIO {

    long long flipEndianInt(long long num);
    std::int64_t flipEndianLongInt(std::int64_t num);
    float flipEndianFloat(float num);
    double flipEndianDouble(double num);
    bool isEOF(std::fstream* fileH);
    bool fileExists(const std::string& filename);
    bool isFormatted(const std::string& filename);
    bool is_number(const std::string& numstr);

    /// Compute the linearly combined summary vector ID number from two
    /// constituents.
    ///
    /// Constituents are *typically* one-based region IDs, but at least one
    /// of the two could be a component ID too.
    long long combineSummaryNumbers(const long long n1, const long long n2);

    /// Split a combined summary vector ID ('NUMS' entry) into its original
    /// two constituent IDs.
    std::tuple<long long,long long> splitSummaryNumber(const long long n);

    std::tuple<long long, long long> block_size_data_binary(eclArrType arrType);
    std::tuple<long long, long long, long long> block_size_data_formatted(eclArrType arrType);

    std::string trimr(const std::string &str1);

    std::uint64_t sizeOnDiskBinary(std::int64_t num, Opm::EclIO::eclArrType arrType, long long elementSize);
    std::uint64_t sizeOnDiskFormatted(const std::int64_t num, Opm::EclIO::eclArrType arrType, long long elementSize);

    void readBinaryHeader(std::fstream& fileH, std::string& tmpStrName,
                      long long& tmpSize, std::string& tmpStrType);

    void readBinaryHeader(std::fstream& fileH, std::string& arrName,
                      std::int64_t& size, Opm::EclIO::eclArrType &arrType, long long& elementSize);

    void readFormattedHeader(std::fstream& fileH, std::string& arrName,
                      std::int64_t &num, Opm::EclIO::eclArrType &arrType, long long& elementSize);

    template<typename T, typename T2>
    std::vector<T> readBinaryArray(std::fstream& fileH, const std::int64_t size, Opm::EclIO::eclArrType type,
                               std::function<T(T2)>& flip, long long elementSize);

    std::vector<long long> readBinaryInteArray(std::fstream &fileH, const std::int64_t size);
    std::vector<float> readBinaryRealArray(std::fstream& fileH, const std::int64_t size);
    std::vector<double> readBinaryDoubArray(std::fstream& fileH, const std::int64_t size);
    std::vector<bool> readBinaryLogiArray(std::fstream &fileH, const std::int64_t size);
    std::vector<size_t> readBinaryRawLogiArray(std::fstream &fileH, const std::int64_t size);
    std::vector<std::string> readBinaryCharArray(std::fstream& fileH, const std::int64_t size);
    std::vector<std::string> readBinaryC0nnArray(std::fstream& fileH, const std::int64_t size, long long elementSize);

    template<typename T>
    std::vector<T> readFormattedArray(const std::string& file_str, const long long size, std::int64_t fromPos,
                                       std::function<T(const std::string&)>& process);

    std::vector<long long> readFormattedInteArray(const std::string& file_str, const std::int64_t size, std::int64_t fromPos);

    std::vector<std::string> readFormattedCharArray(const std::string& file_str, const std::int64_t size,
                                                    std::int64_t fromPos, long long elementSize);

    std::vector<float> readFormattedRealArray(const std::string& file_str, const std::int64_t size, std::int64_t fromPos);
    std::vector<std::string> readFormattedRealRawStrings(const std::string& file_str, const std::int64_t size, std::int64_t fromPos);

    std::vector<bool> readFormattedLogiArray(const std::string& file_str, const std::int64_t size, std::int64_t fromPos);
    std::vector<double> readFormattedDoubArray(const std::string& file_str, const std::int64_t size, std::int64_t fromPos);

}} // namespace Opm::EclIO

#endif // OPM_IO_ECLUTIL_HPP
