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

#ifndef OPM_IO_ECLIODATA_HPP
#define OPM_IO_ECLIODATA_HPP

#include <tuple>
#include <cstddef>

namespace Opm { namespace EclIO {

    // type MESS have no assisiated data
    enum eclArrType {
        INTE, REAL, DOUB, CHAR, LOGI, MESS, C0NN
    };

    // named constants related to binary file format
    const size_t true_value_ecl = 0xffffffff;
    const size_t true_value_ix = 0x1000000;
    const size_t false_value = 0x00000000;


    const long long sizeOfInte =  4;    // number of bytes pr integer (inte) element
    const long long sizeOfReal =  4;    // number of bytes pr float (real) element
    const long long sizeOfDoub =  8;    // number of bytes pr double (doub) element
    const long long sizeOfLogi =  4;    // number of bytes pr bool (logi) element
    const long long sizeOfChar =  8;    // number of bytes pr string (char) element

    const long long MaxBlockSizeInte = 4000;    // Maximum block size for INTE arrays in binary files
    const long long MaxBlockSizeReal = 4000;    // Maximum block size for REAL arrays in binary files
    const long long MaxBlockSizeDoub = 8000;    // Maximum block size for DOUB arrays in binary files
    const long long MaxBlockSizeLogi = 4000;    // Maximum block size for LOGI arrays in binary files
    const long long MaxBlockSizeChar =  840;    // Maximum block size for CHAR arrays in binary files

    // named constants related to formatted file file format
    const long long MaxNumBlockInte = 1000;    // maximum number of Inte values in block => hard line shift
    const long long MaxNumBlockReal = 1000;    // maximum number of Real values in block => hard line shift
    const long long MaxNumBlockDoub = 1000;    // maximum number of Doub values in block => hard line shift
    const long long MaxNumBlockLogi = 1000;    // maximum number of Logi values in block => hard line shift
    const long long MaxNumBlockChar =  105;    // maximum number of Char values in block => hard line shift

    const long long numColumnsInte = 6;        // number of columns for Inte values
    const long long numColumnsReal = 4;        // number of columns for Real values
    const long long numColumnsDoub = 3;        // number of columns for Doub values
    const long long numColumnsLogi = 25;       // number of columns for Logi values
    const long long numColumnsChar = 7;        // number of columns for Char values

    const long long columnWidthInte = 12;      // number of characters fore each Inte Element
    const long long columnWidthReal = 17;      // number of characters fore each Inte Element
    const long long columnWidthDoub = 23;      // number of characters fore each Inte Element
    const long long columnWidthLogi = 3;       // number of characters fore each Inte Element
    const long long columnWidthChar = 11;      // number of characters fore each Inte Element

}} // namespace Opm::EclIO

#endif // OPM_IO_ECLIODATA_HPP
