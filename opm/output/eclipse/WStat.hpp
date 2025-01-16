/*
  Copyright 2021 Equinor ASA.

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
#ifndef WSTAT_HPP
#define WSTAT_HPP

#include <string>

namespace Opm {

namespace WStat {
namespace numeric {
constexpr long long UNKNOWN = 0;
constexpr long long PROD    = 1;
constexpr long long INJ     = 2;
constexpr long long SHUT    = 3;
constexpr long long STOP    = 4;
constexpr long long PSHUT   = 5;
constexpr long long PSTOP   = 6;
}

namespace symbolic {
const std::string UNKNOWN = "UNKNOWN";
const std::string PROD    = "PROD";
const std::string INJ     = "INJ";
const std::string SHUT    = "SHUT";
const std::string STOP    = "STOP";
const std::string PSHUT   = "PSHUT";
const std::string PSTOP   = "PSTOP";
}





}
}




#endif

