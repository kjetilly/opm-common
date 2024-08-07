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

#include <opm/input/eclipse/Schedule/UDQ/UDQInput.hpp>

#include <opm/input/eclipse/Schedule/UDQ/UDQAssign.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQDefine.hpp>

#include <stdexcept>
#include <string>
#include <variant>

namespace Opm {

UDQInput::UDQInput(const UDQIndex&    index_arg,
                   const UDQDefine*   udq_define,
                   const std::string& unit_arg)
    : index     (index_arg)
    , value     (udq_define)
    , m_keyword (udq_define->keyword())
    , m_var_type(udq_define->var_type())
    , m_unit    (unit_arg)
{}

UDQInput::UDQInput(const UDQIndex&    index_arg,
                   const UDQAssign*   udq_assign,
                   const std::string& unit_arg)
    : index     (index_arg)
    , value     (udq_assign)
    , m_keyword (udq_assign->keyword())
    , m_var_type(udq_assign->var_type())
    , m_unit    (unit_arg)
{}

template<>
bool UDQInput::is<UDQAssign>() const
{
    return std::holds_alternative<const UDQAssign*>(this->value);
}

template<>
bool UDQInput::is<UDQDefine>() const
{
    return std::holds_alternative<const UDQDefine*>(this->value);
}

template<>
const UDQAssign& UDQInput::get<UDQAssign>() const
{
    if (this->is<UDQAssign>()) {
        return *std::get<const UDQAssign*>(this->value);
    }

    throw std::runtime_error {
        "Requested UDQ assignment object from non-assignment container"
    };
}

template<>
const UDQDefine& UDQInput::get<UDQDefine>() const
{
    if (this->is<UDQDefine>()) {
        return *std::get<const UDQDefine*>(this->value);
    }

    throw std::runtime_error {
        "Requested UDQ definition object from non-definition container"
    };
}

bool UDQInput::operator==(const UDQInput& other) const
{
    const auto structure_okay =
        (this->value.index() == other.value.index()) &&
        (this->m_keyword == other.m_keyword)         &&
        (this->m_var_type == other.m_var_type)       &&
        (this->m_unit == other.m_unit)
        ;

    if (! structure_okay) { return false; }

    return (this->is<UDQDefine>())
        ? (this->get<UDQDefine>() == other.get<UDQDefine>())
        : (this->get<UDQAssign>() == other.get<UDQAssign>())
        ;
}

} // namespace Opm
