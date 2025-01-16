/*
  Copyright 2014 Statoil ASA.

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

#ifndef BOX_HPP_
#define BOX_HPP_

#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace Opm {
    class DeckRecord;
}

namespace Opm
{
    class Box
    {
    public:
        using IsActive = std::function<bool(const std::size_t globalIdx)>;
        using ActiveIdx = std::function<std::size_t(const std::size_t globalIdx)>;

        struct cell_index
        {
            std::size_t global_index;
            std::size_t active_index;
            std::size_t data_index;

            cell_index(std::size_t g,std::size_t a, std::size_t d)
                : global_index(g)
                , active_index(a)
                , data_index(d)
            {}

            // This constructor should is used by the global_index_list() member
            // which will return a list of *all* the cells in the box. In this
            // case the active_index will be set to the global_index. This is a
            // hack to simplify the treatment of global fields in the FieldProps
            // implementation.
            cell_index(std::size_t g, std::size_t d)
                : global_index(g)
                , active_index(g)
                , data_index(d)
            {}
        };

        explicit Box(const GridDims& gridDims,
                     IsActive        isActive,
                     ActiveIdx       activeIdx);

        Box(const GridDims& gridDims,
            IsActive        isActive,
            ActiveIdx       activeIdx,
            long long i1, long long i2,
            long long j1, long long j2,
            long long k1, long long k2);

        void update(const DeckRecord& deckRecord);
        void reset();

        bool isGlobal() const;
        std::size_t size() const;
        std::size_t getDim(std::size_t idim) const;

        const std::vector<cell_index>& index_list() const;
        const std::vector<cell_index>& global_index_list() const;

        bool operator==(const Box& other) const;
        bool equal(const Box& other) const;

        long long I1() const;
        long long I2() const;
        long long J1() const;
        long long J2() const;
        long long K1() const;
        long long K2() const;

    private:
        GridDims m_globalGridDims_{};
        IsActive m_globalIsActive_{};
        ActiveIdx m_globalActiveIdx_{};

        std::array<std::size_t, 3> m_dims{};
        std::array<std::size_t, 3> m_offset{};

        std::vector<cell_index> m_active_index_list;
        std::vector<cell_index> m_global_index_list;

        void init(long long i1, long long i2, long long j1, long long j2, long long k1, long long k2);
        void initIndexList();
        long long lower(long long dim) const;
        long long upper(long long dim) const;
    };
}


#endif
