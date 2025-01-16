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

#ifndef OPM_IO_EGRID_HPP
#define OPM_IO_EGRID_HPP

#include <opm/io/eclipse/EclFile.hpp>

#include <array>
#include <filesystem>
#include <string>
#include <vector>
#include <map>

namespace Opm { namespace EclIO {

class EGrid : public EclFile
{
public:
    explicit EGrid(const std::string& filename, const std::string& grid_name = "global");

    long long global_index(long long i, long long j, long long k) const;
    long long active_index(long long i, long long j, long long k) const;

    const std::array<long long, 3>& dimension() const { return nijk; }

    std::array<long long, 3> ijk_from_active_index(long long actInd) const;
    std::array<long long, 3> ijk_from_global_index(long long globInd) const;

    void getCellCorners(long long globindex, std::array<double, 8>& X, std::array<double, 8>& Y, std::array<double, 8>& Z);
    void getCellCorners(const std::array<long long, 3>& ijk, std::array<double, 8>& X, std::array<double, 8>& Y, std::array<double, 8>& Z);

    std::vector<std::array<float, 3>> getXYZ_layer(long long layer, bool bottom=false);
    std::vector<std::array<float, 3>> getXYZ_layer(long long layer, const std::array<long long, 4>& box, bool bottom=false);

    long long activeCells() const { return nactive; }
    long long totalNumberOfCells() const { return nijk[0] * nijk[1] * nijk[2]; }

    void load_grid_data();
    void load_nnc_data();
    bool with_mapaxes() const { return m_mapaxes_loaded; }
    void mapaxes_transform(double& x, double& y) const;
    bool is_radial() const { return m_radial; }

    const std::vector<long long>& hostCellsGlobalIndex() const { return host_cells; }
    std::vector<std::array<long long, 3>> hostCellsIJK();

    // zero based: i1,j1,k1, i2,j2,k2, transmisibility
    using NNCentry = std::tuple<long long, long long, long long, long long, long long, long long, float>;
    std::vector<NNCentry> get_nnc_ijk();

    const std::vector<std::string>& list_of_lgrs() const { return lgr_names; }

    const std::array<double, 6>& get_mapaxes() const { return m_mapaxes; }
    const std::string& get_mapunits() const { return m_mapunits; }
    const std::vector<float>& get_coord() const { return coord_array; }
    const std::vector<float>& get_zcorn() const { return zcorn_array; }




private:
    std::filesystem::path inputFileName, initFileName;
    std::string m_grid_name;
    bool m_radial;

    std::array<double, 6> m_mapaxes;
    std::string m_mapunits;
    bool m_mapaxes_loaded;
    std::array<double, 4> origin;
    std::array<double, 2> unit_x;
    std::array<double, 2> unit_y;

    std::array<long long, 3> nijk;
    std::array<long long, 3> host_nijk;

    long long nactive;
    mutable bool m_nncs_loaded;

    std::vector<long long> act_index;
    std::vector<long long> glob_index;

    std::vector<float> coord_array;
    std::vector<float> zcorn_array;

    std::vector<long long> nnc1_array;
    std::vector<long long> nnc2_array;
    std::vector<float> transnnc_array;
    std::vector<long long> host_cells;
    std::map<long long,long long> res;
     
    std::vector<std::string> lgr_names;
    
    long long numres;
    
    long long zcorn_array_index;
    long long coord_array_index;
    long long coordsys_array_index;
    long long actnum_array_index;
    long long nnc1_array_index;
    long long nnc2_array_index;

    std::vector<float> get_zcorn_from_disk(long long layer, bool bottom);

    void getCellCorners(const std::array<long long, 3>& ijk, const std::vector<float>& zcorn_layer,
                        std::array<double, 4>& X, std::array<double, 4>& Y, std::array<double, 4>& Z);

    void mapaxes_init();
    
};

}} // namespace Opm::EclIO

#endif // OPM_IO_EGRID_HPP
