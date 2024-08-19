/*
  Copyright 2014 Andreas Lauser

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
#include "config.h"

#define BOOST_TEST_MODULE ReadSPE11C
#include <boost/test/unit_test.hpp>

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/output/data/Cells.hpp>

#include <opm/io/eclipse/EGrid.hpp>
#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/EclFile.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <fmt/format.h>
#include <fstream>
#include <ios>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <time.h>

#include <tests/WorkArea.hpp>

using namespace Opm;

namespace
{

bool
keywordExists(const std::vector<EclIO::EclFile::EclEntry>& knownVec, const std::string& arrayname)
{
    return std::any_of(knownVec.begin(), knownVec.end(), [&arrayname](const EclIO::EclFile::EclEntry& entry) -> bool {
        return std::get<0>(entry) == arrayname;
    });
}

template <typename T>
T
sum(const std::vector<T>& array)
{
    return std::accumulate(array.begin(), array.end(), T(0));
}

template <typename T, typename U>
void
compareErtData(const std::vector<T>& src, const std::vector<U>& dst, double tolerance)
{
    BOOST_REQUIRE_EQUAL(src.size(), dst.size());

    for (size_t i = 0; i < src.size(); ++i)
        BOOST_CHECK_CLOSE(src[i], dst[i], tolerance);
}

void
compareErtData(const std::vector<int>& src, const std::vector<int>& dst)
{
    BOOST_CHECK_EQUAL_COLLECTIONS(src.begin(), src.end(), dst.begin(), dst.end());
}

void
checkEgridFile(const EclipseGrid& eclGrid, const std::string& gridFilename)
{
    auto egridFile = EclIO::EGrid(gridFilename);
    std::cout << "In file: " << std::endl;
    for (const auto& arrayname : egridFile.arrayNames()) {
        std::cout << arrayname << std::endl;
    }
    std::cout << "#####################" << std::endl;
    // {
    //     std::cout << "Getting coord" << std::endl;
    //     for (const auto& arrayname : egridFile.arrayNames()) {
    //         std::cout << arrayname << std::endl;
    //     }
    //     const auto& coord = egridFile.get<float>("COORD");
    //     const auto& expect = eclGrid.getCOORD();
    //     // compareErtData(expect, coord, 1e-6);
    // }

    // {
    //     const auto& zcorn = egridFile.get<float>("ZCORN");
    //     const auto& expect = eclGrid.getZCORN();
    //     // compareErtData(expect, zcorn, 1e-6);
    // }

    // if (egridFile.hasKey("ACTNUM")) {
    //     const auto& actnum = egridFile.get<int>("ACTNUM");
    //     auto expect = eclGrid.getACTNUM();

    //     if (expect.empty()) {
    //         const auto numCells = eclGrid.getNX() * eclGrid.getNY() * eclGrid.getNZ();
    //         expect.assign(numCells, 1);
    //     }

    //     //        compareErtData(expect, actnum);
    // }
}

void
checkInitFile(const Deck& deck, const std::string& initpath)
{
    EclIO::EclFile initFile {initpath};

    if (initFile.hasKey("PORO")) {
        const auto& poro = initFile.get<float>("PORO");
        const auto& expect = deck["PORO"].back().getSIDoubleData();

        compareErtData(expect, poro, 1e-4);
    }

    if (initFile.hasKey("PERMX")) {
        const auto& expect = deck["PERMX"].back().getSIDoubleData();
        auto permx = initFile.get<float>("PERMX");

        for (auto& kx : permx) {
            kx *= 9.869233e-16;
        }

        compareErtData(expect, permx, 1e-4);
    }

    // These arrays should always be in the INIT file, irrespective of
    // keyword presence in the inut deck.
    BOOST_CHECK_MESSAGE(initFile.hasKey("NTG"), R"(INIT file must have "NTG" array)");
    BOOST_CHECK_MESSAGE(initFile.hasKey("FIPNUM"), R"(INIT file must have "FIPNUM" array)");
    BOOST_CHECK_MESSAGE(initFile.hasKey("SATNUM"), R"(INIT file must have "SATNUM" array)");
    std::array<std::string, 3> multipliers {"MULTX", "MULTY", "MULTZ"};
    // for (const auto& mult : multipliers) {
    //     BOOST_CHECK_MESSAGE(initFile.hasKey(mult), R"(INIT file must have ")" + mult + R"(" array)");
    // }
    std::cout << "init file names" << std::endl;
    for (const auto& arrayname : initFile.arrayNames()) {
        std::cout << arrayname << std::endl;
    }
    std::cout << "############################" << std::endl;

    // for (const auto& prop : simProps) {
    //     BOOST_CHECK_MESSAGE(initFile.hasKey(prop.first), R"(INIT file must have ")" + prop.first + R"(" array)");
    // }
}

void
checkRestartFile(int timeStepIdx, const std::string& rstFilename)
{
    EclIO::ERst rstFile {rstFilename};

    for (int i = 1; i <= timeStepIdx; ++i) {
        if (!rstFile.hasReportStepNumber(i)) {
            std::cout << "skipping " << i << std::endl;
            continue;
        }



        rstFile.loadReportStepNumber(i);

        const auto& knownVec = rstFile.listOfRstArrays(i);

        if (keywordExists(knownVec, "PRESSURE")) {
            std::cout << "Reading pressure " << std::endl;
            const auto& press = rstFile.getRestartData<float>("PRESSURE", i, 0);
            std::cout << "Sum press = " << sum(press) << std::endl;
            std::cout << "press.size() = " << press.size() << std::endl;
        }

        if (keywordExists(knownVec, "SWAT")) {
            std::cout << "Reading swat " << std::endl;
            const auto& swat = rstFile.getRestartData<float>("SWAT", i, 0);

            std::cout << "Sum swat = " << sum(swat) << std::endl;
        }

        if (keywordExists(knownVec, "SGAS")) {
            std::cout << "Reading sgas " << std::endl;
            const auto& sgas = rstFile.getRestartData<float>("SGAS", i, 0);

            std::cout << "Sum sgas = " << sum(sgas) << std::endl;
        }

        if (keywordExists(knownVec, "KRO")) {
            std::cout << "Reading kro " << std::endl;
            const auto& kro = rstFile.getRestartData<float>("KRO", i, 0);
        }

        if (keywordExists(knownVec, "KRG")) {
            std::cout << "Reading krg " << std::endl;
            const auto& krg = rstFile.getRestartData<float>("KRG", i, 0);
        }
    }
}


} // Anonymous namespace

BOOST_AUTO_TEST_CASE(ReadLargeSPE11C)
{
    // TODO: Use filesystem lib
    const std::string basedirectory = "/mnt/external/kjetil/alldata/spe11_full_run";
    const std::string deckfilename = fmt::format("{}/deck/SPE11C.DATA", basedirectory);
    const std::string initFilename = fmt::format("{}/output/SPE11C.INIT", basedirectory);
    const std::string unrstFilename = fmt::format("{}/output/SPE11C.UNRST", basedirectory);
    const size_t nx = 466;
    const size_t ny = 466;
    const size_t nz = 466;
    const size_t totalcells = nx * ny * nz;
    const size_t nxny = nx * ny;

    std::ifstream deckfile(deckfilename);
    std::string deckString((std::istreambuf_iterator<char>(deckfile)), std::istreambuf_iterator<char>());

    const auto deck = Parser().parseFile(deckfilename);
    auto es = EclipseState(deck);
    const auto& eclGrid = es.getInputGrid();

    checkInitFile(deck, initFilename);
    checkEgridFile(eclGrid, initFilename);
    checkRestartFile(30, unrstFilename);
}
