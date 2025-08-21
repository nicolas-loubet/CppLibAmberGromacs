#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <filesystem>
#include <cmath>

// ===================== HELPERS =======================

void checkVectorClose(const std::vector<Real>& v1, const std::vector<Real>& v2, Real tol=1e-6) {
    REQUIRE(v1.size() == v2.size());
    for (size_t i=0;i<v1.size();++i) {
        REQUIRE(v1[i] == Approx(v2[i]).margin(tol));
    }
}

// ===================== AMBER TESTS =======================

TEST_CASE("Configuration - AMBER basic and interaction values", "[Configuration][AMBER]") {
    std::string prmtop = "../files/amber.prmtop";
    std::string pdb    = "../files/amber_frame_1.pdb";

    REQUIRE(std::filesystem::exists(prmtop));
    REQUIRE(std::filesystem::exists(pdb));

    TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
    TopolInfo topol = top_reader->readTopology(prmtop);

    CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    Configuration conf(coord_reader, pdb, topol);

    SECTION("Número de moléculas y límites") {
        REQUIRE(conf.getNMolec() == topol.num_molecules);
        REQUIRE(conf.getBounds().x > 0);
        REQUIRE(conf.getBounds().y > 0);
        REQUIRE(conf.getBounds().z > 0);
    }

    SECTION("Acceso a moléculas y aguas") {
        REQUIRE(conf.getMolec(1).getID() == 1);
        REQUIRE(conf.getWater(1000).isWater());
    }

    SECTION("Valores V1S-V4S de una molécula de referencia") {
        int mol_id = 1000;
        auto per_site = conf.getInteractionsPerSite(mol_id);

        std::vector<Real> expected = { -28.3249, -24.1428, -11.7532, -11.3781 }; 

        checkVectorClose(per_site, expected, 1e-3);

        REQUIRE(conf.v_4S(mol_id) == Approx(expected[3]).margin(1e-3));
    }

    SECTION("findNearby devuelve vecinos correctos") {
        auto neigh = conf.findNearby(1000, 5.0);
        REQUIRE(neigh.size > 0);
        for (int i=0;i<neigh.size;i++){
            REQUIRE(neigh.arr[i] != 1);
        }
    }

    delete top_reader;
    delete coord_reader;
}

// ===================== GROMACS TESTS =======================

TEST_CASE("Configuration - GROMACS interactions", "[Configuration][GROMACS]") {
    std::string top = "../files/gromacs.top";
    std::string gro = "../files/gromacs.gro";

    REQUIRE(std::filesystem::exists(top));
    REQUIRE(std::filesystem::exists(gro));

    TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS);
    TopolInfo topol = top_reader->readTopology(top);

    CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
    Configuration conf(coord_reader, gro, topol);

    SECTION("Número de moléculas y chequeo básico") {
        REQUIRE(conf.getNMolec() == topol.num_molecules);
        REQUIRE(conf.getMolec(1000).isWater());
    }

    SECTION("Valores V1S-V4S de molécula 1000") {
        int mol_id = 1000;
        auto per_site = conf.getInteractionsPerSite(mol_id);

        std::vector<Real> expected = { -35.6644, -29.2759, -26.3631, -25.9072 };

        checkVectorClose(per_site, expected, 1e-3);
        REQUIRE(conf.v_4S(mol_id) == Approx(expected[3]).margin(1e-3));
    }

    delete top_reader;
    delete coord_reader;
}

// ===================== LAMMPS TESTS =======================

TEST_CASE("Configuration - LAMMPS interactions", "[Configuration][LAMMPS]") {
    std::string data = "../files/lammps_system.data";
    std::string traj = "../files/lammps_traj_minimized_9A.xyz";

    REQUIRE(std::filesystem::exists(data));
    REQUIRE(std::filesystem::exists(traj));

    TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS);
    TopolInfo topol = top_reader->readTopology(data);

    CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
    Configuration conf(coord_reader, traj, topol);

    SECTION("Número de moléculas y chequeo básico") {
        REQUIRE(conf.getNMolec() == topol.num_molecules);
    }

    SECTION("Valores V1S-V4S de molécula 1000") {
        int mol_id = 1000;
        auto per_site = conf.getInteractionsPerSite(mol_id);

        std::vector<Real> expected = { -42.6962, -42.6758, -36.0950, -30.1109 };

        checkVectorClose(per_site, expected, 1e-3);
    }

    delete top_reader;
    delete coord_reader;
}
