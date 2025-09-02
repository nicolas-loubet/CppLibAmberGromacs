#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <filesystem>
#include <cmath>
#include <limits>

// ===================== HELPERS =======================

// Compare two vectors within tolerance
static void checkVectorClose(const std::vector<Real>& v1, const std::vector<Real>& v2, Real tol=1e-6) {
	REQUIRE(v1.size() == v2.size());
	for (size_t i=0;i<v1.size();++i) {
		REQUIRE(v1[i] == Approx(v2[i]).margin(tol));
	}
}

// Check increasing order
static void checkIncreasing(const std::vector<Real>& v) {
	for (size_t i=1;i<v.size();++i) {
		REQUIRE(v[i-1] <= v[i] - 1e-12);
	}
}

// ===================== AMBER TESTS (expanded) =======================

TEST_CASE("Configuration - AMBER: basic checks and interaction values (expanded)", "[Configuration][AMBER]") {
	std::string prmtop = "../files/amber.prmtop";
	std::string pdb	= "../files/amber_frame_1.pdb";

	REQUIRE(std::filesystem::exists(prmtop));
	REQUIRE(std::filesystem::exists(pdb));

	TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
	TopolInfo topol = top_reader->readTopology(prmtop);

	CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
	Configuration conf(coord_reader, pdb, topol);

	SECTION("Number of molecules and bounds") {
		REQUIRE(conf.getNMolec() == topol.num_molecules);
		auto b = conf.getBounds();
		REQUIRE(b.x > 0);
		REQUIRE(b.y > 0);
		REQUIRE(b.z > 0);
	}

	SECTION("Accessors: water and molecule") {
		REQUIRE(conf.getMolec(1).getID() == 1);
		REQUIRE(conf.getWater(1000).isWater());
	}

	SECTION("Per-site interactions for molecule 1000 (V1S–V4S) with v_4S cross-check") {
		int mol_id = 1000;
		auto per_site = conf.getInteractionsPerSite(mol_id);

		std::vector<Real> expected = { -28.3249, -24.1428, -11.7532, -11.3781 };
		checkVectorClose(per_site, expected, 1e-3);

		REQUIRE(conf.v_4S(mol_id) == Approx(expected[3]).margin(1e-3));
	}

	SECTION("findNearby returns valid neighbors and excludes self") {
		auto neigh = conf.findNearby(1000, 5.0);
		REQUIRE(neigh.size > 0);
		for (int i=0;i<neigh.size;i++){
			REQUIRE(neigh.arr[i] != 1000);
		}
		delete[] neigh.arr; // free if ArrInt uses heap for arr
	}

	SECTION("vI and getVList internal consistency") {
		int mol_id = 1000;
		// getVList is ascending; spot-check first and fourth values with hand-calculated ones
		Real* vlist = conf.getVList(mol_id, 5.5);
		REQUIRE(vlist != nullptr);
		REQUIRE(vlist[0] == Approx(-23.0847).margin(1e-3));
		REQUIRE(vlist[3] == Approx(-6.3528).margin(1e-3));
		delete[] vlist;
	}

	SECTION("Per-site (flagged) vs non-flagged and waterOnly consistency") {
		int mol_id = 1000;

		bool ww_flag = false;
		auto per_flagged = conf.getInteractionsPerSite(mol_id, ww_flag);
		auto per_all	 = conf.getInteractionsPerSite(mol_id);
		auto per_water   = conf.getInteractionsPerSite_waterOnly(mol_id);

		REQUIRE(per_flagged.size() == 4);
		REQUIRE(per_all.size() == 4);
		REQUIRE(per_water.size() == 4);

		// flagged list is sorted Ascending (per code), per_all is Ascending => check orders and bounds
		checkIncreasing(per_flagged);

		// Each water-only site's magnitude should not exceed "all" (ions+others add up)
		for (size_t i=0;i<4;i++) {
			REQUIRE(std::abs(per_water[i]) <= std::abs(per_flagged[0]) + std::abs(per_flagged[1]) + std::abs(per_flagged[2]) + std::abs(per_flagged[3]) + 1e-6);
		}

		REQUIRE(ww_flag == false);
	}

	delete top_reader;
	delete coord_reader;
}

// ===================== GROMACS TESTS (expanded) =======================

TEST_CASE("Configuration - GROMACS: extended interaction checks", "[Configuration][GROMACS]") {
	std::string top = "../files/gromacs.top";
	std::string gro = "../files/gromacs.gro";

	REQUIRE(std::filesystem::exists(top));
	REQUIRE(std::filesystem::exists(gro));

	TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS);
	TopolInfo topol = top_reader->readTopology(top);

	CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
	Configuration conf(coord_reader, gro, topol);

	SECTION("Basic info and water check") {
		REQUIRE(conf.getNMolec() == topol.num_molecules);
		REQUIRE(conf.getMolec(1000).isWater());
	}

	SECTION("Per-site interactions for molecule 1000 and v_4S") {
		int mol_id = 1000;
		auto per_site = conf.getInteractionsPerSite(mol_id);
		std::vector<Real> expected = { -30.6020, -23.2096, -22.9963, -22.5206 };
		checkVectorClose(per_site, expected, 1e-3);
		REQUIRE(conf.v_4S(mol_id) == Approx(expected[3]).margin(1e-3));
	}

	SECTION("Neighbors and vI ordering") {
		auto neigh = conf.findNearby(1000, 5.0);
		REQUIRE(neigh.size >= 0);
		for (int i=0;i<neigh.size;i++){
			REQUIRE(neigh.arr[i] != 1000);
		}
		delete[] neigh.arr;

		int mol_id = 1000;
		REQUIRE(conf.vI(mol_id,1) <= conf.vI(mol_id,2));
		REQUIRE(conf.vI(mol_id,2) <= conf.vI(mol_id,3));
		REQUIRE(conf.vI(mol_id,3) <= conf.vI(mol_id,4));
	}

	SECTION("waterOnly vs all sites sanity") {
		int mol_id = 1000;
		auto per_all   = conf.getInteractionsPerSite(mol_id);
		auto per_water = conf.getInteractionsPerSite_waterOnly(mol_id);
		REQUIRE(per_all.size() == 4);
		REQUIRE(per_water.size() == 4);
		REQUIRE(per_all[0] == Approx(-30.6020).margin(1e-3));
		REQUIRE(per_all[1] == Approx(-23.2096).margin(1e-3));
		REQUIRE(per_all[2] == Approx(-22.9963).margin(1e-3));
		REQUIRE(per_all[3] == Approx(-22.5206).margin(1e-3));
		REQUIRE(per_water[0] == Approx(-30.6020).margin(1e-3));
		REQUIRE(per_water[1] == Approx(-22.6554).margin(1e-3));
		REQUIRE(per_water[2] == Approx(-21.4925).margin(1e-3));
		REQUIRE(per_water[3] == Approx(-17.7484).margin(1e-3));
	}

	delete top_reader;
	delete coord_reader;
}

// ===================== LAMMPS TESTS (expanded) =======================

TEST_CASE("Configuration - LAMMPS: extended interaction checks", "[Configuration][LAMMPS]") {
	std::string data = "../files/lammps_system.data";
	std::string traj = "../files/lammps_traj_minimized_9A.xyz";

	REQUIRE(std::filesystem::exists(data));
	REQUIRE(std::filesystem::exists(traj));

	TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS);
	TopolInfo topol = top_reader->readTopology(data);

	CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
	Configuration conf(coord_reader, traj, topol);

	SECTION("Number of molecules") {
		REQUIRE(conf.getNMolec() == topol.num_molecules);
	}

	SECTION("Per-site interactions for molecule 1000") {
		int mol_id = 1000;
		auto per_site = conf.getInteractionsPerSite(mol_id);
		std::vector<Real> expected = { -42.6962, -42.6758, -36.0950, -30.1109 };
		checkVectorClose(per_site, expected, 1e-3);
	}

	SECTION("vI monotonicity and v_4S cross-check when accessible") {
		int mol_id = 1000;
		REQUIRE(conf.vI(mol_id,1) <= conf.vI(mol_id,2));
		REQUIRE(conf.vI(mol_id,2) <= conf.vI(mol_id,3));
	}

	delete top_reader;
	delete coord_reader;
}

// ===================== CONFIGURATION BULK TESTS =======================

TEST_CASE("ConfigurationBulk - classification, potentials matrix, arrays and DJ/metrics", "[ConfigurationBulk]") {
	std::string prmtop = "../files/control_bulk_amber.prmtop";
	std::string pdb	= "../files/control_bulk_amber.pdb";

	REQUIRE(std::filesystem::exists(prmtop));
	REQUIRE(std::filesystem::exists(pdb));

	TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
	TopolInfo topol = top_reader->readTopology(prmtop);

	CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
	ConfigurationBulk conf(coord_reader, pdb, topol);

	const int N = conf.getNMolec();
	REQUIRE(N == topol.num_molecules);

	SECTION("createPotentialMatrix/deletePotentialMatrix shape and initialization") {
		Real** pm = conf.createPotentialMatrix();
		// Triangular allocation: row i has length i (0-based), pm[0] is unused
		for (int i=1; i<N/100; ++i) {
			for (int j=0; j<i; ++j) {
				// All entries must be NOT_CLASSIFIED
				REQUIRE(pm[i][j] == Approx(NOT_CLASSIFIED).margin(0.0));
			}
		}
		conf.deletePotentialMatrix(pm);
	}

	SECTION("classifyMolecules assigns allowed labels and D implication") {
		conf.classifyMolecules(4, -12.0);
		// Allowed: D (0) or T0(1)/T1(2)/T2(3)
		int countD = 0, countT = 0;
		for (int i=topol.num_solutes+1;i<=N/100;i++) 
		{
			const Molecule& m = conf.getMolec(i);
			if (!m.isWater()) continue;
			const Water& w = static_cast<const Water&>(m);
			int c = w.getClassification();
			bool is_something= (c == ConfigurationBulk::CLASSIFICATION_D_MOLECULE) ||
							   (c == ConfigurationBulk::CLASSIFICATION_T0_MOLECULE) ||
							   (c == ConfigurationBulk::CLASSIFICATION_T1_MOLECULE) ||
							   (c == ConfigurationBulk::CLASSIFICATION_T2_MOLECULE);
			REQUIRE(is_something);
			if (c == ConfigurationBulk::CLASSIFICATION_D_MOLECULE) {
				countD++;
				// If labeled D, then v4 must be > threshold
				REQUIRE(conf.vI(i,4) > -12.0);
			} else {
				countT++;
				REQUIRE(conf.vI(i,4) < -12.0);
			}
		}
		// In liquid water there should be at least some D and some T2/T0/T1
		REQUIRE(countD > 0);
		REQUIRE(countT > 0);
	}

	SECTION("classifyMolecules_includePentacoordinated assigns allowed labels and TA has neighbor DX") {
		conf.classifyMolecules_includePentacoordinated(4, -12.0);
		int countDX = 0, countTAorTB = 0;
		for (int i=1;i<=N/100;i++) {
			const Molecule& m = conf.getMolec(i);
			if (!m.isWater()) continue;
			const Water& w = static_cast<const Water&>(m);
			int c = w.getClassification();
			bool is_something= (c == ConfigurationBulk::CLASSIFICATION_D_MOLECULE) ||
							   (c == ConfigurationBulk::CLASSIFICATION_T0_MOLECULE) ||
							   (c == ConfigurationBulk::CLASSIFICATION_T1_MOLECULE) ||
							   (c == ConfigurationBulk::CLASSIFICATION_T2_MOLECULE);
			REQUIRE(is_something);
			if (c == ConfigurationBulk::CLASSIFICATION_D3_MOLECULE ||
				c == ConfigurationBulk::CLASSIFICATION_D5_MOLECULE) countDX++;
			else countTAorTB++;
		}
		REQUIRE(countDX > 0);
		REQUIRE(countTAorTB > 0);
	}

	SECTION("isD / isD3 / isD5 / isDX logical relationships on a sample molecule") {
		int mol_id = 400;
		REQUIRE(conf.isD(mol_id) == true);

		for(int i= 400; i < 600; i++) {
			bool d3 = conf.isD3(i);
			bool d5 = conf.isD5(i);
			bool dx = conf.isDX(i);
			
			// DX <=> D3 or D5
			REQUIRE(dx == (d3 || d5));
		}

	}

	SECTION("v_4S_arr matches per-molecule v_4S for a subset and respects start index") {
		// Build neighbor list container (optional) to exercise that code path
		ToolKit::ArrInt* neigh = nullptr;
		Real* arr = conf.v_4S_arr(/*inic_value*/1, /*i_V*/4, /*R_CUT*/5.0, neigh);
		REQUIRE(arr != nullptr);

		// Check a few positions: index 0 corresponds to molecule 1
		REQUIRE(arr[1000-1] == Approx(conf.v_4S(1000)).margin(1e-6)); // molecule 1000
		REQUIRE(arr[1-1]	== Approx(conf.v_4S(1)).margin(1e-6));	// molecule 1

		delete[] arr;
	}

	SECTION("isDJ returns consistent structure (sum_per_site sorted, flags coherent)") {
		int mol_id = 1000;
		auto dj = conf.classifyDefect(mol_id, /*R_CUT*/5.0, /*V_CUT*/-12.0);

		// sum_per_site must be sorted Ascending (per implementation)
		checkIncreasing(dj.sum_per_site);

		// Coherence checks:
		if (dj.is_DJ) {
			// If DJ, we must have recorded both: a lacking site and a bifurcated site (by construction)
			// Magnitudes are finite
			REQUIRE(std::isfinite(dj.lacking_site_potential));
			REQUIRE(std::isfinite(dj.bifurcated_site_potential));
		}
	}

	SECTION("Tanaka and LSI values are finite and physically plausible") {
		// Pick a water molecule pointer
		Water* w = const_cast<Water*>(&conf.getWater(1000));

		// Tanaka ζ: finite; typically negative in liquid, but enforce only finiteness (no NaNs) to avoid overfitting
		Real zeta = conf.Tanaka(w, /*MAX_D_HB*/3.5, /*MAX_A_HB*/30.0);
		REQUIRE(std::isfinite(zeta));

		// LSI is a variance-like measure of neighbor shell gaps => non-negative
		Real lsi = conf.LSI(1000);
		REQUIRE(std::isfinite(lsi));
		REQUIRE(lsi >= -1e-12); // numerical tolerance around 0
	}

	delete top_reader;
	delete coord_reader;
}

// ===================== CONFIGURATION BULK WITH GROMACS =======================

TEST_CASE("ConfigurationBulk - GROMACS: core methods sanity and v4 hand-check", "[ConfigurationBulk][GROMACS]") {
	std::string top = "../files/gromacs.top";
	std::string gro = "../files/gromacs.gro";

	REQUIRE(std::filesystem::exists(top));
	REQUIRE(std::filesystem::exists(gro));

	TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS);
	TopolInfo topol = top_reader->readTopology(top);

	CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
	ConfigurationBulk conf(coord_reader, gro, topol);

	SECTION("v_4S matches hand value for molecule 1000") {
		int mol_id = 1000;
		std::vector<Real> expected = { -30.6020, -23.2096, -22.9963, -22.5206 };
		auto per_site = conf.getInteractionsPerSite(mol_id);
		checkVectorClose(per_site, expected, 1e-3);
		REQUIRE(conf.v_4S(mol_id) == Approx(expected[3]).margin(1e-3));
	}

	SECTION("Classification runs and assigns values") {
		conf.classifyMolecules();
		int counts[4] = {0,0,0,0};
		for (int i=1;i<=conf.getNMolec()/5;++i) {
			const Molecule& m = conf.getMolec(i);
			if (!m.isWater()) continue;
			const Water& w = static_cast<const Water&>(m);
			int c = w.getClassification();
			REQUIRE(c >= 0);
			REQUIRE(c <= 3);
			counts[c]++;
		}
		REQUIRE(counts[ConfigurationBulk::CLASSIFICATION_D_MOLECULE] >= 0);
	}

	delete top_reader;
	delete coord_reader;
}

// ===================== CONFIGURATION BULK WITH LAMMPS =======================

TEST_CASE("ConfigurationBulk - LAMMPS: core methods sanity and v4 hand-check", "[ConfigurationBulk][LAMMPS]") {
	std::string data = "../files/lammps_system.data";
	std::string traj = "../files/lammps_traj_minimized_9A.xyz";

	REQUIRE(std::filesystem::exists(data));
	REQUIRE(std::filesystem::exists(traj));

	TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS);
	TopolInfo topol = top_reader->readTopology(data);

	CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
	ConfigurationBulk conf(coord_reader, traj, topol);

	SECTION("Per-site interactions and v_4S for molecule 1000") {
		int mol_id = 1000;
		std::vector<Real> expected = { -42.6962, -42.6758, -36.0950, -30.1109 };
		auto per_site = conf.getInteractionsPerSite(mol_id);
		checkVectorClose(per_site, expected, 1e-3);
		REQUIRE(conf.v_4S(mol_id) == Approx(expected[3]).margin(1e-3));
	}

	delete top_reader;
	delete coord_reader;
}

// ===================== DEFECT INFO TESTS =======================

TEST_CASE("Configuration - DefectInfo classification", "[Configuration][DefectInfo]") {
    std::string prmtop = "../files/control_bulk_amber.prmtop";
    std::string pdb    = "../files/control_bulk_amber.pdb";

    REQUIRE(std::filesystem::exists(prmtop));
    REQUIRE(std::filesystem::exists(pdb));

    TopologyReader* top_reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
    TopolInfo topol = top_reader->readTopology(prmtop);

    CoordinateReader* coord_reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    Configuration conf(coord_reader, pdb, topol);

    SECTION("DefectInfo flags are mutually consistent") {
        int mol_id = 1000; // pick arbitrary molecule
        auto info = conf.classifyDefect(mol_id, /*R_CUT*/5.0, /*V_CUT*/-12.0);

        // Flags must be coherent
        if (info.is_DJ) {
            REQUIRE_FALSE(info.is_D3);
            REQUIRE_FALSE(info.is_D5);
            REQUIRE(info.lacking_sites > 0);
            REQUIRE(info.bifurcated_sites > 0);
        } else if (info.is_D3) {
            REQUIRE(info.lacking_sites > 0);
            REQUIRE(info.bifurcated_sites == 0);
        } else if (info.is_D5) {
            REQUIRE(info.lacking_sites == 0);
            REQUIRE(info.bifurcated_sites > 0);
        } else {
            REQUIRE(info.lacking_sites == 0);
            REQUIRE(info.bifurcated_sites == 0);
        }

        // Potentials must be finite
        for (auto v : info.sum_per_site) {
            REQUIRE(std::isfinite(v));
        }

        if (info.bifurcated_sites > 0) {
            REQUIRE(std::isfinite(info.bifurcated_site_potential));
            REQUIRE(std::isfinite(info.bifurcated_individual_potentials.first));
            REQUIRE(std::isfinite(info.bifurcated_individual_potentials.second));
        }

        if (info.lacking_sites > 0) {
            REQUIRE(std::isfinite(info.lacking_site_potential));
        }
    }

    SECTION("Population check: some molecules must be D3, D5, and DJ") {
        int countD3 = 0, countD5 = 0, countDJ = 0;
        for (int i = 1; i <= conf.getNMolec()/5; i++) {
            const Molecule& m = conf.getMolec(i);
            if (!m.isWater()) continue;

            auto info = conf.classifyDefect(i);

            if (info.is_DJ) countDJ++;
            else if (info.is_D3) countD3++;
            else if (info.is_D5) countD5++;
        }
        REQUIRE(countD3 > 0);
        REQUIRE(countD5 > 0);
        REQUIRE(countDJ > 0);
    }

    SECTION("Known defect examples must be correctly identified") {
        int id_D3 = 732;
        int id_D5 = 731;
        int id_DJ = 733;

        auto d3_info = conf.classifyDefect(id_D3);
        REQUIRE(d3_info.is_D3);
        REQUIRE_FALSE(d3_info.is_D5);
        REQUIRE_FALSE(d3_info.is_DJ);

        auto d5_info = conf.classifyDefect(id_D5);
        REQUIRE(d5_info.is_D5);
        REQUIRE_FALSE(d5_info.is_D3);
        REQUIRE_FALSE(d5_info.is_DJ);

        auto dj_info = conf.classifyDefect(id_DJ);
        REQUIRE(dj_info.is_DJ);
        REQUIRE_FALSE(dj_info.is_D3);
        REQUIRE_FALSE(dj_info.is_D5);
    }

    delete top_reader;
    delete coord_reader;
}
