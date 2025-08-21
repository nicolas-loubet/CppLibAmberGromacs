#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

// Mock concrete classes for testing abstract interfaces
class MockCoordinateReader : public CoordinateReader {
public:
	bool readCoordinates(const std::string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
		// Mock: return true if filename is non-empty and topol_info is valid
		if(molecs == nullptr) return false;
		return !filename.empty() && topol_info.num_molecules >= 0;
	}
};

class MockTopologyReader : public TopologyReader {
public:
	TopolInfo readTopology(const std::string& filename) const override {
		TopolInfo info;
		info.num_molecules = 100;
		info.num_solutes = 50;
		info.num_solvents = 50;
		info.total_number_of_atoms = 1000;
		info.number_of_each_different_molecule["Water"] = 50;
		info.number_of_atoms_per_different_molecule["Water"] = 3;
		info.name_type["O"] = "Oxygen";
		info.type_Z["Oxygen"] = 8;
		info.type_LJparam["Oxygen"] = {1.0, 3.5};
		info.special_interaction[{"O", "H"}] = {0.5, 2.5};
		info.default_system_bounds = Vector();
		return info;
	}

	// Helper to access protected periodic_table
	int getPeriodicTableValue(const std::string& element) const {
		return periodic_table.at(element);
	}
	size_t getPeriodicTableSize() const {
		return periodic_table.size();
	}
};

// Helper function to create temporary test files
void createTestFiles(const std::string& directory, const std::vector<std::string>& filenames) {
	std::filesystem::create_directory(directory);
	for (const auto& fname : filenames) {
		std::ofstream file(directory + "/" + fname);
		file << "test content";
		file.close();
	}
}

TEST_CASE("ReaderFactory - CoordinateReader creation", "[ReaderFactory][CoordinateReader]") {
	SECTION("Should create an AmberCoordinateReader for AMBER format") {
		CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
		REQUIRE(dynamic_cast<AmberCoordinateReader*>(reader) != nullptr);
		delete reader;
	}

	SECTION("Should create a GromacsCoordinateReader for GROMACS format") {
		CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
		REQUIRE(dynamic_cast<GromacsCoordinateReader*>(reader) != nullptr);
		delete reader;
	}

	SECTION("Should create a LammpsCoordinateReader for LAMMPS format") {
		CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
		REQUIRE(dynamic_cast<LammpsCoordinateReader*>(reader) != nullptr);
		delete reader;
	}

	SECTION("Should throw an error for unsupported coordinate format") {
		REQUIRE_THROWS_AS(ReaderFactory::createCoordinateReader(static_cast<ReaderFactory::ProgramFormat>(999)), std::runtime_error);
	}
}

TEST_CASE("ReaderFactory - TopologyReader creation", "[ReaderFactory][TopologyReader]") {
	SECTION("Should create an AmberTopologyReader for AMBER format") {
		TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
		REQUIRE(dynamic_cast<AmberTopologyReader*>(reader) != nullptr);
		delete reader;
	}

	SECTION("Should create a GromacsTopologyReader for GROMACS format") {
		TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS);
		REQUIRE(dynamic_cast<GromacsTopologyReader*>(reader) != nullptr);
		delete reader;
	}

	SECTION("Should create a LammpsTopologyReader for LAMMPS format") {
		TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS);
		REQUIRE(dynamic_cast<LammpsTopologyReader*>(reader) != nullptr);
		delete reader;
	}

	SECTION("Should throw an error for unsupported topology format") {
		REQUIRE_THROWS_AS(ReaderFactory::createTopologyReader(static_cast<ReaderFactory::ProgramFormat>(999)), std::runtime_error);
	}
}

TEST_CASE("CoordinateReader - Frame management", "[CoordinateReader]") {
	MockCoordinateReader reader;

	SECTION("Should initialize frame at 0") {
		REQUIRE(reader.getFrame() == 0);
	}

	SECTION("Should set frame correctly") {
		reader.setFrame(5);
		REQUIRE(reader.getFrame() == 5);
		reader.setFrame(0);
		REQUIRE(reader.getFrame() == 0);
	}

	SECTION("Should increment frame correctly") {
		REQUIRE(reader.incFrame() == 1);
		REQUIRE(reader.getFrame() == 1);
		REQUIRE(reader.incFrame() == 2);
		REQUIRE(reader.getFrame() == 2);
	}

	SECTION("Should handle negative frames") {
		reader.setFrame(-1);
		REQUIRE(reader.getFrame() == -1);
		REQUIRE(reader.incFrame() == 0);
	}

	SECTION("Should handle large frames") {
		reader.setFrame(std::numeric_limits<int>::max() - 1);
		REQUIRE(reader.getFrame() == std::numeric_limits<int>::max() - 1);
		REQUIRE(reader.incFrame() == std::numeric_limits<int>::max());
	}
}

TEST_CASE("CoordinateReader - File iterator", "[CoordinateReader]") {
	const std::string test_dir = "test_dir";
	const std::string pattern = "frame*.gro";
	const std::vector<std::string> test_files = {"frame1.gro", "frame2.gro", "frame10.gro"};

	std::filesystem::remove_all(test_dir);

	SECTION("Should find and sort files correctly") {
		createTestFiles(test_dir, test_files);
		auto files = CoordinateReader::getFileIterator(test_dir, pattern);
		REQUIRE(files.size() == 3);
		REQUIRE(files[0].first == 1);
		REQUIRE(files[0].second == "frame1.gro");
		REQUIRE(files[1].first == 2);
		REQUIRE(files[1].second == "frame2.gro");
		REQUIRE(files[2].first == 10);
		REQUIRE(files[2].second == "frame10.gro");
	}

	SECTION("Should throw error if no files are found") {
		std::filesystem::create_directory(test_dir);
		REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, pattern), std::runtime_error);
	}

	SECTION("Should handle patterns with no matches") {
		createTestFiles(test_dir, test_files);
		REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, "invalid*.xyz"), std::runtime_error);
	}

	SECTION("Should ignore non-regular files") {
		createTestFiles(test_dir, test_files);
		std::filesystem::create_directory(test_dir + "/subdir");
		auto files = CoordinateReader::getFileIterator(test_dir, pattern);
		REQUIRE(files.size() == 3); // Subdirectory should be ignored
	}

	SECTION("Should handle nonexistent directory") {
		REQUIRE_THROWS_AS(CoordinateReader::getFileIterator("nonexistent_dir", pattern), std::filesystem::filesystem_error);
	}

	SECTION("Should handle filenames without a number") {
		createTestFiles(test_dir, {"frameX.gro"});
		REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, pattern), std::runtime_error);
		REQUIRE_THROWS_WITH(CoordinateReader::getFileIterator(test_dir, pattern), "No files found matching pattern: " + pattern);
	}

	SECTION("Should throw error for pattern without asterisk") {
		createTestFiles(test_dir, test_files);
		REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, "frame.gro"), std::runtime_error);
		REQUIRE_THROWS_WITH(CoordinateReader::getFileIterator(test_dir, "frame.gro"), "The pattern must contain a *");
	}

	std::filesystem::remove_all(test_dir);
}

TEST_CASE("CoordinateReader - Reading coordinates", "[CoordinateReader]") {
	MockCoordinateReader reader;
	TopolInfo topol_info;
	topol_info.num_molecules = 10;
	std::vector<Molecule*> molecs(10, nullptr);
	Vector bounds;

	SECTION("Should read coordinates with valid file") {
		REQUIRE(reader.readCoordinates("valid.gro", topol_info, molecs.data(), bounds) == true);
	}

	SECTION("Should fail with empty file") {
		REQUIRE(reader.readCoordinates("", topol_info, molecs.data(), bounds) == false);
	}

	SECTION("Should handle invalid topology") {
		TopolInfo invalid_topol;
		invalid_topol.num_molecules = -1;
		REQUIRE(reader.readCoordinates("valid.gro", invalid_topol, molecs.data(), bounds) == false);
	}

	SECTION("Should handle nullptr for molecs") {
		REQUIRE(reader.readCoordinates("valid.gro", topol_info, nullptr, bounds) == false);
	}
}

TEST_CASE("TopologyReader - Reading topology", "[TopologyReader]") {
	MockTopologyReader reader;

	SECTION("Should read topology correctly") {
		TopolInfo info = reader.readTopology("dummy.top");
		REQUIRE(info.num_molecules == 100);
		REQUIRE(info.num_solutes == 50);
		REQUIRE(info.num_solvents == 50);
		REQUIRE(info.total_number_of_atoms == 1000);
		REQUIRE(info.number_of_each_different_molecule["Water"] == 50);
		REQUIRE(info.number_of_atoms_per_different_molecule["Water"] == 3);
		REQUIRE(info.name_type["O"] == "Oxygen");
		REQUIRE(info.type_Z["Oxygen"] == 8);
		REQUIRE(info.type_LJparam["Oxygen"].first == Approx(1.0));
		REQUIRE(info.type_LJparam["Oxygen"].second == Approx(3.5));
		REQUIRE(info.special_interaction[{"O", "H"}].first == Approx(0.5));
		REQUIRE(info.special_interaction[{"O", "H"}].second == Approx(2.5));
	}

	SECTION("Should initialize default_system_bounds correctly") {
		TopolInfo info = reader.readTopology("dummy.top");
		REQUIRE(info.default_system_bounds.x == Approx(0.0));
		REQUIRE(info.default_system_bounds.y == Approx(0.0));
		REQUIRE(info.default_system_bounds.z == Approx(0.0));
	}

	SECTION("Should handle empty file") {
		TopolInfo info = reader.readTopology("");
		REQUIRE(info.num_molecules == 100); // Mock behavior, adjust if real implementation differs
	}
}

TEST_CASE("TopolInfo - Constructor and destructor", "[TopolInfo]") {
	SECTION("Should initialize default values correctly") {
		TopolInfo info;
		REQUIRE(info.num_molecules == 0);
		REQUIRE(info.num_solutes == 0);
		REQUIRE(info.num_solvents == 0);
		REQUIRE(info.total_number_of_atoms == 0);
		REQUIRE(info.number_of_each_different_molecule.empty());
		REQUIRE(info.number_of_atoms_per_different_molecule.empty());
		REQUIRE(info.atom_type_name_charge_mass.empty());
		REQUIRE(info.name_type.empty());
		REQUIRE(info.type_Z.empty());
		REQUIRE(info.type_LJparam.empty());
		REQUIRE(info.special_interaction.empty());
	}

	SECTION("Should destruct without errors") {
		TopolInfo* info = new TopolInfo();
		delete info; // Should not crash
	}
}

TEST_CASE("TopolInfo - Stream output", "[TopolInfo]") {
	MockTopologyReader reader;
	TopolInfo info = reader.readTopology("dummy.top");

	SECTION("Should generate stream output with all fields") {
		std::ostringstream oss;
		oss << info;
		std::string output = oss.str();

		REQUIRE(output.find("Number of molecules: 100") != std::string::npos);
		REQUIRE(output.find("Number of solutes: 50") != std::string::npos);
		REQUIRE(output.find("Number of solvents: 50") != std::string::npos);
		REQUIRE(output.find("Total number of atoms: 1000") != std::string::npos);
		REQUIRE(output.find("Water: 50") != std::string::npos);
		REQUIRE(output.find("Water: 3") != std::string::npos);
		REQUIRE(output.find("O: Oxygen") != std::string::npos);
		REQUIRE(output.find("Oxygen: 8") != std::string::npos);
		REQUIRE(output.find("Oxygen: (e=1, s=3.5)") != std::string::npos);
		REQUIRE(output.find("(O,H): (0.5,2.5)") != std::string::npos);
	}

	SECTION("Should handle empty TopolInfo") {
		TopolInfo empty_info;
		std::ostringstream oss;
		oss << empty_info;
		std::string output = oss.str();
		REQUIRE(output.find("Number of molecules: 0") != std::string::npos);
		REQUIRE(output.find("Number of solutes: 0") != std::string::npos);
		REQUIRE(output.find("Number of solvents: 0") != std::string::npos);
		REQUIRE(output.find("Total number of atoms: 0") != std::string::npos);
	}
}

TEST_CASE("TopologyReader - Periodic table", "[TopologyReader]") {
	MockTopologyReader reader;

	SECTION("Should contain expected elements in periodic table") {
		REQUIRE(reader.getPeriodicTableValue("H") == 1);
		REQUIRE(reader.getPeriodicTableValue("O") == 8);
		REQUIRE(reader.getPeriodicTableValue("C") == 6);
		REQUIRE(reader.getPeriodicTableValue("BA") == 56);
		REQUIRE(reader.getPeriodicTableSize() >= 56); // At least 56 elements
	}

	SECTION("Should throw exception for non-existent element") {
		REQUIRE_THROWS_AS(reader.getPeriodicTableValue("XX"), std::out_of_range);
	}
}
