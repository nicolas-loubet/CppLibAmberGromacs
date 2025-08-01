#define CATCH_CONFIG_MAIN // Esto define el main de Catch2
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

TEST_CASE("ReaderFactory - Creación de CoordinateReader", "[ReaderFactory][CoordinateReader]") {
    SECTION("Debería crear un AmberCoordinateReader para el formato AMBER") {
        CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
        REQUIRE(dynamic_cast<AmberCoordinateReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería crear un GromacsCoordinateReader para el formato GROMACS") {
        CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
        REQUIRE(dynamic_cast<GromacsCoordinateReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería crear un LammpsCoordinateReader para el formato LAMMPS") {
        CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
        REQUIRE(dynamic_cast<LammpsCoordinateReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería lanzar un error para un formato de coordenada no soportado") {
        REQUIRE_THROWS_AS(ReaderFactory::createCoordinateReader(static_cast<ReaderFactory::ProgramFormat>(999)), std::runtime_error);
    }
}

TEST_CASE("ReaderFactory - Creación de TopologyReader", "[ReaderFactory][TopologyReader]") {
    SECTION("Debería crear un AmberTopologyReader para el formato AMBER") {
        TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
        REQUIRE(dynamic_cast<AmberTopologyReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería crear un GromacsTopologyReader para el formato GROMACS") {
        TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS);
        REQUIRE(dynamic_cast<GromacsTopologyReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería crear un LammpsTopologyReader para el formato LAMMPS") {
        TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS);
        REQUIRE(dynamic_cast<LammpsTopologyReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería lanzar un error para un formato de topología no soportado") {
        REQUIRE_THROWS_AS(ReaderFactory::createTopologyReader(static_cast<ReaderFactory::ProgramFormat>(999)), std::runtime_error);
    }
}

TEST_CASE("CoordinateReader - Gestión de frames", "[CoordinateReader]") {
    MockCoordinateReader reader;

    SECTION("Debería inicializar el frame en 0") {
        REQUIRE(reader.getFrame() == 0);
    }

    SECTION("Debería establecer el frame correctamente") {
        reader.setFrame(5);
        REQUIRE(reader.getFrame() == 5);
        reader.setFrame(0);
        REQUIRE(reader.getFrame() == 0);
    }

    SECTION("Debería incrementar el frame correctamente") {
        REQUIRE(reader.incFrame() == 1);
        REQUIRE(reader.getFrame() == 1);
        REQUIRE(reader.incFrame() == 2);
        REQUIRE(reader.getFrame() == 2);
    }

    SECTION("Debería manejar frames negativos") {
        reader.setFrame(-1);
        REQUIRE(reader.getFrame() == -1);
        REQUIRE(reader.incFrame() == 0);
    }

    SECTION("Debería manejar frames grandes") {
        reader.setFrame(std::numeric_limits<int>::max() - 1);
        REQUIRE(reader.getFrame() == std::numeric_limits<int>::max() - 1);
        REQUIRE(reader.incFrame() == std::numeric_limits<int>::max());
    }
}

TEST_CASE("CoordinateReader - Iterador de archivos", "[CoordinateReader]") {
    const std::string test_dir = "test_dir";
    const std::string pattern = "frame*.gro";
    const std::vector<std::string> test_files = {"frame1.gro", "frame2.gro", "frame10.gro"};

    std::filesystem::remove_all(test_dir);

    SECTION("Debería encontrar y ordenar archivos correctamente") {
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

    SECTION("Debería lanzar un error si no se encuentran archivos") {
        std::filesystem::create_directory(test_dir);
        REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, pattern), std::runtime_error);
    }

    SECTION("Debería manejar patrones sin coincidencias") {
        createTestFiles(test_dir, test_files);
        REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, "invalid*.xyz"), std::runtime_error);
    }

    SECTION("Debería ignorar archivos no regulares") {
        createTestFiles(test_dir, test_files);
        std::filesystem::create_directory(test_dir + "/subdir");
        auto files = CoordinateReader::getFileIterator(test_dir, pattern);
        REQUIRE(files.size() == 3); // Subdirectory should be ignored
    }

    SECTION("Debería manejar directorio inexistente") {
        REQUIRE_THROWS_AS(CoordinateReader::getFileIterator("nonexistent_dir", pattern), std::filesystem::filesystem_error);
    }

    SECTION("Debería manejar nombres de archivo sin número") {
        createTestFiles(test_dir, {"frameX.gro"});
        REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, pattern), std::runtime_error);
        REQUIRE_THROWS_WITH(CoordinateReader::getFileIterator(test_dir, pattern), "No files found matching pattern: " + pattern);
    }

    SECTION("Debería lanzar un error para patrón sin asterisco") {
        createTestFiles(test_dir, test_files);
        REQUIRE_THROWS_AS(CoordinateReader::getFileIterator(test_dir, "frame.gro"), std::runtime_error);
        REQUIRE_THROWS_WITH(CoordinateReader::getFileIterator(test_dir, "frame.gro"), "The pattern must contain a *");
    }

    std::filesystem::remove_all(test_dir);
}

TEST_CASE("CoordinateReader - Lectura de coordenadas", "[CoordinateReader]") {
    MockCoordinateReader reader;
    TopolInfo topol_info;
    topol_info.num_molecules = 10;
    std::vector<Molecule*> molecs(10, nullptr);
    Vector bounds;

    SECTION("Debería leer coordenadas con archivo válido") {
        REQUIRE(reader.readCoordinates("valid.gro", topol_info, molecs.data(), bounds) == true);
    }

    SECTION("Debería fallar con archivo vacío") {
        REQUIRE(reader.readCoordinates("", topol_info, molecs.data(), bounds) == false);
    }

    SECTION("Debería manejar topología inválida") {
        TopolInfo invalid_topol;
        invalid_topol.num_molecules = -1;
        REQUIRE(reader.readCoordinates("valid.gro", invalid_topol, molecs.data(), bounds) == false);
    }

    SECTION("Debería manejar nullptr para molecs") {
        REQUIRE(reader.readCoordinates("valid.gro", topol_info, nullptr, bounds) == false);
    }
}

TEST_CASE("TopologyReader - Lectura de topología", "[TopologyReader]") {
    MockTopologyReader reader;

    SECTION("Debería leer topología correctamente") {
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

    SECTION("Debería inicializar default_system_bounds correctamente") {
        TopolInfo info = reader.readTopology("dummy.top");
        REQUIRE(info.default_system_bounds.x == Approx(0.0));
        REQUIRE(info.default_system_bounds.y == Approx(0.0));
        REQUIRE(info.default_system_bounds.z == Approx(0.0));
    }

    SECTION("Debería manejar archivo vacío") {
        TopolInfo info = reader.readTopology("");
        REQUIRE(info.num_molecules == 100); // Mock behavior, adjust if real implementation differs
    }
}

TEST_CASE("TopolInfo - Constructor y destructor", "[TopolInfo]") {
    SECTION("Debería inicializar valores por defecto correctamente") {
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

    SECTION("Debería destruir sin errores") {
        TopolInfo* info = new TopolInfo();
        delete info; // Should not crash
    }
}

TEST_CASE("TopolInfo - Salida a stream", "[TopolInfo]") {
    MockTopologyReader reader;
    TopolInfo info = reader.readTopology("dummy.top");

    SECTION("Debería generar salida de stream con todos los campos") {
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

    SECTION("Debería manejar TopolInfo vacío") {
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

TEST_CASE("TopologyReader - Tabla periódica", "[TopologyReader]") {
    MockTopologyReader reader;

    SECTION("Debería contener elementos esperados en la tabla periódica") {
        REQUIRE(reader.getPeriodicTableValue("H") == 1);
        REQUIRE(reader.getPeriodicTableValue("O") == 8);
        REQUIRE(reader.getPeriodicTableValue("C") == 6);
        REQUIRE(reader.getPeriodicTableValue("BA") == 56);
        REQUIRE(reader.getPeriodicTableSize() >= 56); // At least 56 elements
    }

    SECTION("Debería lanzar excepción para elemento no existente") {
        REQUIRE_THROWS_AS(reader.getPeriodicTableValue("XX"), std::out_of_range);
    }
}

#include "MoleculeTests.cpp"
