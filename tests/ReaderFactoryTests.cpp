#define CATCH_CONFIG_MAIN // Esto define el main de Catch2
#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"

TEST_CASE("ReaderFactory - Creación de CoordinateReader", "[ReaderFactory][CoordinateReader]") {
    SECTION("Debería crear un GromacsCoordinateReader para el formato GROMACS") {
        CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
        // Usamos dynamic_cast para verificar que el puntero devuelto es del tipo esperado.
        REQUIRE(dynamic_cast<GromacsCoordinateReader*>(reader) != nullptr);
        delete reader; // Importante liberar la memoria
    }

    SECTION("Debería crear un LammpsCoordinateReader para el formato LAMMPS") {
        CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
        REQUIRE(dynamic_cast<LammpsCoordinateReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería crear un AmberCoordinateReader para el formato AMBER") {
        CoordinateReader* reader = ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
        REQUIRE(dynamic_cast<AmberCoordinateReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería lanzar un error para un formato de coordenada no soportado") {
        // Asumiendo que existe un valor que no está en el enum, o si quisiéramos probar el default case [12]
        // Para este ejemplo, lanzaré una excepción si no es uno de los 3 formatos.
        // Podrías necesitar un valor 'desconocido' en el enum si quieres probar el default.
        // O simplemente verificas que el factory devuelve nullptr o lanza una excepción por defecto si no lo encuentra.
        // Las fuentes indican que lanza std::runtime_error para formatos no soportados [12, 13].
        REQUIRE_THROWS_AS(ReaderFactory::createCoordinateReader(static_cast<ReaderFactory::ProgramFormat>(999)), std::runtime_error);
    }
}

TEST_CASE("ReaderFactory - Creación de TopologyReader", "[ReaderFactory][TopologyReader]") {
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

    SECTION("Debería crear un AmberTopologyReader para el formato AMBER") {
        TopologyReader* reader = ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER);
        REQUIRE(dynamic_cast<AmberTopologyReader*>(reader) != nullptr);
        delete reader;
    }

    SECTION("Debería lanzar un error para un formato de topología no soportado") {
        REQUIRE_THROWS_AS(ReaderFactory::createTopologyReader(static_cast<ReaderFactory::ProgramFormat>(999)), std::runtime_error);
    }
}
