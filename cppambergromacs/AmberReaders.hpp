#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

#include "ReaderFactory.hpp"
#include <string>

class AmberTopologyReader : public ReaderFactory::TopologyReader {
    public:
        Configuration::TopolInfo readTopology(const std::string& filename) const override {
            // Implementación para leer prmtop de AMBER
            return Configuration::TopolInfo();
        }
};

class AmberCoordinateReader : public ReaderFactory::CoordinateReader {
    private:
        std::string filename;
    public:
        AmberCoordinateReader(const std::string& file): filename(file) {}
        bool readCoordinates(Molecule** molecules, int num_molecules, const Configuration::TopolInfo& topol_info) const override {
            // Implementación para leer mdcrd de AMBER
            return true;
        }
};

#endif