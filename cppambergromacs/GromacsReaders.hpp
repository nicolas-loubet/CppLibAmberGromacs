#ifndef GROMACS_READERS_HPP
#define GROMACS_READERS_HPP

#include "ReaderFactory.hpp"
#include <string>

class GromacsTopologyReader : public ReaderFactory::TopologyReader {
    public:
        Configuration::TopolInfo readTopology(const std::string& filename) const override {
            // Implementación para leer topología de GROMACS (.top o .tpr)
            return Configuration::TopolInfo();
        }
};

class GromacsCoordinateReader : public ReaderFactory::CoordinateReader {
    private:
        std::string filename;
    public:
        GromacsCoordinateReader(const std::string& file) : filename(file) {}
        bool readCoordinates(Molecule** molecules, int num_molecules, 
                            const Configuration::TopolInfo& topol_info) const override {
            // Implementación para leer .gro o .trr de GROMACS
            return true;
        }
};

#endif