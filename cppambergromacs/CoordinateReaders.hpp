#ifndef COORDINATEREADER_HPP
#define COORDINATEREADER_HPP

#include "Configuration.hpp"

class AmberCoordinateReader : public Configuration::CoordinateReader {
    private:
        std::string filename;
    public:
        AmberCoordinateReader(const std::string& file) : filename(file) {}
        bool readCoordinates(const Molecule** molecules, const int num_molecules, const Configuration::TopolInfo& topol_info) const override {
            // Implementación para leer archivo de coordenadas de AMBER
            return true;
        }
};
    
class GromacsCoordinateReader : public Configuration::CoordinateReader {
    private:
        std::string filename;
    public:
        GromacsCoordinateReader(const std::string& file) : filename(file) {}
        bool readCoordinates(const Molecule** molecules, const int num_molecules, const Configuration::TopolInfo& topol_info) const override {
            // Implementación para leer archivo .gro o .trr de GROMACS
            return true;
        }
};

#endif
