#ifndef READER_FACTORY_HPP
#define READER_FACTORY_HPP

#include "Configuration.hpp"
#include "CoordinateReaders.hpp"
#include "TopologyReaders.hpp"

class ReaderFactory {
    public:
        static constexpr int AMBER= 0;
        static constexpr int GROMACS= 1;
        static constexpr int LAMMPS= 2;

        static Configuration::CoordinateReader* createCoordinateReader(const int format, const std::string& filename) {
            if(format == AMBER) return new AmberCoordinateReader(filename);
            if(format == GROMACS) return new GromacsCoordinateReader(filename);
            //if(format == LAMMPS) return new LammpsCoordinateReader(filename);
            throw std::runtime_error("Unsupported coordinate format");
        }
    
        static Configuration::TopolInfo createTopologyReader(const int format, const std::string& filename) {
            if(format == AMBER) return TopolReader::readTopolInfoOfAmber(filename);
            if(format == GROMACS) return TopolReader::readTopolInfoOfGromacs(filename);
            //if(format == LAMMPS) return TopolReader::readTopolInfoOfLammps(filename);
            throw std::runtime_error("Unsupported topology format");
        }
};

#endif
