#ifndef READER_FACTORY_HPP
#define READER_FACTORY_HPP

/**
 * Version: June 2025
 * Author: Nicolás Loubet
 */

#include "AmberReaders.hpp"
#include "GromacsReaders.hpp"
// TODO: #include "LammpsReaders.hpp"

class ReaderFactory {
	public:
		enum class ProgramFormat { AMBER, GROMACS, LAMMPS	};

		/**
		 * Factory method to create a CoordinateReader object
		 * @param format The coordinate format (AMBER, GROMACS or LAMMPS)
		 * @param filename The name of the file to read
		 * @return A CoordinateReader object
		 */
		static std::unique_ptr<CoordinateReader> createCoordinateReader(ProgramFormat format, const std::string& filename) {
			switch(format) {
				case ProgramFormat::AMBER:
					return std::make_unique<AmberCoordinateReader>(filename);
				case ProgramFormat::GROMACS:
					return std::make_unique<GromacsCoordinateReader>(filename);
				// TODO: case ProgramFormat::LAMMPS:
				// TODO: 	return std::make_unique<LammpsCoordinateReader>(filename);
				default:
					throw std::runtime_error("Unsupported coordinate format");
			}
		}
	
		/**
		 * Factory method to create a TopologyInfo struct
		 * @param format The topology format (AMBER, GROMACS or LAMMPS)
		 * @param filename The name of the file to read
		 * @return A TopolInfo struct
		 */
		static std::unique_ptr<TopologyReader> createTopologyReader(ProgramFormat format) {
			switch(format) {
				case ProgramFormat::AMBER:
					return std::make_unique<AmberTopologyReader>();
				case ProgramFormat::GROMACS:
					return std::make_unique<GromacsTopologyReader>();
				// TODO: case ProgramFormat::LAMMPS:
				// TODO: 	return std::make_unique<LammpsTopologyReader>();
				default:
					throw std::runtime_error("Unsupported topology format");
			}
		}
};

#endif
