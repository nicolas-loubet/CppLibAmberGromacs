#ifndef READER_FACTORY_HPP
#define READER_FACTORY_HPP

/**
 * Version: June 2025
 * Author: Nicolás Loubet
 */

#include "AmberReaders.hpp"
#include "GromacsReaders.hpp"
#ifndef USE_VECTOR_TOPOLOGY
#include "LammpsReaders.hpp"
#endif

class ReaderFactory {
	public:
		enum class ProgramFormat { AMBER, GROMACS, LAMMPS };

		/**
		 * Factory method to create a CoordinateReader object
		 * @param format The coordinate format (AMBER, GROMACS or LAMMPS)
		 * @param filename The name of the file to read
		 * @return A CoordinateReader object
		 */
		static CoordinateReader* createCoordinateReader(ProgramFormat format) {
			switch(format) {
				case ProgramFormat::AMBER:
					return new AmberCoordinateReader();
				case ProgramFormat::GROMACS:
					return new GromacsCoordinateReader();
#ifndef USE_VECTOR_TOPOLOGY
				case ProgramFormat::LAMMPS:
				 	return new LammpsCoordinateReader();
#endif
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
		static TopologyReader* createTopologyReader(ProgramFormat format) {
			switch(format) {
				case ProgramFormat::AMBER:
					return new AmberTopologyReader();
				case ProgramFormat::GROMACS:
					return new GromacsTopologyReader();
#ifndef USE_VECTOR_TOPOLOGY
				case ProgramFormat::LAMMPS:
				 	return new LammpsTopologyReader();
#endif
				default:
					throw std::runtime_error("Unsupported topology format");
			}
		}
};

#endif // READER_FACTORY_HPP
