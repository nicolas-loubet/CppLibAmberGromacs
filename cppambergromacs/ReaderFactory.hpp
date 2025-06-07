#ifndef READER_FACTORY_HPP
#define READER_FACTORY_HPP

/**
 * Version: May 2025
 * Author: Nicolás Loubet
 */

#include "Configuration.hpp"
#include "AmberReaders.hpp"
#include "GromacsReaders.hpp"
//#include "LammpsReaders.hpp"

class ReaderFactory {
	public:
		enum class ProgramFormat { AMBER, GROMACS, LAMMPS	};

		/**
		 * Abstract class to create a CoordinateReader object (adapter)
		 */
		class CoordinateReader {
			public:
				virtual ~CoordinateReader() = default;
				virtual bool readCoordinates(Molecule** molecules, int num_molecules, const Configuration::TopolInfo& topol_info) const= 0;
		};

		/**
		 * Abstract class to create a TopologyReader object (adapter)
		 */
		class TopologyReader {
			public:
				virtual ~TopologyReader() = default;
				virtual Configuration::TopolInfo readTopology(const std::string& filename) const= 0;
		};

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
				//case ProgramFormat::LAMMPS:
				//	return std::make_unique<LammpsCoordinateReader>(filename);
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
				//case ProgramFormat::LAMMPS:
				//	return std::make_unique<LammpsTopologyReader>();
				default:
					throw std::runtime_error("Unsupported topology format");
			}
		}
};

#endif
