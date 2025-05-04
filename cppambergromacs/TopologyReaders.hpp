#ifndef TOPOLOGYREADER_HPP
#define TOPOLOGYREADER_HPP

/**
 * Version: May 2025
 * Author: Nicolás Loubet
 */

#include "Configuration.hpp"

namespace TopolReader {
	static Configuration::TopolInfo readTopolInfoOfAmber(std::string filename) {
		// Implementación para leer archivo de topología de AMBER
		return Configuration::TopolInfo();
	}

	static Configuration::TopolInfo readTopolInfoOfGromacs(std::string filename) {
		// Implementación para leer archivo de topología de GROMACS
		return Configuration::TopolInfo();
	}
}

#endif