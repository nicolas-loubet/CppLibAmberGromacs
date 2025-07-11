#ifndef CONFIGURATIONLIPID_HPP
#define CONFIGURATIONLIPID_HPP

/**
 * Version: July 2025
 * Author: Ezequiel Cuenca
 */

#include "Configuration.hpp"

/**
 * This class creates a Configuration object, with an array of molecules (specific for lipid systems)
 */
class ConfigurationLipid : public Configuration {
	public:
        ConfigurationLipid(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info) :
            Configuration(coord_reader, filename, topol_info) {}

        array<vector<float>,2> orderParameter() {
            array<vector<float>,2> op;
            //Ejemplo de insertar al comienzo un valor (10.5), hace el "kick" automático
            op[0].insert(op[0].begin(), 10.5f);
            //TODO
            throw runtime_error("Method not implemented");
        }
		
};

#endif
