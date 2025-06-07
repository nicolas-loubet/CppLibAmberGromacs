#ifndef READER_INTERFACES_HPP
#define READER_INTERFACES_HPP

/**
 * Version: June 2025
 * Author: Nicolás Loubet
 */

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <tuple>
#include "Molecules/Molecule.hpp"

/**
 * This struct is an adapter to read the topology file
 */
struct TopolInfo {
    int num_molecules;
    int num_solutes;
    int num_solvents;
    std::map<std::string,int> number_of_each_different_molecule;
    std::map<std::string,int> number_of_atoms_per_different_molecule;
    std::vector<std::map<int,std::tuple<std::string,std::string,float,float>>> atom_type_name_charge_mass;
    std::map<std::string,std::string> name_type;
    std::map<std::string,std::pair<float,float>> type_LJparam; //0=epsilon 1=sigma
    std::map<std::pair<std::string,std::string>,std::pair<float,float>> special_interaction;//keep ij and ji

    TopolInfo(): num_molecules(0), num_solutes(0), num_solvents(0) {}
    ~TopolInfo()= default;
};

/**
 * Interface for coordinate readers
 */
class CoordinateReader {
public:
	virtual ~CoordinateReader()= default;
	virtual bool readCoordinates(Molecule** molecules, int num_molecules, const TopolInfo& topol_info) const= 0;
};

/**
 * Interface for topology readers
 */
class TopologyReader {
public:
	virtual ~TopologyReader()= default;
	virtual TopolInfo readTopology(const std::string& filename) const= 0;
};

#endif
