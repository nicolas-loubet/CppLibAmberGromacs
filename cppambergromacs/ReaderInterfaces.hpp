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
#include "General/ToolKit.hpp"

using namespace std;

/**
 * This struct is an adapter to read the topology file
 */
struct TopolInfo {
    int num_molecules;
    int num_solutes;
    int num_solvents;
    int total_number_of_atoms;
    map<string,int> number_of_each_different_molecule;
    map<string,int> number_of_atoms_per_different_molecule;
    vector<map<int,tuple<string,string,float,float>>> atom_type_name_charge_mass;
    map<string,string> name_type;
    map<string,int> type_Z;
    map<string,pair<float,float>> type_LJparam; //0=epsilon 1=sigma
    map<pair<string,string>,pair<float,float>> special_interaction;//keep ij and ji

    TopolInfo(): num_molecules(0), num_solutes(0), num_solvents(0) {}
    ~TopolInfo()= default;
};

/**
 * Interface for coordinate readers
 */
class CoordinateReader {
public:
    CoordinateReader()= default;
	virtual ~CoordinateReader()= default;

    /**
     * Reads the coordinates file
     * @param filename The name of the coordinates file
     * @param topol_info The topology information
     * @param molecs An empty array of molecule pointers
     * @return True if the coordinates were read successfully
     */
	virtual bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs) const= 0;
};

/**
 * Interface for topology readers
 */
class TopologyReader {
public:
    TopologyReader()= default;
	virtual ~TopologyReader()= default;
	virtual TopolInfo readTopology(const string& filename) const= 0;
};

#endif
