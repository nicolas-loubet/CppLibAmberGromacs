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
#include <regex>
#include <filesystem>
#include "Molecules/Water.hpp"
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
 * Abstract class for coordinate readers
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
     * @param bounds The system bounds to be modified (Vector object still not created)
     * @return True if the coordinates were read successfully
     */
    virtual bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const= 0;

    static vector<pair<int,string>> getFileIterator(const string& directory, const string& pattern) {
        string regex_pattern= regex_replace(pattern, regex("\\*"), "(\\d+)");
        regex file_regex(regex_pattern);
        vector<pair<int,string>> list_files;

        for(const auto& entry: filesystem::directory_iterator(directory)) {
            if(!entry.is_regular_file()) continue;
            string filename= entry.path().filename().string();
            smatch match;
            if(regex_match(filename, match, file_regex))
                list_files.emplace_back(stoi(match[1].str()),filename);
        }

        if(list_files.empty())
            throw runtime_error("No files found matching pattern: " + pattern);

        sort(list_files.begin(), list_files.end(),
             [](const auto& a,const auto& b) { return a.first < b.first; });

        return list_files;
    }
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

/**
 * Output stream operator for TopolInfo
 */
std::ostream& operator<<(std::ostream& os, const TopolInfo& info) {
    os << "TopolInfo:\n";
    os << "  Number of molecules: " << info.num_molecules << "\n";
    os << "  Number of solutes: " << info.num_solutes << "\n";
    os << "  Number of solvents: " << info.num_solvents << "\n";
    os << "  Total number of atoms: " << info.total_number_of_atoms << "\n";

    os << "  Number of each different molecule:\n";
    for(const auto& [molecule, count]: info.number_of_each_different_molecule)
        os << "    " << molecule << ": " << count << "\n";

    os << "  Number of atoms per different molecule:\n";
    for(const auto& [molecule, count]: info.number_of_atoms_per_different_molecule)
        os << "    " << molecule << ": " << count << "\n";

    os << "  Atom type name charge mass:\n";
    for(size_t i = 0; i < info.atom_type_name_charge_mass.size(); i++) {
        os << "    Map " << i << ":\n";
        for(const auto& [index, data]: info.atom_type_name_charge_mass[i]) {
            os << "      " << index << ": (" << std::get<0>(data) << ", "
               << std::get<1>(data) << ", "  << std::get<2>(data) << ", " << std::get<3>(data) << ")\n";
        }
    }

    os << "  Name type:\n";
    for(const auto& [name, type]: info.name_type)
        os << "    " << name << ": " << type << "\n";

    os << "  Type Z:\n";
    for(const auto& [type, z]: info.type_Z)
        os << "    " << type << ": " << z << "\n";

    os << "  Type LJ parameters:\n";
    for (const auto& [type, params]: info.type_LJparam)
        os << "    " << type << ": (e=" << params.first << ", s=" << params.second << ")\n";

    os << "  Special interactions:\n";
    for(const auto& [pair, params]: info.special_interaction)
        os << "    (" << pair.first << "," << pair.second << "): (" << params.first << "," << params.second << ")\n";

    return os;
}

#endif
