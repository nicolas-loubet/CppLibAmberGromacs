#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

/**
 * Version: June 2025
 * Author: Ezequiel Cuenca
 */

#include "ReaderFactory.hpp"
#include "amber_parser.hpp"
#include <string>
#include <fstream>
#include <sstream>

class AmberTopologyReader : public TopologyReader {

    public:
        TopolInfo readTopology(const string& filename) const override {
            // Implementación para leer prmtop de AMBER-EZEQUIEL

            TopolInfo topology= TopolInfo();

            ifstream file(filename);
            map<string, int> positions_of_flags=flag_position(file);
            map<string,int> dict_pointers = read_pointers(file,positions_of_flags["POINTERS"]);

            vector<tuple<string,float>> map_atoms = read_charge(file, dict_pointers, read_atom_name(file, dict_pointers,positions_of_flags["ATOM_NAME"]),positions_of_flags["CHARGE"]);

            vector<int> number_solute_solvent =read_solvent_pointers(file, dict_pointers,positions_of_flags["SOLVENT_POINTERS"]);
            map<string,tuple<int,int>> atoms_per_diff_molecule=read_atoms_per_different_molecule(file, positions_of_flags);
            
            vector<int> atoms_per_molecule = read_atoms_per_molecule(file, dict_pointers,positions_of_flags["SOLVENT_POINTERS"]); 
            vector<int> atom_type_index=read_ati(file,dict_pointers,positions_of_flags["ATOM_TYPE_INDEX"]); 
            vector<float> mass=read_mass(file,dict_pointers,positions_of_flags["MASS"]);
            map<tuple<int,int>,tuple<float,float>> lj_coefficient=read_lj(file,dict_pointers,positions_of_flags["LENNARD_JONES_ACOEF"]);

            topology.num_molecules=number_solute_solvent[1];
			topology.num_solutes=number_solute_solvent[0];
			topology.num_solvents==number_solute_solvent[1]-number_solute_solvent[0];
            
            for (const auto& par : atoms_per_diff_molecule) {
                 //cout << "Clave: " << par.first << ", Atomos: " << get<0>(par.second) << " numero de moleculas: " << get<1>(par.second)<<endl;
                 topology.number_of_atoms_per_different_molecule[par.first]=get<0>(par.second);
                 topology.number_of_each_different_molecule[par.first]=get<1>(par.second);
            }
			//map<string,int> number_of_each_different_molecule;
			//topology.number_of_atoms_per_different_molecule;
            
            map<int, tuple<string, string, float, float>> molecule_atoms;
            
            for (size_t i = 0; i < dict_pointers["NATOM"]; i+=1) {
                stringstream ss;
                ss << atom_type_index[i];
                molecule_atoms[i] = make_tuple(ss.str(), get<0>(map_atoms[i]), get<1>(map_atoms[i]), mass[i]);
            }

            topology.atom_type_name_charge_mass.push_back(molecule_atoms);

			map<string,string> name_type;
            map<string,pair<float,float>> type_LJparam; //0=epsilon 1=sigma
            map<pair<string,string>,pair<float,float>> special_interaction;//keep ij and ji
            //return Configuration::TopolInfo();
            return topology;
        }
};

class AmberCoordinateReader : public CoordinateReader {
    private:
        /**
         * Given an atom name, return its atomic number
         * @param atom_name The name of the atom
         * @return The atomic number
         */
        static int getAtomicNumber(const string& atom_name, const string& atom_type) {
            static const map<string,int> atomic_numbers= {
                {"H",1},{"B",5},{"C",6},{"N",7},{"O",8},
                {"F",9},{"Na",11},{"Mg",12},{"Al",13},
                {"Si",14},{"P",15},{"S",16},{"Cl",17},
                {"K",19},{"Br",35},{"I",53}
            };
            if(atom_name.length() >= 2) {
                string symbol= atom_name.substr(0,2);
                auto it= atomic_numbers.find(symbol);
                if(it != atomic_numbers.end())
                    return it->second;
            }

            string symbol= atom_name.substr(0,1);
            auto it= atomic_numbers.find(symbol);
            if(it != atomic_numbers.end())
                return it->second;

            if(!atom_type.empty()) {
                auto type_it= atomic_numbers.find(atom_type);
                if(type_it != atomic_numbers.end())
                    return type_it->second;
            }

            return 0; //Default value
        }
    public:
        AmberCoordinateReader()= default;

        /**
         * Reads the coordinates file .pdb
         * @param filename The name of the coordinates file (single frame)
         * @param topol_info The topology information (use AmberTopologyReader)
         * @param molecs An empty array of molecule pointers
         * @return True if the coordinates were read successfully
         */
        bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs) const override {
            ifstream f(filename);
            if(!f.is_open()) {
                cout << "Failed to open file " << filename << endl;
                return false;
            }

            int total_atoms= 0;
            for(const auto& pair: topol_info.number_of_each_different_molecule) {
                const string& mol_name= pair.first;
                int num_molecules= pair.second;
                int atoms_per_molecule= topol_info.number_of_atoms_per_different_molecule.at(mol_name);
                total_atoms+= num_molecules * atoms_per_molecule;
            }

            string line;
            vector<Vector> coords;
            vector<string> atom_names;

            while(getline(f,line)) {
                if(line.rfind("ATOM  ",0) == 0 || line.rfind("HETATM",0) == 0) {
                    stringstream ss(line);
                    string record, atom_name, res_name;
                    int atom_id, res_id;
                    float x, y, z;

                    ss >> record >> atom_id >> atom_name >> res_name >> res_id >> x >> y >> z;
                    if(ss.fail()) {
                        f.close();
                        return false;
                    }

                    coords.emplace_back(Vector(x,y,z));
                    atom_names.push_back(atom_name);
                }
            }
            f.close();

            int atom_idx= 0;
            int molec_idx= 0;
            for(const auto& pair: topol_info.number_of_each_different_molecule) {
                const string& mol_name= pair.first;
                int num_molecules= pair.second;
                int num_atoms_per_molecule= topol_info.number_of_atoms_per_different_molecule.at(mol_name);

                for(int i= 0; i < num_molecules; i++,molec_idx++) {
                    Atom* atoms= new Atom[num_atoms_per_molecule];

                    for(int j= 0; j < num_atoms_per_molecule; j++,atom_idx++) {
                        const auto& atom_data= topol_info.atom_type_name_charge_mass[atom_idx];
                        const auto& [type,name,charge,mass]= atom_data.at(atom_idx);
                        const auto& [epsilon,sigma]= topol_info.type_LJparam.at(type);

                        int Z= getAtomicNumber(atom_names[atom_idx],type);
                        atoms[j]= Atom(coords[atom_idx], atom_idx+1, mass, charge, epsilon, sigma, Z);
                    }

                    molecs[molec_idx]= new Molecule(molec_idx+1, atoms, num_atoms_per_molecule);
                }
            }

            return true;
        }
};

#endif
