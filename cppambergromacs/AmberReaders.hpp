#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

/**
 * Version: June 2025
 * Author: Ezequiel Cuenca
 */

#include "ReaderInterfaces.hpp"
#include <string>
#include <fstream>
#include <sstream>

class AmberTopologyReader : public TopologyReader {
    private:
        /**
         * Given a file, find the position of the flag in the file.
         * The returned map contains the flag name as the key and the position as the value.
         * The position is the byte position of the first character after the flag.
         * @param file The file to search
         * @return A map with the position of the flag.
         */
        static map<string,int> flag_position(ifstream &file)
        {
            string line= "";
            map<string,int> flag;
            while(getline(file, line)) {
                if(line.find("%FLAG") != string::npos) {
                    line= ToolKit::strip(line.substr(5, 80));
                    flag[line]= file.tellg();
                    flag[line]-= 81;
                }
            }
            return(flag);
        }

    public:
        AmberTopologyReader()= default;
        TopolInfo readTopology(const std::string& filename) const override {
            // Implementación para leer prmtop de AMBER-EZEQUIEL
            return TopolInfo();
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
