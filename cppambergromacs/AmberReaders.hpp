#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

/**
 * Version: June 2025
 * Author: Ezequiel Cuenca
 */

#include "ReaderInterfaces.hpp"
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
            vector<string> atom_names = read_atom_name(file, dict_pointers,positions_of_flags["ATOM_NAME"]);
            vector<tuple<string,float>> map_atoms = read_charge(file, dict_pointers, atom_names ,positions_of_flags["CHARGE"]);
            vector<int> number_solute_solvent =read_solvent_pointers(file, dict_pointers,positions_of_flags["SOLVENT_POINTERS"]);
            vector<string> ati_to_amber_type=atom_type_index_to_amber_type(file,  dict_pointers , positions_of_flags);        
            vector<int> atomic_number=read_atomic_number(file,dict_pointers,positions_of_flags["ATOMIC_NUMBER"]); 
            map<string,tuple<int,int>> atoms_per_diff_molecule=read_atoms_per_different_molecule(file, positions_of_flags);
            vector<int> atoms_per_molecule = read_atoms_per_molecule(file, dict_pointers,positions_of_flags["SOLVENT_POINTERS"]); 
            vector<int> atom_type_index=read_ati(file,dict_pointers,positions_of_flags["ATOM_TYPE_INDEX"]); 
            vector<float> mass=read_mass(file,dict_pointers,positions_of_flags["MASS"]);
            map<pair<string,string>,pair<float,float>> lj_coefficient=read_lj(file,dict_pointers,positions_of_flags["LENNARD_JONES_ACOEF"],  ati_to_amber_type);
            map<string,pair<float,float>> lj_diagonal= read_lj_diagonal(lj_coefficient);
            map<string,int> type_atomic_z = type_atomic_number(dict_pointers, atom_type_index,atomic_number,ati_to_amber_type);
            map<string,string> name_types = read_name_type(dict_pointers, atom_names,atom_type_index,ati_to_amber_type);

            topology.num_molecules=number_solute_solvent[1];
			topology.num_solutes=number_solute_solvent[0];
			topology.num_solvents=number_solute_solvent[1]-number_solute_solvent[0];
            
            for (const auto& par : atoms_per_diff_molecule) {
                 //cout << "Clave: " << par.first << ", Atomos: " << get<0>(par.second) << " numero de moleculas: " << get<1>(par.second)<<endl;
                 topology.number_of_atoms_per_different_molecule[par.first]=get<0>(par.second);
                 topology.number_of_each_different_molecule[par.first]=get<1>(par.second);
            }
			
            topology.total_number_of_atoms=dict_pointers["NATOM"];
            map<int, tuple<string, string, float, float>> molecule_atoms;
            vector<map<int, tuple<string, string, float, float>>> atom_type_name_charge_mass;
            int _k=0;

            for(size_t _j=0; _j<number_solute_solvent[1];_j++)
            {
                for (size_t i = 0; i < atoms_per_molecule[_j]; i+=1) {

                    molecule_atoms[i] = make_tuple(ati_to_amber_type[atom_type_index[_k]], get<0>(map_atoms[_k]), get<1>(map_atoms[_k]), mass[_k]);
                    _k++;
                }
                topology.atom_type_name_charge_mass.push_back(molecule_atoms); //molecula numero atomo iniciando en 1 cada molecula
            }
            topology.type_Z=type_atomic_z;
			topology.name_type=name_types; //atom name string  // type es un int
            topology.type_LJparam=lj_diagonal; //0=epsilon 1=sigma //string de numero valor1 valor 2
            topology.special_interaction=lj_coefficient;//keep ij and ji
            //return Configuration::TopolInfo();
            return topology;
            //agregar z agregar el numero de atomos totales
        }
};

class AmberCoordinateReader : public CoordinateReader {
    private:
        /**
         * Given an atom name, return its atomic number
         * @param atom_name The name of the atom
         * @return The atomic number
         */
        /*
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
        }*/
    public:
        AmberCoordinateReader()= default;

        /**
         * Reads the coordinates file .pdb
         * @param filename The name of the coordinates file (single frame)
         * @param topol_info The topology information (use AmberTopologyReader)
         * @param molecs An empty array of molecule pointers
         * @return True if the coordinates were read successfully
         */
        bool readCoordinates_slow(const string& filename, const TopolInfo& topol_info, Molecule** molecs) const {
            ifstream f(filename);
            if(!f.is_open()) {
                cout << "Failed to open file " << filename << endl;
                return false;
            }

            string line;
            vector<Vector> coords;
            map<int,vector<int>> atoms_each_order_molecule;

            while(getline(f,line)) {
                if(line.rfind("ATOM  ",0) == 0 || line.rfind("HETATM",0) == 0) {
                    string record, atom_name, res_name;
                    int atom_id, res_id;
                    float x, y, z;
                    record=ToolKit::strip(line.substr(0, 6));
                    atom_id=stoi(ToolKit::strip(line.substr(6, 5)));
                    atom_name=ToolKit::strip(line.substr(11, 5));
                    res_name=ToolKit::strip(line.substr(16, 4));
                    res_id=stoi(ToolKit::strip(line.substr(20, 6)));
                    x=stof(ToolKit::strip(line.substr(26, 12)));
                    y=stof(ToolKit::strip(line.substr(38, 8)));
                    z=stof(ToolKit::strip(line.substr(46, 8)));
                 
                    atoms_each_order_molecule[res_id].push_back(atom_id);
                    coords.emplace_back(Vector(x,y,z));
                }
            }
            f.close();

            int atom_idx= 0;
            int molec_idx= 0;
            for (const auto& mol_pair : atoms_each_order_molecule)
            {
                const auto& atom_data= topol_info.atom_type_name_charge_mass[mol_pair.first];
                int num_atoms_per_molecule=mol_pair.second.size();
                Atom* atoms= new Atom[num_atoms_per_molecule];

                for(int j= 0; j < num_atoms_per_molecule; j++,atom_idx++) {
                    
                    const auto& [type,name,charge,mass]= atom_data.at(j);
                    const auto& [epsilon,sigma]= topol_info.type_LJparam.at(type);

                    int Z=topol_info.type_Z.at(type);
                    atoms[j]= Atom(coords[atoms_each_order_molecule[mol_pair.first][j]], atom_idx+1, mass, charge, epsilon, sigma, Z);
                }

                molecs[molec_idx]= new Molecule(mol_pair.first, atoms, num_atoms_per_molecule);
                molec_idx+=1;
            }
            return true;
        }

        bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs) const override{
            ifstream f(filename);
            if(!f.is_open()) {
                cout << "Failed to open file " << filename << endl;
                return false;
            }

            string line;
            int number_of_atoms=0;
            for(int i= 0; i < topol_info.num_molecules; i++)
            {
                const auto& atom_data= topol_info.atom_type_name_charge_mass[i];
                Atom* atoms= new Atom[atom_data.size()];

                number_of_atoms=0;
                while(getline(f,line)) {
                    if(line.rfind("TER  ",0) == 0) {break;}
                    if(line.rfind("ATOM  ",0) == 0 || line.rfind("HETATM",0) == 0) {

                        float x=stof(ToolKit::strip(line.substr(26, 12)));
                        float y=stof(ToolKit::strip(line.substr(38, 8)));
                        float z=stof(ToolKit::strip(line.substr(46, 8)));
                        
                        const auto& [type,name,charge,mass]= atom_data.at(number_of_atoms);
                        const auto& [epsilon,sigma]= topol_info.type_LJparam.at(type);
                        int Z=topol_info.type_Z.at(type);

                        atoms[number_of_atoms]= Atom(Vector(x,y,z), number_of_atoms+1, mass, charge, epsilon, sigma, Z);
                        number_of_atoms+=1;
                    }
                }
                molecs[i]= new Molecule(i+1, atoms, number_of_atoms);
            }
            f.close();

            return true;
        }
};

#endif
