#ifndef GROMACS_READERS_HPP
#define GROMACS_READERS_HPP

/**
 * Version: June 2025
 * Author: Nicolás Loubet
 */

#include "ReaderInterfaces.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

class GromacsTopologyReader : public TopologyReader {
    private:
        static inline const unordered_map<string,int> periodic_table= {
            {"H" , 1},                                                                                                                                                                                 {"He", 2},
            {"Li", 3}, {"Be", 4},                                                                                                               {"B" , 5}, {"C" , 6}, {"N" , 7}, {"O" , 8}, {"F" , 9}, {"Ne",10},
            {"Na",11}, {"Mg",12},                                                                                                               {"Al",13}, {"Si",14}, {"P" ,15}, {"S" ,16}, {"Cl",17}, {"Ar",18},
            {"K" ,19}, {"Ca",20}, {"Sc",21}, {"Ti",22}, {"V" ,23}, {"Cr",24}, {"Mn",25}, {"Fe",26}, {"Co",27}, {"Ni",28}, {"Cu",29}, {"Zn",30}, {"Ga",31}, {"Ge",32}, {"As",33}, {"Se",34}, {"Br",35}, {"Kr",36},
            {"Rb",37}, {"Sr",38}, {"Y" ,39}, {"Zr",40}, {"Nb",41}, {"Mo",42}, {"Tc",43}, {"Ru",44}, {"Rh",45}, {"Pd",46}, {"Ag",47}, {"Cd",48}, {"In",49}, {"Sn",50}, {"Sb",51}, {"Te",52}, {"I" ,53}, {"Xe",54},
            {"Cs",55}, {"Ba",56}, // ...
        };

        /**
         * Given a file, read the molecule type in the next line.
         * @param file The file to read
         * @return The molecule type
         */
        static string readMolecType(ifstream &file) {
            string line= "";
            while(getline(file, line))
                if(line[0] != ';')
                    return ToolKit::strip(line.substr(0,8));
            return "ERROR";
        }

        /**
         * Given a file, find the position of the flag in the file.
         * The returned map contains the flag name as the key and the position as the value.
         * The position is the byte position of the first character after the flag.
         * @param file The file to search
         * @return A map with the position of the flag.
         */
        static map<string,int> flagPositions(ifstream &file) {
            string line= "";
            map<string,int> flag;
            string molec_type= "";

            while(getline(file, line)) {
                size_t pos_flag= line.find('[');
                if(pos_flag == string::npos) continue;
                if(line[0] == ';') continue;
                
                line= line.substr(pos_flag, line.find(']')-pos_flag+1);
                if(line == "[ moleculetype ]") {
                    molec_type= "_"+readMolecType(file);
                    continue;
                } else if(line == "[ system ]") molec_type= "";
                line+= molec_type;
                flag[line]= file.tellg();
            }
            return flag;
        }

        /**
         * Given a file, read the system flag.
         * @param file The file to read
         * @param position The position of the flag
         * @return The system flags mapped string -> int
         */
        static map<string,int> readSystemFlag(ifstream &file, int position) {
            map<string,int> molecules;
            string line= "";
            file.clear();
            file.seekg(position);

            while(getline(file, line)) {
                if(line[0] == ';') continue;
                string molec_name= ToolKit::strip(line.substr(0,8));
                int num_molec= stoi(line.substr(9,line.length()-9));
                molecules[molec_name]= num_molec;
            }
            return molecules;
        }

        /**
         * Given a map, sum the number of molecules in the map.
         * @param molecules The map with the number of molecules
         * @return The total number of molecules
         */
        static int sumMoleculesInMap(map<string,int> molecules) {
            int sum= 0;
            for(auto it= molecules.begin(); it != molecules.end(); it++)
                sum+= it->second;
            return sum;
        }

        /**
         * Given a map, sum the number of molecules considered solvent.
         * @param molecules The map with the number of molecules
         * @return The total number of molecules considered solvent
         */
        static int sumMoleculesConsideredSolvent(map<string,int> molecules) {
            int sum= 0;
            for(auto it= molecules.begin(); it != molecules.end(); it++)
                if((it->first.find("SOL") != string::npos) || (it->first.find("WAT") != string::npos))
                    sum+= it->second;
            return sum;
        }

        /**
         * Given a file, read the atoms flags.
         * @param file The file to read
         * @param position The position of the flag
         * @return The atom info mapped int -> tuple (type, atom, charge, mass)
         */
        static map<int,tuple<string,string,float,float>> readAtomsFlags(ifstream &file, int position) {
            map<int,tuple<string,string,float,float>> atoms;
            string line= "";
            file.clear();
            file.seekg(position);

            while(getline(file, line)) {
                if(line[0] == ';') continue;
                if(ToolKit::strip(line) == "") break;

                stringstream ss(line.substr(0,line.find(";")));

                int nr,resi,cgnr;
                string type,res,atom;
                float charge,mass;

                ss >> nr >> type >> resi >> res >> atom >> cgnr >> charge >> mass;
                atoms[nr]= make_tuple(type, atom, charge, mass);
            }
            return atoms;
        }

        static int getZFromName(string name) {
            if(name.length() >= 2) {
                auto it= periodic_table.find(name.substr(0,2));
                if(it != periodic_table.end())
                    return it->second;
            }
            if(name.length() >= 1) {
                auto it= periodic_table.find(name.substr(0,1));
                if(it != periodic_table.end())
                    return it->second;
            }
            return -1;
        }

        static map<string,pair<float,float>> readLJParameters(ifstream &file, int position) {
            map<string,pair<float,float>> parameters;
            string line= "";
            file.clear();
            file.seekg(position);

            while(getline(file, line)) {
                if(line[0] == ';') continue;
                if(ToolKit::strip(line) == "") break;

                stringstream ss(line.substr(0,line.find(";")));

                string type1,type2;
                char ptype;
                float sigma,epsilon,mass,q;

                ss >> type1 >> type2 >> mass >> q >> ptype >> sigma >> epsilon;
                if(!ss.fail()) {
                    parameters[type1]= make_pair(epsilon,sigma);
                    continue;
                }
                
                ss >> type1 >> mass >> q >> ptype >> sigma >> epsilon;
                if(!ss.fail()) {
                    parameters[type1]= make_pair(epsilon,sigma);
                    continue;
                }
                throw runtime_error("Error reading LJ parameters from GROMACS topology file.");
            }
            return parameters;
        }

        static map<pair<string,string>,pair<float,float>> readSpecialInteractions(ifstream &file, int position) {
            map<pair<string,string>,pair<float,float>> parameters;
            string line= "";
            file.clear();
            file.seekg(position);

            while(getline(file, line)) {
                if(line[0] == ';') continue;
                if(ToolKit::strip(line) == "") break;

                stringstream ss(line.substr(0,line.find(";")));

                string type1,type2;
                int func;
                float sigma,epsilon;

                ss >> type1 >> type2 >> func >> sigma >> epsilon;
                if(!ss.fail()) {
                    parameters[make_pair(type1,type2)]= make_pair(epsilon,sigma);
                    continue;
                }
                throw runtime_error("Error reading LJ special parameters from GROMACS topology file.");
            }
            return parameters;
        }


    public:
        GromacsTopologyReader()= default;

        TopolInfo readTopology(const std::string& filename) const override {
            TopolInfo ti= TopolInfo();

            string line;
            ifstream f(filename);

            if(!f.is_open()) {
                std::cerr << "Topology not found" << std::endl;
                return ti;
            }

            map<string,int> flags= flagPositions(f);
            ti.number_of_each_different_molecule= readSystemFlag(f, flags["[ molecules ]"]);
            ti.num_molecules= sumMoleculesInMap(ti.number_of_each_different_molecule);
            ti.num_solvents= sumMoleculesConsideredSolvent(ti.number_of_each_different_molecule);
            ti.num_solutes= ti.num_molecules-ti.num_solvents;
            
            for(auto it_molec= ti.number_of_each_different_molecule.begin(); it_molec != ti.number_of_each_different_molecule.end(); it_molec++) {
                map<int,tuple<string,string,float,float>> atoms= readAtomsFlags(f, flags["[ atoms ]_"+it_molec->first]);
                ti.number_of_atoms_per_different_molecule[it_molec->first]= atoms.size();
                ti.atom_type_name_charge_mass.push_back(atoms);
                for(auto it_atom= atoms.begin(); it_atom != atoms.end(); it_atom++) {
                    string type= get<0>(it_atom->second);
                    string name= it_molec->first+":"+get<1>(it_atom->second);
                    ti.name_type[name]= type;
                    ti.type_Z[type]= getZFromName(name);
                }
            }

            ti.type_LJparam= readLJParameters(f, flags["[ atomtypes ]"]);
            if(flags.find("[ nonbond_params ]") == flags.end())
                ti.special_interaction= map<pair<string,string>,pair<float,float>>();
            else
                ti.special_interaction= readSpecialInteractions(f, flags["[ nonbond_params ]"]);

            return ti;
        }
};

class GromacsCoordinateReader : public CoordinateReader {
    public:
        GromacsCoordinateReader()= default;
        bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
            // Implementación para leer .gro o .trr de GROMACS
            // Los datos se guardan en el Molecule**
            return true;
        }
};

#endif
