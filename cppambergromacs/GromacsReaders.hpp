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

class GromacsTopologyReader : public TopologyReader {
    private:
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

        /**
         * Given a string, return the atomic number from the periodic table.
         * @param name The name of the atom
         * @return The atomic number
         */
        static int getZFromName(string name) {
            transform(name.begin(), name.end(), name.begin(), ::toupper);

            if(name.length() >= 2)
                if(periodic_table.count(name.substr(0,2)))
                    return periodic_table.at(name.substr(0,2));
            if(name.length() >= 1)
                if(periodic_table.count(name.substr(0,1)))
                    return periodic_table.at(name.substr(0,1));
            return -1;
        }

        /**
         * Given a file, read the LJ parameters.
         * @param file The file to read
         * @param position The position of the flag
         * @return The LJ parameters mapped type -> pair (epsilon, sigma)
         */
        static map<string,tuple<float,float,float,float>> readLJFlagFully(ifstream &file, int position) {
            map<string,tuple<float,float,float,float>> parameters;
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
                    parameters[type1]= make_tuple(mass,q,epsilon,sigma);
                    continue;
                }
                
                ss >> type1 >> mass >> q >> ptype >> sigma >> epsilon;
                if(!ss.fail()) {
                    parameters[type1]= make_tuple(mass,q,epsilon,sigma);
                    continue;
                }
                throw runtime_error("Error reading LJ parameters from GROMACS topology file.");
            }
            return parameters;
        }

        static void checkMass(map<int,tuple<string,string,float,float>> atoms, map<string,tuple<float,float,float,float>> params) {
            for(auto it= atoms.begin(); it != atoms.end(); it++) {
                string type, name; float charge, mass;
                tie(type,name,charge,mass)= it->second;
                if(mass > 0.8f) continue;

                float mass_LJ, charge_LJ, s, e;
                if(params.count(type) > 0) continue;
                tie(mass_LJ,charge_LJ,e,s)= params.at(type);
                
                if(mass_LJ < 0.8f) continue;
                get<3>(it->second)= mass_LJ;
            }
        }

        /**
         * Given a file, read the LJ special parameters.
         * @param file The file to read
         * @param position The position of the flag
         * @return The LJ special parameters mapped pair (type1,type2) -> pair (epsilon,sigma)
         */
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

        /**
         * Given a file, read the topology.
         * @param filename The name of the file
         * @return A TopolInfo object
         */
        TopolInfo readTopology(const string& filename) const override {
            TopolInfo ti= TopolInfo();

            string line;
            ifstream f(filename);

            if(!f.is_open()) {
                cerr << "Topology not found" << endl;
                return ti;
            }

            map<string,int> flags= flagPositions(f);
            ti.number_of_each_different_molecule= readSystemFlag(f, flags["[ molecules ]"]);
            ti.num_molecules= sumMoleculesInMap(ti.number_of_each_different_molecule);
            ti.num_solvents= sumMoleculesConsideredSolvent(ti.number_of_each_different_molecule);
            ti.num_solutes= ti.num_molecules-ti.num_solvents;

            map<string,tuple<float,float,float,float>> parameters_LJflag_full= readLJFlagFully(f, flags["[ atomtypes ]"]);
            for(auto it= parameters_LJflag_full.begin(); it != parameters_LJflag_full.end(); it++)
                ti.type_LJparam[it->first]= make_pair(get<2>(it->second),get<3>(it->second));
            
            for(auto it_molec= ti.number_of_each_different_molecule.begin(); it_molec != ti.number_of_each_different_molecule.end(); it_molec++) {
                map<int,tuple<string,string,float,float>> atoms= readAtomsFlags(f, flags["[ atoms ]_"+it_molec->first]);
                checkMass(atoms, parameters_LJflag_full);
                ti.number_of_atoms_per_different_molecule[it_molec->first]= atoms.size();
                ti.total_number_of_atoms+= atoms.size()*it_molec->second;
                ti.atom_type_name_charge_mass.push_back(atoms);
                for(auto it_atom= atoms.begin(); it_atom != atoms.end(); it_atom++) {
                    string type= get<0>(it_atom->second);
                    string name= it_molec->first+":"+get<1>(it_atom->second);
                    ti.name_type[name]= type;
                    if(!ti.type_Z.count(type))
                        ti.type_Z[type]= getZFromName(get<1>(it_atom->second));
                }
            }

            if(flags.find("[ nonbond_params ]") == flags.end())
                ti.special_interaction= map<pair<string,string>,pair<float,float>>();
            else
                ti.special_interaction= readSpecialInteractions(f, flags["[ nonbond_params ]"]);

            return ti;
        }
};

class GromacsCoordinateReader : public CoordinateReader {
    private:
        /**
         * Create a new Molecule object
         * @param molec_name The name of the molecule
         * @param i_molec The index of the molecule
         * @param molecs The array of Molecule objects
         * @param atom_list The array of Atom objects
         * @param number_of_atom_in_list The number of atoms in the atom_list
         */
        static void createNewMolecule(string molec_name, int i_molec, Molecule** molecs, Atom* atom_list, int number_of_atom_in_list) {
            if(molec_name == "SOL" || molec_name == "WAT") {
                molecs[i_molec-1]= new Water(i_molec, atom_list, number_of_atom_in_list);
            } else {
                molecs[i_molec-1]= new Molecule(i_molec, atom_list, number_of_atom_in_list);
            }
        }

        /**
         * Check if a new Molecule object must be created
         * @param i_molec The index of the molecule
         * @param molec_name The name of the molecule
         * @param atom_list The array of Atom objects
         * @param number_of_atom_in_list The number of atoms in the atom_list
         * @param topol_info The topology information
         * @param molecs The array of Molecule objects
         * @param previous_molec_name The name of the previous molecule
         * @param previous_different_molec_id The index of the previous different molecule
         * @param previous_molec_id The index of the previous molecule
         */
        static void checkIfNewMolecule(int i_molec, string molec_name, Atom*& atom_list, int& number_of_atom_in_list, const TopolInfo& topol_info, Molecule** molecs, string& previous_molec_name, int& previous_different_molec_id, int& previous_molec_id) {
            if(atom_list != nullptr) {
                createNewMolecule(previous_molec_name, i_molec-1, molecs, atom_list, number_of_atom_in_list);
                number_of_atom_in_list= 0;
            }

            previous_molec_id= i_molec;
            if(molec_name != previous_molec_name) {
                previous_different_molec_id++;
                previous_molec_name= molec_name;
            }

            atom_list= new Atom[topol_info.number_of_atoms_per_different_molecule.at(molec_name)];
        }

        /**
         * Read an Atom object
         * @param line The line to read in .gro format
         * @param topol_info The topology information
         * @param molecs The array of Molecule objects
         * @param previous_molec_name The name of the previous molecule
         * @param previous_different_molec_id The index of the previous different molecule
         * @param previous_molec_id The index of the previous molecule
         * @param atom_list The array of Atom objects
         * @param number_of_atom_in_list The number of atoms in the atom_list
         * @return The Atom object
         */
        static Atom readAtom(string line, const TopolInfo& topol_info, Molecule** molecs, string& previous_molec_name, int& previous_different_molec_id, int& previous_molec_id, Atom*& atom_list, int& number_of_atom_in_list) {
            int i_molec= stoi(line.substr(0,5));
            string molec_name= ToolKit::strip(line.substr(5,5));
            string atom_name= ToolKit::strip(line.substr(10,5));
            int i_atom= stoi(line.substr(15,5));
            float x= stof(line.substr(20,8));
            float y= stof(line.substr(28,8));
            float z= stof(line.substr(36,8));

            if(i_molec != previous_molec_id)
                checkIfNewMolecule(i_molec, molec_name, atom_list, number_of_atom_in_list, topol_info, molecs, previous_molec_name, previous_different_molec_id, previous_molec_id);

            string type, name; float q, mass, e, s;
            tie(type,name,q,mass)= topol_info.atom_type_name_charge_mass[previous_different_molec_id].at(number_of_atom_in_list+1);
            int Z= topol_info.type_Z.at(type);
            tie(e,s)= topol_info.type_LJparam.at(type);

            return Atom(Vector(x*10,y*10,z*10), i_atom, mass, q, e, s, Z);
        }

        /**
         * Read the bounds from the .gro file
         * @param line The line to read
         * @return The bounds as a Vector object
         */
        static Vector readBounds(string line) {
            float x,y,z;
            stringstream ss(line);
            ss >> x >> y >> z;
            return Vector(x*10,y*10,z*10);
        }

    public:
        GromacsCoordinateReader()= default;

        /**
         * Read the coordinates from a .gro file
         * @param filename The name of the file
         * @param topol_info The topology information
         * @param molecs The array of Molecule objects to be created
         * @param bounds The bounds of the system to be read
         * @return True if the coordinates were read successfully
         */
        bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
            ifstream f(filename);
            string line;

            getline(f, line); // Title
            getline(f, line);
            int natoms= stoi(line);

            string previous_molec_name= "ERRORMOLECULE";
            int previous_different_molec_id= -1;
            int previous_molec_id= -1;
            Atom* atom_list= nullptr;
            int number_of_atom_in_list= 0;

            for(int i= 0; i < natoms; i++) {
                getline(f, line);
                atom_list[number_of_atom_in_list++]= readAtom(line, topol_info, molecs, previous_molec_name, previous_different_molec_id, previous_molec_id, atom_list, number_of_atom_in_list);
            }

            if(atom_list != nullptr) // Add last molecule read
                createNewMolecule(previous_molec_name, previous_molec_id, molecs, atom_list, number_of_atom_in_list);

            getline(f, line);
            bounds= readBounds(line);

            return true;
        }
};

#endif
