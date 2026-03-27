#ifndef LAMMPS_READERS_HPP
#define LAMMPS_READERS_HPP

/**
* Version: June 2025
* Author: Nicolás Loubet
*/

#ifdef USE_VECTOR_TOPOLOGY
static_assert(false, "USE_VECTOR_TOPOLOGY is not compatible with LammpsReaders. "
                     "Disable it or exclude LammpsReaders from your build.");
#endif

#include "ReaderInterfaces.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include "General/ToolKit.hpp" 
#include "Molecules/Water.hpp"

class LammpsTopologyReader : public TopologyReader {
	private:

		/**
		 * Accumulates the masses of the atoms in the LAMMPS data file.
		 */
		struct MolAccumulator {
			vector<Real> masses;
			int atom_count= 0;
		};

		/**
		 * Finds the positions of the sections in the LAMMPS data file.
		 * @param file The LAMMPS data file.
		 * @return A map of section names to their positions in the file.
		 */
		inline map<string,streampos> findSectionPositions(ifstream &file) const {
			map<string,streampos> section_positions;
			string line;
			file.clear();
			file.seekg(0);

			while(getline(file,line)) {
				if(line.find("Masses") != string::npos && line.find("Masses") == 0) {
					section_positions["Masses"]= file.tellg();
				} else if(line.find("Atoms") != string::npos && line.find("Atoms") == 0) {
					section_positions["Atoms"]= file.tellg();
				} else if(line.find("Pair Coeff") != string::npos) {
					section_positions["Pair Coeff"]= file.tellg();
				} else if(line.find("Special Interactions") != string::npos) {
					section_positions["Special Interactions"]= file.tellg();
				} else if(line.find("xlo xhi") != string::npos) {
					section_positions["BoxDimensions"]= file.tellg();
				}
			}
			return section_positions;
		}

		/**
		 * Read the masses from the 'Masses' section of the LAMMPS data file.
		 * @param file The LAMMPS data file.
		 * @param position The position of the 'Masses' section in the file.
		 * @return A map of atom type IDs to their masses.
		 */
		inline map<string,pair<Real,string>> readMasses(ifstream& file, streampos position) const {
			map<string,pair<Real,string>> masses;
			file.clear();
			file.seekg(position);
			string line;
			while(getline(file,line)) {
				string stripped_line= ToolKit::strip(line);

				if(stripped_line.empty() || stripped_line[0] == '#') continue;
				if(stripped_line.find("Atoms") != string::npos) break;

				string line_before_comment= line.substr(0,line.find("#"));
				string line_after_comment= line.substr(line.find("#")+1);

				stringstream ss(line_before_comment);
				int atom_type_id;
				Real mass_val;
				ss >> atom_type_id >> mass_val;
				string type_name= "at"+to_string(atom_type_id);
				if(ss.fail()) {
					throw runtime_error("Error reading masses from LAMMPS data file.");
				}

				if(line_after_comment.at(0) == ' ') line_after_comment= line_after_comment.substr(1);
				string atom_name= "Atom"+line_after_comment.substr(0,line_after_comment.find(" "))+"_"+to_string(atom_type_id);
				masses[type_name]= make_pair(mass_val,atom_name);
			}
			return masses;
		}

		/**
		 * Initialize the topology information.
		 * @return The initialized topology information.
		 */
		inline static TopolInfo initTopology() {
			TopolInfo topology;
			topology.total_number_of_atoms= 0;
			topology.num_molecules= 0;
			topology.num_solvents= 0;
			topology.num_solutes= 0;

			topology.number_of_each_different_molecule= map<string,int>();
			topology.number_of_atoms_per_different_molecule= map<string,int>();
			topology.atom_type_name_charge_mass= vector<map<int,tuple<string,string,Real,Real>>>();
			topology.name_type= map<string,string>();
			topology.type_Z= map<string,int>();

			topology.atom_type_name_charge_mass.push_back(map<int,tuple<string,string,Real,Real>>());
			return topology;
		}

		/**
		 * Read the atom data from the 'Atoms' section of the LAMMPS data file.
		 * @param line The line containing the atom data.
		 * @return A tuple containing the atom ID, molecule ID, atom type, charge, x, y, and z coordinates.
		 */
		inline static tuple<int,int,int,Real,Real,Real,Real> readAtomData(string& line) {
			stringstream ss(line);
			int atom_id, mol_id, atom_type_int;
			Real q, x, y, z;
			string atom_style_type;
			ss >> atom_id >> mol_id >> atom_type_int >> q >> x >> y >> z;
			
			if(ss.fail()) {
				throw runtime_error("Error reading atoms from LAMMPS data file.");
			}
			return make_tuple(atom_id,mol_id,atom_type_int,q,x,y,z);
		}

		/**
		 * Find the atom name and mass for a given atom type.
		 * @param atom_type_int The atom type ID.
		 * @param atomtype_mass_name A map of atom type IDs to their masses and names.
		 * @return A tuple containing the atom type, atom name, and mass.
		 */
		inline static tuple<string,string,Real> findAtomNameMass(int atom_type_int, map<string,pair<Real,string>> atomtype_mass_name) {
			string type_str= "at"+to_string(atom_type_int);
			Real mass= 0.0;
			string atom_name_str= "ERROR";
			if(atomtype_mass_name.count(type_str)) {
				mass= atomtype_mass_name.at(type_str).first;
				atom_name_str= atomtype_mass_name.at(type_str).second;
			} else {
				throw runtime_error("Missing mass for atom type " + type_str + ".");
			}
			return make_tuple(type_str,atom_name_str,mass);
		}

		/**
		 * Check if a molecule is a water molecule.
		 * @param mol The molecule information.
		 * @return True if the molecule was a water molecule, false otherwise.
		 */
		inline static bool isWaterMolecule(const MolAccumulator& mol) {
			if(mol.atom_count < 3 || mol.atom_count > 5) return false;

			int nO= 0, nH= 0;
			for(Real m: mol.masses) {
				if(int(m + 0.2) == 16)
					nO++;
				else if(int(m + 0.2) == 1)
					nH++;
			}
			return (nO == 1 && nH == 2);
		}

		/**
		 * Read the atoms from the 'Atoms' section of the LAMMPS data file.
		 * @param file The file stream to read from.
		 * @param position The position in the file to start reading from.
		 * @param atomtype_mass_name A map of atom type IDs to their masses and names.
		 * @return The topology information.
		 */
		inline TopolInfo readAtoms(ifstream& file, streampos position, const map<string,pair<Real,string>>& atomtype_mass_name) const {
			int internal_mol_id= 1; // internal mol_id, starts at 1 and increments with each new molecule

			file.clear();
			file.seekg(position);
			string line;
			getline(file,line); // comment

			TopolInfo topology= initTopology();

			map<int,MolAccumulator> mol_acc; // Accumulates the masses of the atoms in each molecule
			map<int,int> mol_atom_counter; // Counts the number of atoms in each molecule

			while(getline(file,line)) {
				if(ToolKit::strip(line).empty()) break;

				int atom_id, mol_id, atom_type_int;
				Real q, x, y, z;
				tie(atom_id,mol_id,atom_type_int,q,x,y,z)= readAtomData(line);

				string type_str, atom_name_str;
				Real mass;
				tie(type_str,atom_name_str,mass)= findAtomNameMass(atom_type_int, atomtype_mass_name);

				topology.atom_type_name_charge_mass[0][atom_id]= make_tuple(type_str, atom_name_str, q, mass);
				topology.name_type[atom_name_str]= type_str;
				topology.total_number_of_atoms++;

				mol_acc[mol_id].masses.push_back(mass);
				mol_acc[mol_id].atom_count++;
			}

			string water_molecule_type= "mtWAT";
			bool water_type_defined= false;

			for(const auto& [mol_id, mol]: mol_acc) {
				topology.num_molecules++;

				if(isWaterMolecule(mol)) {
					if(!water_type_defined) {
						topology.number_of_each_different_molecule[water_molecule_type]= 0;
						topology.number_of_atoms_per_different_molecule[water_molecule_type]= mol.atom_count;
						water_type_defined= true;
					}
					topology.num_solvents++;
					topology.number_of_each_different_molecule[water_molecule_type]++;
					topology.name_type["Mol"+to_string(internal_mol_id++)]= water_molecule_type;
				}
				else {
					topology.num_solutes++;
					string mt= "mt"+to_string(mol_id);
					topology.number_of_each_different_molecule[mt]= 1;
					topology.number_of_atoms_per_different_molecule[mt]= mol.atom_count;
					topology.name_type["Mol"+to_string(internal_mol_id++)]= mt;
				}
			}
			return topology;
		}

		/**
		 * Get the keys of a map that start with a certain prefix.
		 * @param map_to_search The map to search.
		 * @param prefix The prefix to search for.
		 * @return A vector of keys that start with the prefix.
		 */
		inline static vector<string> getKeysStartingWith(map<string,string>& map_to_search, const string& prefix) {
			vector<string> atom_keys;
			for(const auto& [name,type]: map_to_search)
				if(name.size() >= 4 && name.substr(0,4) == prefix)
					atom_keys.push_back(name.substr(4));
			return atom_keys;
		}

		/**
		 * Map the atomic numbers of the atoms.
		 * @param atom_keys A vector of atom keys.
		 * @param map_name_type A map of atom names to their types.
		 * @return A map of atom names to their atomic numbers.
		 */
		inline static map<string,int> mapZValues(vector<string> atom_keys, map<string,string> map_name_type) {
			map<string,int> Z_values;
			for(string& k: atom_keys) {
				string key= k;
				transform(key.begin(), key.end(), key.begin(), ::toupper);

				if(periodic_table.count(key)) {
					Z_values[map_name_type.at("Atom"+key)]= periodic_table.at(key);
					continue;
				}

				// Remove numeric digits
				string key_without_digits= key;
				key_without_digits.erase(remove_if(key_without_digits.begin(), key_without_digits.end(), ::isdigit), key_without_digits.end());
				key_without_digits.erase(remove(key_without_digits.begin(), key_without_digits.end(), '_'), key_without_digits.end());
				if(periodic_table.count(key_without_digits)) {
					Z_values[map_name_type.at("Atom"+k)]= periodic_table.at(key_without_digits);
					continue;
				}

				Z_values[key]= -1;
			}
			return Z_values;
		}

		/**
		 * Reads the LJ parameters from the 'Pair Coeffs' section of the LAMMPS data file.
		 * @param file The LAMMPS data file.
		 * @param position The position of the 'Pair Coeffs' section in the file.
		 * @return A map of atom type IDs to their epsilon and sigma values.
		 */
		inline map<string,pair<Real,Real>> readLJParameters(ifstream& file, streampos position) const {
			map<string,pair<Real,Real>> lj_params;
			file.clear();
			file.seekg(position);
			string line;

			while(getline(file,line)) {
				if(ToolKit::strip(line) == "") break;

				string line_without_comment= line.substr(1,line.find("##")-1);
				stringstream ss(line_without_comment);
				int atom_type_id;Real epsilon, sigma;
				ss >> atom_type_id >> epsilon >> sigma;

				if(ss.fail()) {
					throw runtime_error("Error reading LJ parameters from LAMMPS data file.");
				}

				lj_params["at"+to_string(atom_type_id)]= make_pair(epsilon,sigma);
			}
			return lj_params;
		}

		/**
		 * Reads the LJ parameters from the 'Special Interactions' section of the LAMMPS data file.
		 * @param file The LAMMPS data file.
		 * @param position The position of the 'Special Interactions' section in the file.
		 * @return A map of atom type IDs (two types) to their epsilon and sigma values.
		 */
		inline map<pair<string,string>,pair<Real,Real>> readSpecialInteractions(ifstream &file, int position) const {
			map<pair<string,string>,pair<Real,Real>> special_interactions;
			file.clear();
			file.seekg(position);
			string line;

			while(getline(file,line)) {
				if(ToolKit::strip(line) == "") break;

				string line_without_comment= line.substr(1,line.find("##")-1);
				stringstream ss(line_without_comment);
				string type1,type2;
				int func;
				Real sigma,epsilon;

				ss >> type1 >> type2 >> func >> epsilon >> sigma;
				if(!ss.fail()) {
					special_interactions[make_pair("at"+type1,"at"+type2)]= make_pair(epsilon,sigma);
					special_interactions[make_pair("at"+type2,"at"+type1)]= make_pair(epsilon,sigma);
					continue;
				}
				throw runtime_error("Error reading LJ special parameters from LAMMPS data file.");
			}
			return special_interactions;
		}

		/**
		 * Reads the box dimensions from the 'BoxDimensions' section of the LAMMPS data file.
		 * @param file The LAMMPS data file.
		 * @return The box dimensions.
		 */
		inline Vector readBoxDimensions(ifstream& file) const {
			file.clear();
			file.seekg(0);
			string line;
			Real xlo=0, xhi=0, ylo=0, yhi=0, zlo=0, zhi=0;

			while(getline(file,line)) {
				if(line.find("xlo xhi") != string::npos) {
					stringstream ss(line);
					ss >> xlo >> xhi;
					break; 
				}
			}
			while(getline(file,line)) {
				if(line.find("ylo yhi") != string::npos) {
					stringstream ss(line);
					ss >> ylo >> yhi;
					break;
				}
			}
			while(getline(file,line)) {
				if(line.find("zlo zhi") != string::npos) {
					stringstream ss(line);
					ss >> zlo >> zhi;
					break;
				}
			}
			return Vector(xhi - xlo, yhi - ylo, zhi - zlo);
		}

	public:
		LammpsTopologyReader()= default;

		/**
		 * Reads the topology information from a LAMMPS data file.
		 * @param filename The name of the LAMMPS data file.
		 * @return A TopolInfo struct containing the topology information.
		 */
		TopolInfo readTopology(const string& filename) const override {
			ifstream file(filename);
			if(!file.is_open()) {
				throw runtime_error("Error opening LAMMPS data file.");
			}

			map<string,streampos> positions_of_sections= findSectionPositions(file);

			map<string,pair<Real,string>> atomtype_mass_name;
			if(positions_of_sections.count("Masses"))
				atomtype_mass_name= readMasses(file, positions_of_sections["Masses"]);

			TopolInfo topology;
			if(positions_of_sections.count("Atoms")) {
				topology= readAtoms(file, positions_of_sections["Atoms"], atomtype_mass_name);
				topology.type_Z= mapZValues(getKeysStartingWith(topology.name_type,"Atom"),topology.name_type);
			}

			if(positions_of_sections.count("Pair Coeff"))
				topology.type_LJparam= readLJParameters(file, positions_of_sections["Pair Coeff"]);
			
			if(positions_of_sections.count("Special Interactions"))
				topology.special_interaction= readSpecialInteractions(file, positions_of_sections["Special Interactions"]);

			if(positions_of_sections.count("BoxDimensions"))
				topology.default_system_bounds= readBoxDimensions(file);

			file.close();
			return topology;
		}
};

class LammpsCoordinateReader : public CoordinateReader {
	public:
		LammpsCoordinateReader() {}
		~LammpsCoordinateReader()= default;

		/**
		 * Reads the coordinates from a LAMMPS data file. Iteration over different frames is done automatically.
		 * @param filename The name of the LAMMPS data file.
		 * @param topol_info The topology information.
		 * @param molecs The array of Molecule pointers to be filled with the coordinates.
		 * @param bounds The bounds of the system.
		 * @return True if the coordinates were successfully read, false otherwise.
		 */
		bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
			ifstream f(filename);
			if(!f.is_open()) {
				cout << "Error: Failed to open file " << filename << endl;
				return false;
			}
			if(molecs == nullptr) return false;

			string line;
			int num_atoms_in_frame= 0;

			if(getline(f,line)) {
				num_atoms_in_frame= stoi(ToolKit::strip(line));
			} else {
				cout << "Error: Failed to read number of atoms in frame " << i_frame << endl;
				f.close();
				return false;
			}
			
			for(int i_read= 0; i_read <= i_frame*(num_atoms_in_frame+2); i_read++)
				getline(f,line); // Skip to frame wanted (also kipping comments)

			vector<Vector> coords(num_atoms_in_frame);
			vector<string> atom_symbols(num_atoms_in_frame);

			for(int i= 0; i < num_atoms_in_frame; i++) {
				if(getline(f,line)) {
					stringstream ss(line);
					string symbol; Real x, y, z;
					ss >> symbol >> x >> y >> z;
					atom_symbols[i]= symbol;
					coords[i]= Vector(x,y,z);
				} else {
					cout << "Error: Failed to read coordinates for frame " << i_frame << "." << endl;
					f.close();
					return false;
				}
			}

			int global_atom_idx= 0;
			for(int i= 0; i < topol_info.num_molecules; i++) {
				int num_atoms_per_molecule= topol_info.number_of_atoms_per_different_molecule.at(topol_info.name_type.at("Mol"+to_string(i+1)));
				Atom* atoms_array_for_molecule= new Atom[num_atoms_per_molecule];

				bool is_water= topol_info.name_type.at("Mol"+to_string(i+1)) == "mtWAT";

				for(int j= 0; j < num_atoms_per_molecule; j++,global_atom_idx++) {
					const auto& [type,atom_name,q,mass]= topol_info.atom_type_name_charge_mass[0].at(global_atom_idx+1);

					Real e= 0.0, s= 0.0;
					if(topol_info.type_LJparam.count(type)) {
						e= topol_info.type_LJparam.at(type).first;
						s= topol_info.type_LJparam.at(type).second;
					} else {
						cout << "Error: LJ parameters not found for type " << type << endl;
						f.close();
						return false;
					}

					int Z= 0;
					if(topol_info.type_Z.count(type)) {
						Z= topol_info.type_Z.at(type);
					} else {
						throw runtime_error("Error: Z value not found for type " + type);
					}

					atoms_array_for_molecule[j]= Atom(coords[global_atom_idx], j+1, mass, q, e, s, Z, type);
				}
				
				if(is_water) {
					molecs[i]= new Water(i+1, atoms_array_for_molecule, num_atoms_per_molecule);
				} else {
					molecs[i]= new Molecule(i+1, atoms_array_for_molecule, num_atoms_per_molecule);
				}

			}

			f.close();
			bounds= topol_info.default_system_bounds;
			return true;
		}
};

#endif // LAMMPS_READERS_HPP
