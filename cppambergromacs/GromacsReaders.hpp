#ifndef GROMACS_READERS_HPP
#define GROMACS_READERS_HPP

/**
 * Version: July 2025
 * Author: Nicolás Loubet
 */

#include "ReaderInterfaces.hpp"
#include <string>
#include <fstream>
#include <sstream>

class GromacsTopologyReader : public TopologyReader {
	private:
		/**
		 * Reads the molecule type name from the next non-comment line.
		 */
		static string readMolecType(ifstream& file) {
			string line;
			while(getline(file, line))
				if(line[0] != ';')
					return ToolKit::strip(line.substr(0, 6));
			return "ERROR";
		}

		/**
		 * Scans the file and records byte offsets for each [...] section flag.
		 * Molecule-specific sections are suffixed with "_<MoleculeName>".
		 */
		static map<string, int> flagPositions(ifstream& file) {
			map<string, int> flag;
			string line, molec_type;

			while(getline(file, line)) {
				size_t lp= line.find('[');
				if(lp == string::npos || line[0] == ';') continue;

				line= line.substr(lp, line.find(']') - lp + 1);
				if(line == "[ moleculetype ]") {
					molec_type= "_" + readMolecType(file);
					continue;
				} else if(line == "[ system ]") {
					molec_type= "";
				}
				flag[line + molec_type]= file.tellg();
			}
			return flag;
		}

		/**
		 * Reads the [ molecules ] section.
		 * Returns {name -> count} and the insertion-order vector of names.
		 */
		static pair<map<string,int>, vector<string>> readSystemFlag(ifstream& file, int position) {
			map<string, int> molecules;
			vector<string> order;
			string line;

			file.clear();
			file.seekg(position);
			while(getline(file, line)) {
				if(line.empty() || line[0] == ';') continue;
				string name= ToolKit::strip(line.substr(0, 8));
				molecules[name]= stoi(line.substr(9));
				order.push_back(name);
			}
			return {molecules, order};
		}

		static int sumMoleculesInMap(const map<string,int>& molecules) {
			int sum= 0;
			for(const auto& [name, count] : molecules) sum+= count;
			return sum;
		}

#ifdef DIFF_RESID_MOLECULE
		/**
		 * Reads the [ molecules ] section preserving the exact per-line order and
		 * counts (unlike readSystemFlag, which collapses everything into a map and
		 * loses both order and duplicate blocks). This is the ground truth used to
		 * detect molecule boundaries in the .gro by atom count.
		 */
		static vector<pair<string,int>> readMoleculeSequence(ifstream& file, int position) {
			vector<pair<string,int>> sequence;
			string line;
			file.clear();
			file.seekg(position);
			while(getline(file, line)) {
				if(line.empty() || line[0] == ';') continue;
				string name= ToolKit::strip(line.substr(0, 8));
				int count= stoi(line.substr(9));
				sequence.emplace_back(name, count);
			}
			return sequence;
		}
#endif

		static int sumMoleculesConsideredSolvent(const map<string,int>& molecules) {
			int sum= 0;
			for(const auto& [name, count] : molecules)
				if(name.find("SOL") != string::npos || name.find("WAT") != string::npos)
					sum+= count;
			return sum;
		}

		/**
		 * Reads the [ atoms ] section for one molecule type.
		 * Returns atom index -> (type, name, charge, mass).
		 */
#ifdef USE_VECTOR_TOPOLOGY
		static vector<tuple<string,string,Real,Real>> readAtomsFlags(ifstream& file, int position) {
			vector<tuple<string,string,Real,Real>> atoms;
#else
		static map<int,tuple<string,string,Real,Real>> readAtomsFlags(ifstream& file, int position) {
			map<int,tuple<string,string,Real,Real>> atoms;
#endif
			string line;
			file.clear();
			file.seekg(position);

			while(getline(file, line)) {
				if(line[0] == ';') continue;
				if(ToolKit::strip(line).empty()) break;

				stringstream ss(line.substr(0, line.find(';')));
				int nr, resi, cgnr;
				string type, res, atom;
				Real charge, mass;
				ss >> nr >> type >> resi >> res >> atom >> cgnr >> charge >> mass;

#ifdef USE_VECTOR_TOPOLOGY
				atoms.emplace_back(type, atom, charge, mass);
#else
				atoms[nr]= make_tuple(type, atom, charge, mass);
#endif
			}
			return atoms;
		}

		static int getZFromName(string name) {
			transform(name.begin(), name.end(), name.begin(), ::toupper);
			if(name.length() >= 2 && periodic_table.count(name.substr(0, 2)))
				return periodic_table.at(name.substr(0, 2));
			if(name.length() >= 1 && periodic_table.count(name.substr(0, 1)))
				return periodic_table.at(name.substr(0, 1));
			return -1;
		}

		/**
		 * Reads the [ atomtypes ] section.
		 * Returns type -> (mass, charge, epsilon, sigma[Å]).
		 * Handles both 7-column (with bond type) and 6-column formats.
		 */
		static map<string,tuple<Real,Real,Real,Real>> readLJFlagFully(ifstream& file, int position) {
			map<string, tuple<Real,Real,Real,Real>> parameters;
			string line;
			file.clear();
			file.seekg(position);

			while(getline(file, line)) {
				if(line[0] == ';') continue;
				if(ToolKit::strip(line).empty()) break;

				stringstream ss(line.substr(0, line.find(';')));
				string type1, type2;
				char ptype;
				Real sigma, epsilon, mass, q;

				// Try 7-column format: type bondtype mass q ptype sigma epsilon
				ss >> type1 >> type2 >> mass >> q >> ptype >> sigma >> epsilon;
				if(!ss.fail()) {
					parameters[type1]= {mass, q, epsilon, sigma * 10};
					continue;
				}

				// Try 6-column format: type mass q ptype sigma epsilon
				ss.clear();
				ss.str(line.substr(0, line.find(';')));
				ss >> type1 >> mass >> q >> ptype >> sigma >> epsilon;
				if(!ss.fail()) {
					parameters[type1]= {mass, q, epsilon, sigma * 10};
					continue;
				}

				throw runtime_error("Error reading LJ parameters from GROMACS topology file.");
			}
			return parameters;
		}

		/**
		 * Patches atoms whose mass is ~0 in [ atoms ] using the value from [ atomtypes ].
		 */
#ifdef USE_VECTOR_TOPOLOGY
		static void checkMass(vector<tuple<string,string,Real,Real>>& atoms, const map<string,tuple<Real,Real,Real,Real>>& params) {
			for(auto& [type, name, charge, mass] : atoms) {
				if(mass > 0.8) continue;
				auto it= params.find(type);
				if(it == params.end()) continue;
				Real mass_LJ= get<0>(it->second);
				if(mass_LJ > 0.8) mass= mass_LJ;
			}
		}
#else
		static void checkMass(map<int,tuple<string,string,Real,Real>>& atoms, const map<string,tuple<Real,Real,Real,Real>>& params) {
			for(auto& [idx, data] : atoms) {
				auto& [type, name, charge, mass]= data;
				if(mass > 0.8) continue;
				auto it= params.find(type);
				if(it == params.end()) continue;
				Real mass_LJ= get<0>(it->second);
				if(mass_LJ > 0.8) mass= mass_LJ;
			}
		}
#endif

		/**
		 * Reads the [ nonbond_params ] section.
		 * Returns (type1, type2) -> (epsilon, sigma).
		 */
		static map<pair<string,string>,pair<Real,Real>> readSpecialInteractions(ifstream& file, int position) {
			map<pair<string,string>, pair<Real,Real>> parameters;
			string line;
			file.clear();
			file.seekg(position);

			while(getline(file, line)) {
				if(line[0] == ';') continue;
				if(ToolKit::strip(line).empty()) break;

				stringstream ss(line.substr(0, line.find(';')));
				string type1, type2;
				int func;
				Real sigma, epsilon;
				ss >> type1 >> type2 >> func >> sigma >> epsilon;
				if(!ss.fail()) {
					parameters[{type1, type2}]= {epsilon, sigma};
					continue;
				}
				throw runtime_error("Error reading LJ special parameters from GROMACS topology file.");
			}
			return parameters;
		}

	public:
		GromacsTopologyReader()= default;

		TopolInfo readTopology(const string& filename) const override {
			TopolInfo ti;
			ifstream f(filename);
			if(!f.is_open()) { cerr << "Topology not found" << endl; return ti; }

			map<string,int> flags= flagPositions(f);

			auto [mol_counts, mol_order]= readSystemFlag(f, flags["[ molecules ]"]);
			ti.number_of_each_different_molecule= mol_counts;
			ti.num_molecules= sumMoleculesInMap(mol_counts);
			ti.num_solvents = sumMoleculesConsideredSolvent(mol_counts);
			ti.num_solutes  = ti.num_molecules - ti.num_solvents;

			auto lj_full= readLJFlagFully(f, flags["[ atomtypes ]"]);
			for(const auto& [type, params] : lj_full)
				ti.type_LJparam[type]= {get<2>(params), get<3>(params)};

			for(const auto& molec_name : mol_order) {
				auto atoms= readAtomsFlags(f, flags["[ atoms ]_" + molec_name]);
				checkMass(atoms, lj_full);

				ti.number_of_atoms_per_different_molecule[molec_name]= atoms.size();
				ti.total_number_of_atoms+= atoms.size() * ti.number_of_each_different_molecule.at(molec_name);
				ti.atom_type_name_charge_mass.push_back(move(atoms));

#ifdef DIFF_RESID_MOLECULE
				if(!ti.name_to_diffid.count(molec_name))
					ti.name_to_diffid[molec_name]= (int)ti.atom_type_name_charge_mass.size() - 1;
#endif

#ifdef USE_VECTOR_TOPOLOGY
				for(const auto& [type, name, charge, mass] : ti.atom_type_name_charge_mass.back()) {
					ti.name_type[molec_name + ":" + name]= type;
					if(!ti.type_Z.count(type)) ti.type_Z[type]= getZFromName(name);
				}
#else
				for(const auto& [idx, data] : ti.atom_type_name_charge_mass.back()) {
					const auto& [type, name, charge, mass]= data;
					ti.name_type[molec_name + ":" + name]= type;
					if(!ti.type_Z.count(type)) ti.type_Z[type]= getZFromName(name);
				}
#endif
			}

			ti.special_interaction= (flags.count("[ nonbond_params ]"))
				? readSpecialInteractions(f, flags["[ nonbond_params ]"])
				: map<pair<string,string>,pair<Real,Real>>{};

#ifdef DIFF_RESID_MOLECULE
			ti.molecule_sequence= readMoleculeSequence(f, flags["[ molecules ]"]);
#endif

			return ti;
		}
};

// =============================================================================
//  COORDINATE READER
// =============================================================================

class GromacsCoordinateReader : public CoordinateReader {
	private:
		static void createNewMolecule(const string& molec_name, int i_molec, Molecule** molecs, Atom* atom_list, int n_atoms) {
			if(molec_name == "SOL" || molec_name == "WAT") molecs[i_molec-1]= new Water(i_molec, atom_list, n_atoms);
			else                                           molecs[i_molec-1]= new Molecule(i_molec, atom_list, n_atoms);
		}

		static void checkIfNewMolecule(int i_molec, const string& molec_name, Atom*& atom_list, int& n_atoms, const TopolInfo& topol_info,
		                               Molecule** molecs, string& prev_molec_name, int& prev_diff_id, int& prev_molec_id) {
			if(atom_list != nullptr) {
				if(i_molec < prev_molec_id) i_molec= prev_molec_id + 1;
				createNewMolecule(prev_molec_name, i_molec - 1, molecs, atom_list, n_atoms);
				n_atoms= 0;
			}
			prev_molec_id= i_molec;
			if(molec_name != prev_molec_name) {
				prev_diff_id++;
				prev_molec_name= molec_name;
			}

			auto it= topol_info.number_of_atoms_per_different_molecule.find(molec_name);
			atom_list= new Atom[it->second];
		}

#ifdef DIFF_RESID_MOLECULE
		/**
		 * Tracks progress through topol_info.molecule_sequence while scanning the
		 * .gro file, so molecule boundaries are decided by "did we read as many
		 * atoms as the topology says this molecule has", not by resname/resid.
		 */
		struct MoleculeCursor {
			size_t block_idx= 0;         // index into topol_info.molecule_sequence
			int    instance_in_block= 0; // which instance within that block
			int    molecule_index= 0;    // 0-based count of molecules created so far
			string current_name;
			int    current_diff_id= -1;
			int    expected_atoms= 0;
		};

		static void startNextMolecule(MoleculeCursor& cur, const TopolInfo& topol_info) {
			cur.current_name= topol_info.molecule_sequence[cur.block_idx].first;
			auto it_id= topol_info.name_to_diffid.find(cur.current_name);
			cur.current_diff_id= it_id->second;
			auto it_n= topol_info.number_of_atoms_per_different_molecule.find(cur.current_name);
			cur.expected_atoms= it_n->second;
		}

		static void advanceMoleculeCursor(MoleculeCursor& cur, const TopolInfo& topol_info) {
			cur.molecule_index++;
			cur.instance_in_block++;
			if(cur.instance_in_block >= topol_info.molecule_sequence[cur.block_idx].second) {
				cur.block_idx++;
				cur.instance_in_block= 0;
			}
		}

		/**
		 * Same purpose as readAtom(), but decides molecule membership purely from
		 * atom counts (via MoleculeCursor) instead of resname/resid. resid/resname
		 * are still parsed here only for diagnostic logging.
		 */
		static Atom readAtomByCount(const string& line, const TopolInfo& topol_info, MoleculeCursor& cur,
		                            Atom*& atom_list, int& n_atoms, int atom_line_index) {
			int    i_molec   = stoi(line.substr(0,  5));
			string molec_name= ToolKit::strip(line.substr(5,  5));
			string atom_name = ToolKit::strip(line.substr(10, 5));
			int    i_atom    = stoi(line.substr(15, 5));
			Real   x         = RealParser(line.substr(20, 8));
			Real   y         = RealParser(line.substr(28, 8));
			Real   z         = RealParser(line.substr(36, 8));

			if(atom_list == nullptr) {
				startNextMolecule(cur, topol_info);
				atom_list= new Atom[cur.expected_atoms];
				n_atoms= 0;
			}

			const auto& molec_atoms= topol_info.atom_type_name_charge_mass[cur.current_diff_id];
			int lookup_idx=
#ifdef USE_VECTOR_TOPOLOGY
				n_atoms;
#else
				n_atoms + 1;
#endif

			const auto& [type, name, q, mass]= molec_atoms.at(lookup_idx);
			const auto& [e,s]= topol_info.type_LJparam.at(type);
			int Z= topol_info.type_Z.at(type);
			return Atom(Vector(x*10, y*10, z*10), i_atom, mass, q, e, s, Z, type);
		}
#endif

		static Atom readAtom(const string& line, const TopolInfo& topol_info, Molecule** molecs, string& prev_molec_name, int& prev_diff_id,
		                     int& prev_molec_id, Atom*& atom_list, int& n_atoms, int atom_line_index) {
			int    i_molec   = stoi(line.substr(0,  5));
			string molec_name= ToolKit::strip(line.substr(5,  5));
			string atom_name = ToolKit::strip(line.substr(10, 5));
			int    i_atom    = stoi(line.substr(15, 5));
			Real   x         = RealParser(line.substr(20, 8));
			Real   y         = RealParser(line.substr(28, 8));
			Real   z         = RealParser(line.substr(36, 8));

			if(i_molec != prev_molec_id % 100000)
				checkIfNewMolecule(i_molec, molec_name, atom_list, n_atoms, topol_info, molecs, prev_molec_name, prev_diff_id, prev_molec_id);

			const auto& molec_atoms= topol_info.atom_type_name_charge_mass[prev_diff_id];
			int lookup_idx=
#ifdef USE_VECTOR_TOPOLOGY
				n_atoms;
#else
				n_atoms + 1;
#endif

			const auto& [type, name, q, mass]= molec_atoms.at(lookup_idx);
			const auto& [e,s]= topol_info.type_LJparam.at(type);
			int Z= topol_info.type_Z.at(type);
			return Atom(Vector(x*10, y*10, z*10), i_atom, mass, q, e, s, Z, type);
		}

		static Vector readBounds(const string& line) {
			Real x, y, z;
			stringstream ss(line);
			ss >> x >> y >> z;
			return Vector(x*10, y*10, z*10);
		}

	public:
		GromacsCoordinateReader()= default;

		bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
			ifstream f(filename);
			if(!f.is_open()) { cout << "Error: Failed to open file " << filename << endl; return false; }
			if(!molecs) return false;

			string line;
			getline(f, line); // title
			getline(f, line);
			int natoms= stoi(line);

#ifdef DIFF_RESID_MOLECULE
			MoleculeCursor cur;
			Atom* atom_list= nullptr;
			int n_atoms= 0;

			for(int i= 0; i < natoms; i++) {
				getline(f, line);
				Atom a= readAtomByCount(line, topol_info, cur, atom_list, n_atoms, i);
				atom_list[n_atoms++]= a;
				if(n_atoms == cur.expected_atoms) {
					createNewMolecule(cur.current_name, cur.molecule_index + 1, molecs, atom_list, n_atoms);
					advanceMoleculeCursor(cur, topol_info);
					atom_list= nullptr;
					n_atoms= 0;
				}
			}
#else
			string prev_molec_name= "ERRORMOLECULE";
			int prev_diff_id= -1, prev_molec_id= -1;
			Atom* atom_list= nullptr;
			int n_atoms= 0;

			for(int i= 0; i < natoms; i++) {
				getline(f, line);
				Atom a= readAtom(line, topol_info, molecs, prev_molec_name, prev_diff_id, prev_molec_id, atom_list, n_atoms, i);
				atom_list[n_atoms++]= a;
			}
			if(atom_list != nullptr)
				createNewMolecule(prev_molec_name, prev_molec_id, molecs, atom_list, n_atoms);
#endif

			getline(f, line);
			bounds= readBounds(line);
			return true;
		}
};

#endif // GROMACS_READERS_HPP
