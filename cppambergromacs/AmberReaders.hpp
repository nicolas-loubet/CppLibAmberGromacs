#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

/**
 * Version: July 2025
 * Author: Ezequiel Cuenca
 */

#include "ReaderInterfaces.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

class AmberTopologyReader: public TopologyReader {
	private:

		struct FieldFormat {
			char type;  // 'I', 'E', 'F', 'A'
			int  count; // fields per line
			int  width; // characters per field
		};

		/**
		 * Parses a %FORMAT(...) line into a FieldFormat descriptor.
		 */
		inline FieldFormat parse_format(const string& format_line) const {
			FieldFormat f{'?', 0, 0};
			// Find the content inside parentheses
			size_t lp= format_line.find('(');
			size_t rp= format_line.rfind(')');
			if(lp == string::npos || rp == string::npos || rp <= lp) return f;
			string inner= format_line.substr(lp+1, rp-lp-1);

			// Remove leading/trailing spaces
			size_t s= inner.find_first_not_of(" \t");
			if(s == string::npos) return f;
			inner= inner.substr(s);

			// Identify type character (first alphabetic char)
			size_t ti= 0;
			while(ti < inner.size() && !isalpha((unsigned char)inner[ti])) ti++;
			if(ti >= inner.size()) return f;

			f.type= toupper((unsigned char)inner[ti]);
			f.count= (ti > 0) ? stoi(inner.substr(0, ti)) : 1;
			string after= inner.substr(ti+1);
			size_t dot= after.find('.');
			string width_str= (dot != string::npos) ? after.substr(0, dot) : after;
			// strip non-digits
			width_str.erase(remove_if(width_str.begin(), width_str.end(),
				[](char c){ return !isdigit((unsigned char)c); }), width_str.end());
			if(!width_str.empty()) f.width= stoi(width_str);

			return f;
		}

		// =====================================================================
		//  GENERIC SECTION READERS
		// =====================================================================

		/**
		 * Seeks past a %FLAG <flag_name> line and reads the following %FORMAT line.
		 * Leaves the file positioned at the first data line.
		 * Returns the parsed FieldFormat.
		 */
		inline FieldFormat seek_to_section(ifstream& file, int position, const string& flag_name) const {
			file.clear();
			file.seekg(position);
			string line;
			while(getline(file, line)) {
				if(line.find("%FLAG " + flag_name) != string::npos) {
					// Next line should be %FORMAT
					if(getline(file, line) && line.find("%FORMAT") != string::npos)
						return parse_format(line);
					break;
				}
			}
			return {'?', 0, 0};
		}

		/**
		 * Reads N integer values from the current file position using fixed-width fields.
		 */
		inline vector<int> read_n_ints(ifstream& file, int n, int width) const {
			vector<int> result;
			result.reserve(n);
			string line;
			while(getline(file, line) && (int)result.size() < n) {
				if(line.find("%FLAG") != string::npos) break;
				for(size_t i= 0; i + (size_t)width <= line.size() && (int)result.size() < n; i += width) {
					string chunk= ToolKit::strip(line.substr(i, width));
					if(!chunk.empty()) result.push_back(stoi(chunk));
				}
			}
			return result;
		}

		/**
		 * Reads N integer values using whitespace splitting (fallback for POINTERS).
		 */
		inline vector<int> read_n_ints_ws(ifstream& file, int n) const {
			vector<int> result;
			result.reserve(n);
			string line;
			while(getline(file, line) && (int)result.size() < n) {
				if(line.find("%FLAG") != string::npos) break;
				stringstream ss(line);
				int v;
				while(ss >> v && (int)result.size() < n) result.push_back(v);
			}
			return result;
		}

		/**
		 * Reads N real values from the current file position using fixed-width fields.
		 */
		inline vector<Real> read_n_reals(ifstream& file, int n, int width) const {
			vector<Real> result;
			result.reserve(n);
			string line;
			while(getline(file, line) && (int)result.size() < n) {
				if(line.find("%FLAG") != string::npos) break;
				for(size_t i= 0; i + (size_t)width <= line.size() && (int)result.size() < n; i += width) {
					string chunk= ToolKit::strip(line.substr(i, width));
					if(!chunk.empty()) result.push_back(RealParser(chunk));
				}
			}
			return result;
		}

		/**
		 * Reads N string values from the current file position using fixed-width fields.
		 */
		inline vector<string> read_n_strings(ifstream& file, int n, int width) const {
			vector<string> result;
			result.reserve(n);
			string line;
			while(getline(file, line) && (int)result.size() < n) {
				if(line.find("%FLAG") != string::npos) break;
				for(size_t i= 0; i + (size_t)width <= line.size() && (int)result.size() < n; i += width)
					result.push_back(ToolKit::strip(line.substr(i, width)));
			}
			return result;
		}

		// =====================================================================
		//  FLAG POSITION MAP
		// =====================================================================

		/**
		 * Scans the file and records the byte offset of each %FLAG line.
		 */
		inline map<string, streampos> flag_position(ifstream& file) const {
			map<string, streampos> positions;
			file.clear();
			file.seekg(0);
			string line;
			while(getline(file, line)) {
				if(line.find("%FLAG") != string::npos) {
					string name= ToolKit::strip(line.substr(5));
					streampos pos= file.tellg();
					positions[name]= (streampos)((long long)pos - (long long)line.size() - 1);
					if(positions[name] < (streampos)0) positions[name]= 0;
				}
			}
			return positions;
		}

		// =====================================================================
		//  SECTION READERS (all use seek_to_section + generic readers)
		// =====================================================================

		inline map<string,int> read_pointers(ifstream& file, streampos position) const {
			static const string names[32]= {
				"NATOM","NTYPES","NBONH","MBONA","NTHETH","MTHETA","NPHIH","MPHIA","NHPARM","NPARM",
				"NNB","NRES","NBONA","NTHETA","NPHIA","NUMBND","NUMANG","NPTRA","NATYP","NPHB",
				"IFPERT","NBPER","NGPER","NDPER","MBPER","MGPER","MDPER","IFBOX","NMXRS","IFCAP",
				"NUMEXTRA","NCOPY"
			};
			FieldFormat fmt= seek_to_section(file, (int)position, "POINTERS");
			vector<int> vals= (fmt.width > 0) ? read_n_ints(file, 32, fmt.width) : read_n_ints_ws(file, 32);

			map<string,int> ptrs;
			for(int i= 0; i < (int)vals.size() && i < 32; i++)
				ptrs[names[i]]= vals[i];
			return ptrs;
		}

		inline vector<string> read_atom_name(ifstream& file, int natom, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "ATOM_NAME");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse ATOM_NAME format");
			return read_n_strings(file, natom, fmt.width);
		}

		inline vector<tuple<string,Real>> read_charge(ifstream& file, int natom, const vector<string>& atom_names, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "CHARGE");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse CHARGE format");
			vector<Real> raw= read_n_reals(file, natom, fmt.width);

			vector<tuple<string,Real>> result;
			result.reserve(natom);
			for(int i= 0; i < natom; i++)
				result.emplace_back(atom_names[i], raw[i] / 18.2223);
			return result;
		}

		inline vector<Real> read_mass(ifstream& file, int natom, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "MASS");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse MASS format");
			return read_n_reals(file, natom, fmt.width);
		}

		inline vector<int> read_atomic_number(ifstream& file, int natom, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "ATOMIC_NUMBER");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse ATOMIC_NUMBER format");
			return read_n_ints(file, natom, fmt.width);
		}

		inline vector<int> read_ati(ifstream& file, int natom, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "ATOM_TYPE_INDEX");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse ATOM_TYPE_INDEX format");
			return read_n_ints(file, natom, fmt.width);
		}

		inline vector<string> read_amber_atom_type(ifstream& file, int natom, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "AMBER_ATOM_TYPE");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse AMBER_ATOM_TYPE format");
			return read_n_strings(file, natom, fmt.width);
		}

		/**
		 * Builds the ati_to_amber_type lookup: index = atom_type_index (1-based) -> amber type string.
		 */
		inline vector<string> build_ati_to_amber_type(int ntypes, int natom, const vector<int>& ati, const vector<string>& atype) const {
			vector<string> result(ntypes + 1);
			for(int i= 0; i < natom; i++)
				result[ati[i]]= atype[i];
			return result;
		}

		inline map<string,string> read_name_type(int natom, const vector<string>& atom_names, const vector<int>& ati, const vector<string>& ati_to_amber_type) const {
			map<string,string> result;
			for(int i= 0; i < natom; i++)
				result[atom_names[i]]= ati_to_amber_type[ati[i]];
			return result;
		}

		/**
		 * Reads SOLVENT_POINTERS section.
		 * Returns {IPTRES, NSPM, NSOLV} where:
		 *   IPTRES = last residue of solute
		 *   NSPM   = total number of molecules
		 *   NSOLV  = index of first solvent molecule (1-based)
		 */
		inline vector<int> read_solvent_pointers(ifstream& file, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "SOLVENT_POINTERS");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse SOLVENT_POINTERS format");
			return read_n_ints(file, 3, fmt.width);
		}

		/**
		 * Reads ATOMS_PER_MOLECULE - must be called after SOLVENT_POINTERS so we know nspm (= total number of molecules = NSPM from SOLVENT_POINTERS).
		 */
		inline vector<int> read_atoms_per_molecule(ifstream& file, int nspm, streampos position) const {
			FieldFormat fmt= seek_to_section(file, (int)position, "ATOMS_PER_MOLECULE");
			if(fmt.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse ATOMS_PER_MOLECULE format");
			return read_n_ints(file, nspm, fmt.width);
		}

		/**
		 * Reads LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF, converts A/B to epsilon/sigma
		 * Returns a map keyed by (amber_type_i, amber_type_j).
		 */
		inline map<pair<string,string>,pair<Real,Real>> read_lj(ifstream& file, int ntypes, const vector<string>& ati_to_amber_type, streampos position_acoef, streampos position_bcoef) const {
			const int n_pairs= (ntypes*(ntypes+1))/2;

			FieldFormat fmt_a= seek_to_section(file, (int)position_acoef, "LENNARD_JONES_ACOEF");
			if(fmt_a.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse LENNARD_JONES_ACOEF format");
			vector<Real> acoef= read_n_reals(file, n_pairs, fmt_a.width);

			FieldFormat fmt_b= seek_to_section(file, (int)position_bcoef, "LENNARD_JONES_BCOEF");
			if(fmt_b.width <= 0) throw runtime_error("AmberTopologyReader: cannot parse LENNARD_JONES_BCOEF format");
			vector<Real> bcoef= read_n_reals(file, n_pairs, fmt_b.width);

			// The flat array stores the lower triangle. Flat index for pair (i,k) with k<=i:  i*(i-1)/2 + (k-1)  (0-based)
			map<pair<string,string>,pair<Real,Real>> result;
			for(int i= 1; i <= ntypes; i++) {
				for(int k= 1; k <= ntypes; k++) {
					int lo= min(i,k), hi= max(i,k);
					int idx= hi*(hi-1)/2 + (lo-1);
					Real a= acoef[idx];
					Real b= bcoef[idx];

					Real sigma= 0.0, epsilon= 0.0;
					if(a != 0.0 && b != 0.0) {
						sigma  = pow(a / b, 1.0 / 6.0);
						epsilon= (b*b) / (4.0*a) * 4.184; // kcal -> kJ/mol
					}
					// Guard against any remaining nan/inf (e.g. B=0 but A!=0)
					if(!isfinite(sigma))   sigma=   0.0;
					if(!isfinite(epsilon)) epsilon= 0.0;

					auto key= make_pair(ati_to_amber_type[i], ati_to_amber_type[k]);
					result[key]= {epsilon, sigma};
				}
			}
			return result;
		}

		inline map<string,pair<Real,Real>> read_lj_diagonal(const map<pair<string,string>,pair<Real,Real>>& lj_map) const {
			map<string,pair<Real,Real>> result;
			for(const auto& [key, val]: lj_map)
				if(key.first == key.second)
					result[key.first]= val;
			return result;
		}

		inline map<string,int> type_atomic_number(int natom, const vector<int>& ati, const vector<int>& atomic_number, const vector<string>& ati_to_amber_type) const {
			map<string,int> result;
			for(int i= 0; i < natom; i++)
				result[ati_to_amber_type[ati[i]]]= atomic_number[i];
			return result;
		}

	public:
		/**
		 * Reads the AMBER topology file (.parm7 / .prmtop) and returns a TopolInfo object.
		 */
		TopolInfo readTopology(const string& filename) const override {
			TopolInfo topology;
			ifstream file(filename);
			if(!file.is_open()) throw runtime_error("AmberTopologyReader: cannot open file " + filename);

			map<string,streampos> pos= flag_position(file);

			map<string,int> ptrs= read_pointers(file, pos.at("POINTERS"));
			const int natom  = ptrs["NATOM"];
			const int ntypes = ptrs["NTYPES"];

			vector<string> atom_names         = read_atom_name       (file, natom, pos.at("ATOM_NAME"));
			vector<string> atype_per_atom     = read_amber_atom_type (file, natom, pos.at("AMBER_ATOM_TYPE"));
			vector<int>    ati                = read_ati             (file, natom, pos.at("ATOM_TYPE_INDEX"));
			vector<int>    atomic_number      = read_atomic_number   (file, natom, pos.at("ATOMIC_NUMBER"));
			vector<Real>   mass               = read_mass            (file, natom, pos.at("MASS"));
			vector<tuple<string,Real>> charges= read_charge(file, natom, atom_names, pos.at("CHARGE"));

			vector<string> ati_to_amber_type= build_ati_to_amber_type(ntypes, natom, ati, atype_per_atom);

			vector<int> solvent_ptrs = read_solvent_pointers(file, pos.at("SOLVENT_POINTERS"));
			const int nspm           = solvent_ptrs[1]; // total molecules
			const int nsolv_first    = solvent_ptrs[2]; // 1-based index of first solvent molecule
			vector<int> atoms_per_molecule= read_atoms_per_molecule(file, nspm, pos.at("ATOMS_PER_MOLECULE"));

			map<pair<string,string>,pair<Real,Real>> lj_all= read_lj(file, ntypes, ati_to_amber_type, pos.at("LENNARD_JONES_ACOEF"), pos.at("LENNARD_JONES_BCOEF"));
			map<string,pair<Real,Real>> lj_diagonal= read_lj_diagonal(lj_all);

			map<string,int>    type_Z    = type_atomic_number(natom, ati, atomic_number, ati_to_amber_type);
			map<string,string> name_type = read_name_type(natom, atom_names, ati, ati_to_amber_type);

			topology.num_molecules         = nspm;
			topology.num_solutes           = nsolv_first - 1;
			topology.num_solvents          = nspm - topology.num_solutes;
			topology.total_number_of_atoms = natom;
			topology.type_Z                = type_Z;
			topology.name_type             = name_type;
			topology.type_LJparam          = lj_diagonal;
			topology.special_interaction   = lj_all;

			int atom_idx= 0;

#ifdef USE_VECTOR_TOPOLOGY
			topology.atom_type_name_charge_mass.reserve(nspm);
			for(int mol= 0; mol < nspm; mol++) {
				vector<tuple<string,string,Real,Real>> mol_atoms;
				mol_atoms.reserve(atoms_per_molecule[mol]);
				for(int a= 0; a < atoms_per_molecule[mol]; a++, atom_idx++)
					mol_atoms.emplace_back(
						ati_to_amber_type[ati[atom_idx]],
						get<0>(charges[atom_idx]),
						get<1>(charges[atom_idx]),
						mass[atom_idx]
					);
				topology.atom_type_name_charge_mass.push_back(move(mol_atoms));
			}
#else
			for(int mol= 0; mol < nspm; mol++) {
				map<int,tuple<string,string,Real,Real>> mol_atoms;
				for(int a= 0; a < atoms_per_molecule[mol]; a++, atom_idx++)
					mol_atoms[a]= make_tuple(
						ati_to_amber_type[ati[atom_idx]],
						get<0>(charges[atom_idx]),
						get<1>(charges[atom_idx]),
						mass[atom_idx]
					);
				topology.atom_type_name_charge_mass.push_back(move(mol_atoms));
			}
#endif
			return topology;
		}
};

// =============================================================================
//  COORDINATE READER
// =============================================================================

class AmberCoordinateReader: public CoordinateReader {
	public:
		AmberCoordinateReader()= default;

#ifndef USE_FAST_READER
		/**
		 * Reads a PDB coordinate file, handling molecules split across multiple residues.
		 * Each TER record advances the molecule counter.
		 */
		bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
			ifstream f(filename);
			if(!f.is_open()) { cout << "Failed to open file " << filename << endl; return false; }
			if(!molecs) return false;

			string line;
			getline(f, line);
			if(line.rfind("CRYST1", 0) == 0)
				bounds= Vector(RealParser(line.substr(6,9)), RealParser(line.substr(15,9)), RealParser(line.substr(24,9)));
			f.clear();
			f.seekg(0);

			vector<Vector> coords;
			vector<string> mol_resname(topol_info.num_molecules);
			map<int, vector<int>> mol_atom_indices;
			int mol_idx= 0;

			while(getline(f, line)) {
				if(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) {
					mol_resname[mol_idx]= ToolKit::strip(line.substr(16, 4));
					mol_atom_indices[mol_idx].push_back((int)coords.size());
					coords.emplace_back(
						RealParser(ToolKit::strip(line.substr(30,8))),
						RealParser(ToolKit::strip(line.substr(38,8))),
						RealParser(ToolKit::strip(line.substr(46,8)))
					);
				}
				if(line.rfind("TER", 0) == 0) mol_idx++;
			}
			f.close();

			int global_atom_idx= 0;
			int out_mol_idx= 0;
			for(const auto& [mol_id, atom_indices]: mol_atom_indices) {
				const auto& atom_data= topol_info.atom_type_name_charge_mass[mol_id];
				int n= (int)atom_indices.size();
				Atom* atoms= new Atom[n];
				for(int j= 0; j < n; j++) {
					const auto& [type, name, charge, mass]= atom_data.at(j);
					const auto& [epsilon, sigma]= topol_info.type_LJparam.at(type);
					int Z= topol_info.type_Z.at(type);
					atoms[j]= Atom(coords[atom_indices[j]], global_atom_idx+1, mass, charge, epsilon, sigma, Z, type);
					global_atom_idx++;
				}
				if(mol_resname[out_mol_idx] == "WAT") molecs[out_mol_idx]= new Water(mol_id+1, atoms, n);
				else                                  molecs[out_mol_idx]= new Molecule(mol_id+1, atoms, n);
				out_mol_idx++;
			}
			return true;
		}

#else // USE_FAST_READER

		/**
		 * Assumes the PDB is perfectly ordered and each ATOM...TER block maps
		 * exactly to one molecule in topology order.
		 */
		bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
			ifstream f(filename);
			if(!f.is_open()) { cout << "Failed to open file " << filename << endl; return false; }

			string line;
			getline(f, line);
			if(line.rfind("CRYST1", 0) == 0)
				bounds= Vector(RealParser(line.substr(6,9)), RealParser(line.substr(15,9)), RealParser(line.substr(24,9)));

			for(int i= 0; i < topol_info.num_molecules; i++) {
				const auto& atom_data= topol_info.atom_type_name_charge_mass[i];
				Atom* atoms= new Atom[(int)atom_data.size()];
				int n= 0;
				string resname;

				while(getline(f, line)) {
					if(line.rfind("TER", 0) == 0) break;
					if(line.rfind("ATOM  ", 0) != 0 && line.rfind("HETATM", 0) != 0) continue;

					resname= ToolKit::strip(line.substr(16, 4));
					Real x= RealParser(ToolKit::strip(line.substr(30,8)));
					Real y= RealParser(ToolKit::strip(line.substr(38,8)));
					Real z= RealParser(ToolKit::strip(line.substr(46,8)));
					const auto& [type, name, charge, mass]= atom_data.at(n);
					const auto& [epsilon, sigma]= topol_info.type_LJparam.at(type);
					int Z= topol_info.type_Z.at(type);
					atoms[n]= Atom(Vector(x,y,z), n+1, mass, charge, epsilon, sigma, Z, type);
					n++;
				}
				if(resname == "WAT")  molecs[i]= new Water(i+1, atoms, n);
				else                  molecs[i]= new Molecule(i+1, atoms, n);
			}
			f.close();
			return true;
		}

#endif // USE_FAST_READER
};

#endif // AMBER_READERS_HPP
