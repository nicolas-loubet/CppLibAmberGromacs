#ifndef TOPOLOGYREADER_HPP
#define TOPOLOGYREADER_HPP

/**
 * Version: May 2025
 * Author: Nicolás Loubet
 */

#include "Configuration.hpp"
#include <fstream>
#include <sstream>

namespace TopolReader {
	/**
	 * Removes all spaces from the string
	 * @param _s The string to strip
	 * @return The stripped string
	 */
	string strip(string _s)
	{
		_s.erase(remove(_s.begin(), _s.end(), ' '), _s.end());
		return(_s);
	}

	/**
	 * Given a file, find the position of the flag in the file.
	 * The returned map contains the flag name as the key and the position as the value.
	 * The position is the byte position of the first character after the flag.
	 * @param file The file to search
	 * @return A map with the position of the flag.
	 */
	map<string, int> flag_position_AMBER(ifstream &file)
	{
		string line= "";
		map<string,int> flag;
		while(getline(file, line)) {
			if(line.find("%FLAG") != string::npos) {
				line=strip(line.substr(5, 80));
				flag[line]= file.tellg();
				flag[line]-= 81;
			}
		}
		return(flag);
	}

	string readMolecType(ifstream &file) {
		string line= "";
		while(getline(file, line))
			if(line[0] != ';')
				return strip(line.substr(0,8));
		return "ERROR";
	}

	/**
	 * Given a file, find the position of the flag in the file.
	 * The returned map contains the flag name as the key and the position as the value.
	 * The position is the byte position of the first character after the flag.
	 * @param file The file to search
	 * @return A map with the position of the flag.
	 */
	map<string,int> flag_position_GROMACS(ifstream &file) {
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

	map<string,int> readSystemFlag(ifstream &file, int position) {
		map<string,int> molecules;
		string line= "";
		file.clear();
		file.seekg(position);

		while(getline(file, line)) {
			if(line[0] == ';') continue;
			string molec_name= strip(line.substr(0,8));
			int num_molec= stoi(line.substr(9,line.length()-9));
			molecules[molec_name]= num_molec;
		}
		return molecules;
	}

	int sumMoleculesInMap(map<string,int> molecules) {
		int sum= 0;
		for(auto it= molecules.begin(); it != molecules.end(); it++)
			sum+= it->second;
		return sum;
	}

	int sumMoleculesConsideredSolvent(map<string,int> molecules) {
		int sum= 0;
		for(auto it= molecules.begin(); it != molecules.end(); it++)
			if((it->first.find("SOL") != string::npos) || (it->first.find("WAT") != string::npos))
				sum+= it->second;
		return sum;
	}

	map<int,tuple<string,string,float,float>> readAtomsFlags(ifstream &file, int position) {
		map<int,tuple<string,string,float,float>> atoms;
		string line= "";
		file.clear();
		file.seekg(position);

		while(getline(file, line)) {
			if(line[0] == ';') continue;
			if(strip(line) == "") break;

			stringstream ss(line.substr(0,line.find(";")));

			int nr,resi,cgnr;
			string type,res,atom;
			float charge,mass;

			ss >> nr >> type >> resi >> res >> atom >> cgnr >> charge >> mass;
			atoms[nr]= make_tuple(type, atom, charge, mass);
		}
		return atoms;
	}

	static Configuration::TopolInfo readTopolInfoOfAmber(std::string filename) {
		// Implementación para leer archivo de topología de AMBER
		return Configuration::TopolInfo();
	}

	static Configuration::TopolInfo readTopolInfoOfGromacs(std::string filename) {
		Configuration::TopolInfo ti= Configuration::TopolInfo();

		string line;
		ifstream f(filename);

		if(!f.is_open()) {
			std::cerr << "Topology not found" << std::endl;
			return ti;
		}

		map<string,int> flags= flag_position_GROMACS(f);
		ti.number_of_each_different_molecule= readSystemFlag(f, flags["[ molecules ]"]);
		ti.num_molecules= sumMoleculesInMap(ti.number_of_each_different_molecule);
		ti.num_solvents= sumMoleculesConsideredSolvent(ti.number_of_each_different_molecule);
		ti.num_solutes= ti.num_molecules-ti.num_solvents;
		
		for(auto it_molec= ti.number_of_each_different_molecule.begin(); it_molec != ti.number_of_each_different_molecule.end(); it_molec++) {
			map<int,tuple<string,string,float,float>> atoms= readAtomsFlags(f, flags["[ atoms ]_"+it_molec->first]);
			ti.number_of_atoms_per_different_molecule[it_molec->first]= atoms.size();
			ti.atom_type_name_charge_mass.push_back(atoms);
			for(auto it_atom= atoms.begin(); it_atom != atoms.end(); it_atom++)
				ti.name_type[it_molec->first+":"+get<1>(it_atom->second)]= get<0>(it_atom->second);
		}

		//Falta leer 	map<string,pair<float,float>> type_LJparam
		//y 			map<pair<string,string>,pair<float,float>> special_interaction;

		return ti;
	}
}

#endif