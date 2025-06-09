#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

#include "ReaderFactory.hpp"
#include <string>
#include "amber_parser.hpp"

class AmberTopologyReader : public ReaderFactory::TopologyReader {
    public:
        Configuration::TopolInfo readTopology(const std::string& filename) const override {
            // Implementación para leer prmtop de AMBER-EZEQUIEL
            Configuration::TopolInfo topology;

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
                 //std::cout << "Clave: " << par.first << ", Atomos: " << std::get<0>(par.second) << " numero de moleculas: " << std::get<1>(par.second)<<std::endl;
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

class AmberCoordinateReader : public ReaderFactory::CoordinateReader {
    private:
        std::string filename;
    public:
        AmberCoordinateReader(const std::string& file): filename(file) {}
        bool readCoordinates(Molecule** molecules, int num_molecules, const Configuration::TopolInfo& topol_info) const override {
            // Implementación para leer mdcrd de AMBER
            return true;
        }
};

#endif