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
        inline map<string, int> flag_position(ifstream &file) const
        {
            string line="";
            map<string, int> flag;
            bool inFlag=false;
            while (getline(file, line)) {
                
                if (line.find("%FLAG") != string::npos) {
                    line=ToolKit::strip(line.substr(5, 80));
                    inFlag = true;
                    flag[line]=file.tellg();
                    flag[line]-=200;
                    continue;
                }
            }
            return(flag);
        }

        inline map<string,int> read_pointers(ifstream& file_prmtop,int position) const
        {
            string list_of_names[32]= {
                "NATOM","NTYPES","NBONH","MBONA","NTHETH","MTHETA","NPHIH","MPHIA","NHPARM","NPARM",
                "NNB","NRES","NBONA","NTHETA","NPHIA","NUMBND","NUMANG","NPTRA","NATYP","NPHB",
                "IFPERT","NBPER","NGPER","NDPER","MBPER","MGPER","MDPER","IFBOX","NMXRS","IFCAP",
                "NUMEXTRA","NCOPY"
            };

            string line;
            map<string,int> dict_pointers;
            bool inPointers = false;
            int _i = 0;
            file_prmtop.clear();
            file_prmtop.seekg(position);
            

            while (getline(file_prmtop, line)) {
                if (line.find("%FLAG POINTERS") != string::npos) {
                    inPointers = true;
                    continue;
                }
                
                if (inPointers) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }
                    // Read and split variable values
                    stringstream ss(line);
                    string this_value;
                    while (ss >> this_value) {
                        dict_pointers[list_of_names[_i]] = stoi(this_value);
                        _i++;
                    }
                    if (_i>=31) {
                        // Finished reading the POINTERS section
                        break;
                    }
                }
            }
            return dict_pointers;
        }

        inline vector<string> read_atom_name(ifstream &file, map<string,int>& dict_pointers,int position) const
        {
            string line;
            bool inFlag = false;
            file.clear();
            file.seekg(position);
            vector<string> atoms(dict_pointers["NATOM"]);

            int _j = 0;
            while (getline(file, line)) {
                if (line.find("%FLAG ATOM_NAME") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }
                    
                    // Read and split variable values
                    stringstream ss(line);

                    for (size_t i = 0; i + 3 < line.length(); i += 4) {
                        atoms[_j]=ToolKit::strip(line.substr(i, 4));
                        ++_j;
                    }
                    if (_j>=dict_pointers["NATOM"]) {
                        // Finished reading the ATOM_NAME section
                        break;
                    }
                }
            }
            return(atoms);
        }

        inline vector<tuple<string,float>> read_charge(ifstream &file, map<string,int>& dict_pointers, vector<string> atom_names,int position) const
        {
            string line = "";
            bool inFlag = false;
            vector<tuple<string,float>> charges;
            int _j = 0;
            file.clear();
            file.seekg(position);
            while (getline(file, line)) {
                if (line.find("%FLAG CHARGE") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }
                    stringstream ss(line);

                    for (size_t i = 0; i + 15 < line.length(); i += 16) {
                        charges.push_back(make_tuple(atom_names[_j], stof(line.substr(i, 16))/18.2223));
                        ++_j;
                    }
                    
                    if (_j>=dict_pointers["NATOM"]) {
                        // Finished reading the CHARGE section
                        break;
                    }
                }
            }
            
            return(charges);
        }

        inline vector<int> read_solvent_pointers(ifstream &file, map<string,int>& dict_pointers,int position) const
        {
            string line = "";
            int nspm= 0;
            vector<int> number_solute_solvent(3);
            file.clear();
            file.seekg(position);
            while (getline(file, line)) {
                if (line.find("%FLAG SOLVENT_POINTERS") != string::npos) {
                    getline(file, line);
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        getline(file, line);
                    }
                    
                    number_solute_solvent[0]=stoi(line.substr(0, 8));
                    number_solute_solvent[1]=stoi(line.substr(8, 8));
                    number_solute_solvent[2]=stoi(line.substr(16, 8));
                    break;
                }
            }
            return(number_solute_solvent);
        }

        inline vector<string> atom_type_index_to_amber_type(ifstream &file, map<string,int>& dict_pointers , map<string, int> position) const
        {
            vector<string> ati_to_amber_type(dict_pointers["NTYPES"]+1);

            string line="";
            bool inFlag = false;
            vector<int> ati(dict_pointers["NATOM"]);
            vector<string> atype(dict_pointers["NATOM"]);

            int _j=0;
            file.clear();
            file.seekg(position["AMBER_ATOM_TYPE"]);
            while (getline(file, line)) {
                
                if (line.find("%FLAG AMBER_ATOM_TYPE") != string::npos) {
                    inFlag = true;
                    break;
                }
            }

            if (inFlag) {
                while (getline(file, line)) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Saltamos la línea del formato
                        continue;
                    }
                    stringstream ss(line);
                    
                    
                    for (size_t i = 0; i + 3 < line.length(); i += 4) {
                        atype[_j]=ToolKit::strip(line.substr(i, 4));
                        ++_j;
                    }
                    if (_j>=dict_pointers["NATOM"]) {
                        // Terminó la sección de ATOM_TYPE_INDEX
                        break;
                    }
                }
            }
            _j=0;
            file.clear();
            file.seekg(position["ATOM_TYPE_INDEX"]);
            while (getline(file, line)) {
                
                if (line.find("%FLAG ATOM_TYPE_INDEX") != string::npos) {
                    inFlag = true;
                    break;
                }
            }

            if (inFlag) {
                while (getline(file, line)) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Saltamos la línea del formato
                        continue;
                    }
                    stringstream ss(line);
                    
                    
                    for (size_t i = 0; i + 7 < line.length(); i += 8) {
                        ati_to_amber_type[stoi(ToolKit::strip(line.substr(i, 8)))]=atype[_j];
                        ++_j;
                    }
                    if (_j>=dict_pointers["NATOM"]) {
                        // Terminó la sección de ATOM_TYPE_INDEX
                        break;
                    }
                }
            }
            return(ati_to_amber_type);
        }

        inline map<string,string> read_name_type(map<string,int>& dict_pointers, vector<string> atom_names, vector<int> atom_type_index, vector<string> ati_to_amber_type) const
        {
            map<string,string> name_type;
            for (int i=0;i<dict_pointers["NATOM"];++i)
            {
                name_type[atom_names[i]]=ati_to_amber_type[atom_type_index[i]];
            }
            return(name_type);
        }


        inline map<string,tuple<int,int>> read_atoms_per_different_molecule(ifstream &file, map<string, int> position) const
        {
            string line = "";
            int nspm= 0;
            file.clear();
            file.seekg(position["SOLVENT_POINTERS"]);
            while (getline(file, line)) {
                if (line.find("%FLAG SOLVENT_POINTERS") != string::npos) {
                    getline(file, line);
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        getline(file, line);
                    }
                    nspm= stoi(line.substr(8, 8));
                    
                    break;
                }
            }
            int _j=0;
            vector<int> atoms_per_molecule(nspm);
            bool inFlag = false;

            while (getline(file, line)) {
                if (line.find("%FLAG ATOMS_PER_MOLECULE") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }

                    for (size_t i = 0; i + 7 < line.length(); i += 8) {
                        atoms_per_molecule[_j]=stoi(line.substr(i, 8));
                        ++_j;
                    }
                    if (_j>=nspm) {
                        // Finished reading the SOLVENT POINTERS section
                        break;
                    }
                }
            }

            file.clear();
            file.seekg(position["RESIDUE_LABEL"]);
            int _k=0;
            vector<string> label_per_molecule(nspm);
            map<string,tuple<int,int>> label_and_number_atom;

            inFlag = false;
            
            
            while (getline(file, line)) {
                if (line.find("%FLAG RESIDUE_LABEL") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }

                    for (size_t i = 0; i + 3 < line.length(); i += 4) {
                        label_per_molecule[_k]=ToolKit::strip(line.substr(i, 4));
                        string name = label_per_molecule[_k] + to_string(atoms_per_molecule[_k]);
                        label_and_number_atom[name]=make_tuple(atoms_per_molecule[_k], 0);
                        ++_k;
                    } 
                    if (_k>=nspm) {
                        // Finished reading the SOLVENT POINTERS section
                        break;
                    }
                }
            }

            file.clear();
            file.seekg(position["RESIDUE_LABEL"]);
            int _l=0;
            inFlag = false;
            
            
            while (getline(file, line)) {
                if (line.find("%FLAG RESIDUE_LABEL") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }

                    for (size_t i = 0; i + 3 < line.length(); i += 4) {
                        label_per_molecule[_l]=ToolKit::strip(line.substr(i, 4));
                        stringstream ss;
                        ss << atoms_per_molecule[_l];
                        string number = ss.str();
                        string name = label_per_molecule[_l] + number;
                        label_and_number_atom[name]=make_tuple(atoms_per_molecule[_l], (get<1>(label_and_number_atom[name]))+1);
                        ++_l;
                    } 
                    if (_l>=nspm) {
                        // Finished reading the SOLVENT POINTERS section
                        break;
                    }
                }
            }
            return(label_and_number_atom);
        }


        inline vector<int> read_atoms_per_molecule(ifstream &file, map<string,int>& dict_pointers,int position) const
        {
            string line = "";
            int nspm= 0;
            file.clear();
            file.seekg(position);
            while (getline(file, line)) {
                if (line.find("%FLAG SOLVENT_POINTERS") != string::npos) {
                    getline(file, line);
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        getline(file, line);
                    }
                    nspm= stoi(line.substr(8, 8));
                    
                    break;
                }
            }

            int _j=0;
            vector<int> atoms_per_molecule(nspm);
            bool inFlag = false;

            while (getline(file, line)) {
                if (line.find("%FLAG ATOMS_PER_MOLECULE") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }

                    for (size_t i = 0; i + 7 < line.length(); i += 8) {
                        atoms_per_molecule[_j]=stoi(line.substr(i, 8));
                        ++_j;
                    }
                    if (_j>=nspm) {
                        // Finished reading the SOLVENT POINTERS section
                        break;
                    }
                }
            }
            return(atoms_per_molecule);
        }

        inline vector<int> read_atomic_number(ifstream &file, map<string,int>& dict_pointers,int position) const
        {
            string line="";
            bool inFlag = false;
            vector<int> atomic_number(dict_pointers["NATOM"]);
            int _j=0;
            file.clear();
            file.seekg(position);
            while (getline(file, line)) {
                
                if (line.find("%FLAG ATOMIC_NUMBER") != string::npos) {
                    inFlag = true;
                    break;
                }
            }

            if (inFlag) {
                while (getline(file, line)) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Saltamos la línea del formato
                        continue;
                    }
                    stringstream ss(line);
                    
                    
                    for (size_t i = 0; i + 7 < line.length(); i += 8) {
                        atomic_number[_j]=stoi(ToolKit::strip(line.substr(i, 8)));
                        ++_j;
                    }
                    if (_j>=dict_pointers["NATOM"]) {
                        // Terminó la sección de ATOMIC_NUMBER
                        break;
                    }
                }
            }

            return(atomic_number);
        }

        inline vector<int> read_ati(ifstream &file, map<string,int>& dict_pointers,int position) const
        {
            string line="";
            bool inFlag = false;
            vector<int> ati(dict_pointers["NATOM"]);
            int _j=0;
            file.clear();
            file.seekg(position);
            while (getline(file, line)) {
                
                if (line.find("%FLAG ATOM_TYPE_INDEX") != string::npos) {
                    inFlag = true;
                    break;
                }
            }

            if (inFlag) {
                while (getline(file, line)) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Saltamos la línea del formato
                        continue;
                    }
                    stringstream ss(line);
                    
                    for (size_t i = 0; i + 7 < line.length(); i += 8) {
                        ati[_j]=stoi(ToolKit::strip(line.substr(i, 8)));
                        ++_j;
                    }
                    if (_j>=dict_pointers["NATOM"]) {
                        // Terminó la sección de ATOM_TYPE_INDEX
                        break;
                    }
                }
            }

            return(ati);
        }

        inline vector<float> read_mass(ifstream &file, map<string,int>& dict_pointers,int position) const
        {
            string line="";
            bool inFlag = false;
            vector<float> mass(dict_pointers["NATOM"]); 
            int _j=0;
            file.clear();
            file.seekg(position);
            while (getline(file, line)) {
                
                if (line.find("%FLAG MASS") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Saltamos la línea del formato
                        continue;
                    }
                    stringstream ss(line);
                    
                    
                    for (size_t i = 0; i + 15 < line.length(); i += 16) {
                        mass[_j]=stof(ToolKit::strip(line.substr(i, 16)));
                        ++_j;
                    }
                    if (_j>=dict_pointers["NATOM"]) {
                        // Terminó la sección de ATOM_TYPE_INDEX
                        break;
                    }
                }
            }
            
            return(mass);
        }

        inline map<pair<string,string>,pair<float,float>> read_lj(ifstream &file, map<string,int>& dict_pointers,int position,vector<string> amber_type) const
        {
            string line="";
            bool inFlag = false;
            map<tuple<int,int>,tuple<float,float>> lj_coefficient;
            file.clear();
            file.seekg(position);
            int _j=0;
            int _i=1;//atomo 1
            int _k=1;//atomo 2

            const int cut_value= (dict_pointers["NTYPES"]*(dict_pointers["NTYPES"]+1))/2;
            
            while (getline(file, line)) {
                
                if (line.find("%FLAG LENNARD_JONES_ACOEF") != string::npos) {
                    inFlag = true;
                    
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        
                        continue;
                    }
                    
                    // Read and split variable values
                    stringstream ss(line);
                    
                    for (size_t i = 0; i + 15< line.length(); i += 16) {
                        float coef_a=stof(ToolKit::strip(line.substr(i, 16)));
                        
                        //lj_coefficient[{_i,_k}]={coef_a,0.0};
                        get<0>(lj_coefficient[{_i,_k}])=coef_a;
                        get<0>(lj_coefficient[{_k,_i}])=coef_a;
                        _k+=1;
                        if(_k>_i){_k=1; _i+=1;}
                        
                        _j+=1;
                    }
                    if (_j>=cut_value) {
                        // Finished reading the ATOM_NAME section
                        break;
                    }
                }
            }
            _j=0;
            _i=1;
            _k=1;
            while (getline(file, line)) {
                if (line.find("%FLAG LENNARD_JONES_BCOEF") != string::npos) {
                    inFlag = true;
                    continue;
                }

                if (inFlag) {
                    if (line.find("%FORMAT") != string::npos) {
                        // Skipping format line
                        continue;
                    }
                    
                    // Read and split variable values
                    stringstream ss(line);
                    
                    for (size_t i = 0; i + 15< line.length(); i += 16) {
                        float coef_b=stof(ToolKit::strip(line.substr(i, 16)));

                        get<1>(lj_coefficient[{_i,_k}])=coef_b;
                        get<1>(lj_coefficient[{_k,_i}])=coef_b;

                        _k+=1;
                        if(_k>_i){_k=1; _i+=1;}
                        
                        _j+=1;
                    }
                    if (_j>=(dict_pointers["NTYPES"]*(dict_pointers["NTYPES"]+1))/2) {
                        // Finished reading the ATOM_NAME section
                        break;
                    }
                }
            }

            float sigma=0;
            float epsilon=0;
            float a=0;
            float b=0;
            for(int _i=1;_i<=dict_pointers["NTYPES"];++_i)
                {
                for(int _k=1;_k<=dict_pointers["NTYPES"];++_k)
                    {
                        a=get<0>(lj_coefficient[{_i,_k}]);
                        b=get<1>(lj_coefficient[{_i,_k}]);
                        sigma=pow(a / b, 1.0 / 6.0);
                        epsilon=pow(b,2)/(4*a);
                        get<0>(lj_coefficient[{_i,_k}])=epsilon*4.184; // kcal -> kJ
                        get<1>(lj_coefficient[{_i,_k}])=sigma;
                    }
                }
            map<pair<string,string>,pair<float,float>> lj_coefficient_in_amber_type;
            for(int _i=1;_i<=dict_pointers["NTYPES"];++_i)
                {
                for(int _k=1;_k<=dict_pointers["NTYPES"];++_k)
                    {
                        pair<string,string> lj_type=make_pair(amber_type[_i],amber_type[_k]);
                        pair<float,float> lj_coef=make_pair(get<0>(lj_coefficient[{_i,_k}]),get<1>(lj_coefficient[{_i,_k}]));
                        lj_coefficient_in_amber_type[lj_type]=lj_coef;
                    }
                }
            return(lj_coefficient_in_amber_type);
        }

        inline map<string,pair<float,float>> read_lj_diagonal(map<pair<string,string>,pair<float,float>> lj_map) const
        {
            map<string,pair<float,float>> lj_coefficient;
            for (const auto& par : lj_map) 
            {
                if(get<0>(par.first)==get<1>(par.first))
                {
                    lj_coefficient[get<0>(par.first)]=par.second;
                }
            }
            return(lj_coefficient);
        }

        inline map<string,int> type_atomic_number(map<string,int>& dict_pointers,vector<int> atom_type_index, vector<int> atomic_number, vector<string> ati_to_amber_type) const
        {
            map<string,int> type_atomic_number_index;
            for (int i=0;i<dict_pointers["NATOM"];++i)
            {
                type_atomic_number_index[ati_to_amber_type[atom_type_index[i]]]=atomic_number[i];
            }
            return(type_atomic_number_index);
        }

    public:
        TopolInfo readTopology(const string& filename) const override {
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
            
            for (const auto& pair : atoms_per_diff_molecule) {
                topology.number_of_atoms_per_different_molecule[pair.first]=get<0>(pair.second);
                topology.number_of_each_different_molecule[pair.first]=get<1>(pair.second);
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
            return topology;
        }
};

class AmberCoordinateReader : public CoordinateReader {
    public:
        AmberCoordinateReader()= default;

        /**
         * Reads the coordinates file .pdb
         * @param filename The name of the coordinates file (single frame)
         * @param topol_info The topology information (use AmberTopologyReader)
         * @param molecs An empty array of molecule pointers
         * @return True if the coordinates were read successfully
         */
        bool readCoordinates_slow(const string& filename, const TopolInfo& topol_info, Molecule* molecs) const {
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

                molecs[molec_idx]= Molecule(mol_pair.first, atoms, num_atoms_per_molecule);
                molec_idx+=1;
            }
            return true;
        }

        bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs, Vector& bounds) const override {
            ifstream f(filename);
            if(!f.is_open()) {
                cout << "Failed to open file " << filename << endl;
                return false;
            }
            string line;
            int number_of_atoms=0;
            string is_water="";
            
            int water_type=1; // Default type is TIP3P water 
            for (const auto& par : topol_info.type_LJparam)  // search OW in type_LJparam only once per PDB file
            {
                if(par.first=="OW")
                {
                    float epsilon_water=get<0>(topol_info.type_LJparam.at("OW")); //obtain epsilon of water
                    int epsilon_value = static_cast<int>(epsilon_water*1000); //switch only accepts int or enum so epsilon is converted to int
                    switch(epsilon_value) //switch is faster thas nested ifs
                    {
                        case 635:
                        case 636: //TIP3P
                            water_type=1;
                        break;
                        case 774: //TIP4P/2005
                            water_type=2;
                        break;
                        case 650: //SPC/E
                            water_type=3;
                        break;
                        case 790: //TIP5P-2018
                            water_type=4;
                        break;
                        default:
                            throw std::runtime_error("Unsupported water model");
                    }
                    break;
                }
            }
            
            getline(f,line);
            if(line.rfind("CRYST1",0) == 0)
                bounds= Vector( stof(line.substr(6,9)) , stof(line.substr(15,9)) , stof(line.substr(24,9)) );
                
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
                        is_water = ToolKit::strip(line.substr(16, 4));
                        const auto& [type,name,charge,mass]= atom_data.at(number_of_atoms);
                        const auto& [epsilon,sigma]= topol_info.type_LJparam.at(type);
                        int Z=topol_info.type_Z.at(type);

                        atoms[number_of_atoms]= Atom(Vector(x,y,z), number_of_atoms+1, mass, charge, epsilon, sigma, Z);
                        number_of_atoms+=1;
                    }
                }
                if(is_water != "WAT")
                {
                    molecs[i] = new Molecule(i+1, atoms, number_of_atoms);
                    continue;
                }

                switch(water_type)
                {
                    case 1: //TIP3P
                        molecs[i] = new TIP3PWater(i+1, atoms);
                    break;
                    case 2: //TIP4P/2005
                        molecs[i] = new TIP4PWater(i+1, atoms);
                    break;
                    case 3: //SPC/E
                        molecs[i] = new SPCEWater(i+1, atoms);
                    break;
                    case 4: //TIP5P-2018
                        molecs[i] = new TIP5PWater(i+1, atoms);
                    break;
                }
            }
            f.close();

            return true;
        }
};

#endif
