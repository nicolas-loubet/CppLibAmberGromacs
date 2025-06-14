#ifndef AMBER_PARSER_HPP
#define AMBER_PARSER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <algorithm>
#include <ctime>
#include <cmath>
#include "General/ToolKit.hpp"

inline std::map<std::string,int> read_pointers(std::ifstream& file_prmtop,int position)
{
    std::string list_of_names[32]= {
        "NATOM","NTYPES","NBONH","MBONA","NTHETH","MTHETA","NPHIH","MPHIA","NHPARM","NPARM",
        "NNB","NRES","NBONA","NTHETA","NPHIA","NUMBND","NUMANG","NPTRA","NATYP","NPHB",
        "IFPERT","NBPER","NGPER","NDPER","MBPER","MGPER","MDPER","IFBOX","NMXRS","IFCAP",
        "NUMEXTRA","NCOPY"
    };

    std::string line;
    std::map<std::string,int> dict_pointers;
    bool inPointers = false;
    int _i = 0;
    file_prmtop.clear();
    file_prmtop.seekg(position);
    

    while (getline(file_prmtop, line)) {
        if (line.find("%FLAG POINTERS") != std::string::npos) {
            inPointers = true;
            continue;
        }
        
        if (inPointers) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }
            // Read and split variable values
            std::stringstream ss(line);
            std::string this_value;
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


inline std::map<std::string,std::string> read_name_type(std::map<std::string,int>& dict_pointers, std::vector<std::string> atom_names, std::vector<int> atom_type_index, std::vector<std::string> ati_to_amber_type)
{
    std::map<std::string,std::string> name_type;
    for (int i=0;i<dict_pointers["NATOM"];++i)
    {
        name_type[atom_names[i]]=ati_to_amber_type[atom_type_index[i]];
    }
    return(name_type);
}

inline std::vector<std::string> read_atom_name(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position)
{
    std::string line;
    bool inFlag = false;
    file.clear();
    file.seekg(position);
    std::vector<std::string> atoms(dict_pointers["NATOM"]);

    int _j = 0;
    while (getline(file, line)) {
        if (line.find("%FLAG ATOM_NAME") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }
            
            // Read and split variable values
            std::stringstream ss(line);

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


inline std::vector<std::tuple<std::string,float>> read_charge(std::ifstream &file, std::map<std::string,int>& dict_pointers, std::vector<std::string> atom_names,int position)
{
    std::string line = "";
    bool inFlag = false;
    std::vector<std::tuple<std::string,float>> charges;
    int _j = 0;
    file.clear();
    file.seekg(position);
    while (getline(file, line)) {
        if (line.find("%FLAG CHARGE") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }
            std::stringstream ss(line);

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
inline std::vector<int> read_solvent_pointers(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position)
{
    std::string line = "";
    int nspm= 0;
    std::vector<int> number_solute_solvent(3);
    file.clear();
    file.seekg(position);
    while (getline(file, line)) {
        if (line.find("%FLAG SOLVENT_POINTERS") != std::string::npos) {
            getline(file, line);
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                getline(file, line);
            }
            number_solute_solvent[0]=stoi(line.substr(0, 8));
            number_solute_solvent[1]=stoi(line.substr(8, 8));
            number_solute_solvent[2]=stoi(line.substr(16, 8));
            break;
        }
    }
    //std::cout << "IPTRES: "<< number_solute_solvent[0] << " NSPM: "<< number_solute_solvent[1] <<" NSPSOL: "<< number_solute_solvent[2] << std::endl; 
    //std::cout << "NSPM:" << nspm << std::endl; 
    return(number_solute_solvent);
}

inline std::map<std::string,std::tuple<int,int>> read_atoms_per_different_molecule(std::ifstream &file, std::map<std::string, int> position)
{
    std::string line = "";
    int nspm= 0;
    file.clear();
    file.seekg(position["SOLVENT_POINTERS"]);
    while (getline(file, line)) {
        if (line.find("%FLAG SOLVENT_POINTERS") != std::string::npos) {
            getline(file, line);
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                getline(file, line);
            }
            nspm= stoi(line.substr(8, 8));
            
            break;
        }
    }

    int _j=0;
    std::vector<int> atoms_per_molecule(nspm);
    bool inFlag = false;

    while (getline(file, line)) {
        if (line.find("%FLAG ATOMS_PER_MOLECULE") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
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
    std::vector<std::string> label_per_molecule(nspm);
    std::map<std::string,std::tuple<int,int>> label_and_number_atom;
    //std::map<std::string,int> label_and_number_atom;

    //std::cout << position["RESIDUE_LABEL"] << std::endl;    

    inFlag = false;
    
    
    while (getline(file, line)) {
        if (line.find("%FLAG RESIDUE_LABEL") != std::string::npos) {
            inFlag = true;
            //std::cout << "hola" << std::endl;

            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }

            for (size_t i = 0; i + 3 < line.length(); i += 4) {
                label_per_molecule[_k]=ToolKit::strip(line.substr(i, 4));
                std::stringstream ss;
                ss << atoms_per_molecule[_k];
                std::string number = ss.str();
                std::string name = label_per_molecule[_k] + number;
                //std::cout << name<< std::endl;
                //label_and_number_atom[name]=atoms_per_molecule[_k];
                label_and_number_atom[name]=std::make_tuple(atoms_per_molecule[_k], 0);
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
        if (line.find("%FLAG RESIDUE_LABEL") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }

            for (size_t i = 0; i + 3 < line.length(); i += 4) {
                label_per_molecule[_l]=ToolKit::strip(line.substr(i, 4));
                std::stringstream ss;
                ss << atoms_per_molecule[_l];
                std::string number = ss.str();
                std::string name = label_per_molecule[_l] + number;
                label_and_number_atom[name]=std::make_tuple(atoms_per_molecule[_l], (std::get<1>(label_and_number_atom[name]))+1);
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


inline std::vector<int> read_atoms_per_molecule(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position)
{
    std::string line = "";
    int nspm= 0;
    file.clear();
    file.seekg(position);
    while (getline(file, line)) {
        if (line.find("%FLAG SOLVENT_POINTERS") != std::string::npos) {
            getline(file, line);
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                getline(file, line);
            }
            nspm= stoi(line.substr(8, 8));
            
            break;
        }
    }
    //std::cout << "NSPM:" << nspm << std::endl; 

    int _j=0;
    std::vector<int> atoms_per_molecule(nspm);
    bool inFlag = false;

    while (getline(file, line)) {
        if (line.find("%FLAG ATOMS_PER_MOLECULE") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }

            for (size_t i = 0; i + 7 < line.length(); i += 8) {
                atoms_per_molecule[_j]=stoi(line.substr(i, 8));
                //std::cout<<"atoms_per_molecule[_j]: "<< atoms_per_molecule[_j] <<std::endl;
                ++_j;
            }
            if (_j>=nspm) {
                // Finished reading the SOLVENT POINTERS section
                break;
            }
        }
    }
    /*
    int atom_index=0;
    std::vector<int> atom_number_in_molec(dict_pointers["NATOM"]);
    
    for(int molecule_index=0; molecule_index<nspm;++molecule_index)
    {
        for(int _k=0; _k<atoms_per_molecule[molecule_index]; ++_k)
        {
            atom_number_in_molec[atom_index]= molecule_index;
            ++atom_index;
        }
    }
    return(atom_number_in_molec);*/
    return(atoms_per_molecule);
}

inline std::vector<int> read_atomic_number(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position)
{
    std::string line="";
    bool inFlag = false;
    std::vector<int> atomic_number(dict_pointers["NATOM"]);
    int _j=0;
    file.clear();
    file.seekg(position);
    while (getline(file, line)) {
        
        if (line.find("%FLAG ATOMIC_NUMBER") != std::string::npos) {
            inFlag = true;
            break;
            //continue;
        }
    }

    if (inFlag) {
        while (getline(file, line)) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Saltamos la línea del formato
                continue;
            }
            std::stringstream ss(line);
            
            
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
    //DESCOMENTAR ESTA
    /*
    for (size_t i = 0; i < dict_pointers["NATOM"]; ++i) {
        int _number = ati[i];
        std::cout << "Atom " << i << ": " << "atom_type_index" <<":"<< _number << std::endl; 
    }*/

    return(atomic_number);
}

inline std::vector<std::string> atom_type_index_to_amber_type(std::ifstream &file, std::map<std::string,int>& dict_pointers , std::map<std::string, int> position)
{
    std::vector<std::string> ati_to_amber_type(dict_pointers["NTYPES"]+1);

    std::string line="";
    bool inFlag = false;
    std::vector<int> ati(dict_pointers["NATOM"]);
    std::vector<std::string> atype(dict_pointers["NATOM"]);

    int _j=0;
    file.clear();
    file.seekg(position["AMBER_ATOM_TYPE"]);
    while (getline(file, line)) {
        
        if (line.find("%FLAG AMBER_ATOM_TYPE") != std::string::npos) {
            inFlag = true;
            break;
            //continue;
        }
    }

    if (inFlag) {
        while (getline(file, line)) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Saltamos la línea del formato
                continue;
            }
            std::stringstream ss(line);
            
            
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
        
        if (line.find("%FLAG ATOM_TYPE_INDEX") != std::string::npos) {
            inFlag = true;
            break;
            //continue;
        }
    }

    if (inFlag) {
        while (getline(file, line)) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Saltamos la línea del formato
                continue;
            }
            std::stringstream ss(line);
            
            
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

inline std::vector<int> read_ati(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position)
{
    std::string line="";
    bool inFlag = false;
    std::vector<int> ati(dict_pointers["NATOM"]);
    int _j=0;
    file.clear();
    file.seekg(position);
    while (getline(file, line)) {
        
        if (line.find("%FLAG ATOM_TYPE_INDEX") != std::string::npos) {
            inFlag = true;
            break;
            //continue;
        }
    }

    if (inFlag) {
        while (getline(file, line)) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Saltamos la línea del formato
                continue;
            }
            std::stringstream ss(line);
            
            
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
    //DESCOMENTAR ESTA
    /*
    for (size_t i = 0; i < dict_pointers["NATOM"]; ++i) {
        int _number = ati[i];
        std::cout << "Atom " << i << ": " << "atom_type_index" <<":"<< _number << std::endl; 
    }*/

    return(ati);
}

inline std::vector<float> read_mass(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position)
{
    std::string line="";
    bool inFlag = false;
    std::vector<float> mass(dict_pointers["NATOM"]); 
    int _j=0;
    file.clear();
    file.seekg(position);
    while (getline(file, line)) {
        
        if (line.find("%FLAG MASS") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Saltamos la línea del formato
                continue;
            }
            std::stringstream ss(line);
            
            
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
    //DESCOMENTAR ESTA
    /*
    for (size_t i = 0; i < dict_pointers["NATOM"]; ++i) {
        float _number = mass[i];
        cout << "Atom " << i << ": " << "atom mass" <<":"<< _number << endl; 
    }
    */
    
    return(mass);
}

inline std::map<std::pair<std::string,std::string>,std::pair<float,float>> read_lj(std::ifstream &file, std::map<std::string,int>& dict_pointers,int position,std::vector<std::string> amber_type)
{
    std::string line="";
    bool inFlag = false;
    std::map<std::tuple<int,int>,std::tuple<float,float>> lj_coefficient;
    file.clear();
    file.seekg(position);
    int _j=0;
    int _i=1;//atomo 1
    int _k=1;//atomo 2
    
    while (getline(file, line)) {
        
        if (line.find("%FLAG LENNARD_JONES_ACOEF") != std::string::npos) {
            inFlag = true;
            
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                
                continue;
            }
            
            // Read and split variable values
            std::stringstream ss(line);
               
            for (size_t i = 0; i + 15< line.length(); i += 16) {
                float coef_a=stof(ToolKit::strip(line.substr(i, 16)));
                
                //lj_coefficient[{_i,_k}]={coef_a,0.0};
                std::get<0>(lj_coefficient[{_i,_k}])=coef_a;
                std::get<0>(lj_coefficient[{_k,_i}])=coef_a;
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
    _j=0;
    _i=1;
    _k=1;
    while (getline(file, line)) {
        if (line.find("%FLAG LENNARD_JONES_BCOEF") != std::string::npos) {
            inFlag = true;
            continue;
        }

        if (inFlag) {
            if (line.find("%FORMAT") != std::string::npos) {
                // Skipping format line
                continue;
            }
            
            // Read and split variable values
            std::stringstream ss(line);
               
            for (size_t i = 0; i + 15< line.length(); i += 16) {
                float coef_b=stof(ToolKit::strip(line.substr(i, 16)));

                std::get<1>(lj_coefficient[{_i,_k}])=coef_b;
                std::get<1>(lj_coefficient[{_k,_i}])=coef_b;

                //cout << get<1>(lj_coefficient[{_i,_k}]) << endl;

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
    /*
    for(int _i=1;_i<=4;++_i)
        {
        for(int _k=1;_k<=_i;++_k)
            {
                
                cout << "pair: " << _i << ": " << _k << "--" << get<0>(lj_coefficient[{_i,_k}])  <<":"<< get<1>(lj_coefficient[{_i,_k}]) << endl; 
            }
        }*/
    float sigma=0;
    float epsilon=0;
    float a=0;
    float b=0;
    for(int _i=1;_i<=dict_pointers["NTYPES"];++_i)
        {
        for(int _k=1;_k<=dict_pointers["NTYPES"];++_k)
            {
                a=std::get<0>(lj_coefficient[{_i,_k}]);
                b=std::get<1>(lj_coefficient[{_i,_k}]);
                sigma=pow(a / b, 1.0 / 6.0);
                epsilon=pow(b,2)/(4*a);
                // get<0>(lj_coefficient[{_i,_k}])=sigma; 
                // get<1>(lj_coefficient[{_i,_k}])=epsilon;
                std::get<0>(lj_coefficient[{_i,_k}])=epsilon; 
                std::get<1>(lj_coefficient[{_i,_k}])=sigma;
            }
        }
    std::map<std::pair<std::string,std::string>,std::pair<float,float>> lj_coefficient_in_amber_type;
    for(int _i=1;_i<=dict_pointers["NTYPES"];++_i)
        {
        for(int _k=1;_k<=dict_pointers["NTYPES"];++_k)
            {
                std::pair<std::string,std::string> lj_type=std::make_pair(amber_type[_i],amber_type[_k]);
                std::pair<float,float> lj_coef=std::make_pair(std::get<0>(lj_coefficient[{_i,_k}]),std::get<1>(lj_coefficient[{_i,_k}]));
                lj_coefficient_in_amber_type[lj_type]=lj_coef;
            }
        }
    return(lj_coefficient_in_amber_type);
}

inline std::map<std::string,std::pair<float,float>> read_lj_diagonal(std::map<std::pair<std::string,std::string>,std::pair<float,float>> lj_map)
{
    std::map<std::string,std::pair<float,float>> lj_coefficient;
    for (const auto& par : lj_map) 
    {
        if(std::get<0>(par.first)==std::get<1>(par.first))
        {
            lj_coefficient[std::get<0>(par.first)]=par.second;
        }
    }
    return(lj_coefficient);
}

inline std::map<std::string,int> type_atomic_number(std::map<std::string,int>& dict_pointers,std::vector<int> atom_type_index, std::vector<int> atomic_number, std::vector<std::string> ati_to_amber_type)
{
    std::map<std::string,int> type_atomic_number_index;
    for (int i=0;i<dict_pointers["NATOM"];++i)
    {
        type_atomic_number_index[ati_to_amber_type[atom_type_index[i]]]=atomic_number[i];
    }
    return(type_atomic_number_index);
}
    
inline std::map<std::string, int> flag_position(std::ifstream &file)
{
    std::string line="";
    std::map<std::string, int> flag;
    bool inFlag=false;
    while (getline(file, line)) {
        
        if (line.find("%FLAG") != std::string::npos) {
            line=ToolKit::strip(line.substr(5, 80));
            inFlag = true;
            flag[line]=file.tellg();
            flag[line]-=91;
            //cout<<line<<":"<<flag[line]<<endl;
            continue;
        }
    }
    return(flag);
}


#endif // AMBER_PARSER_HPP