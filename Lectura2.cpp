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

using namespace std;

string strip(string _s)
{
    _s.erase(remove(_s.begin(), _s.end(), ' '), _s.end());
    return(_s);
}


map<string,int> read_pointers(ifstream& file_prmtop,int position)
{
    ifstream file_pointer_names("amber_pointers.txt");
    
    string line;
    vector<string> list_of_names;

    while (getline(file_pointer_names, line)) {
        // Read and split variable names
        stringstream ss(line);
        string this_name;
        while (ss >> this_name) {
            list_of_names.push_back(this_name);
        }
    }
    file_pointer_names.close();

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
    /*for (size_t i = 0; i < list_of_names.size(); ++i) {
        cout << "Variable " << list_of_names[i] << ": " << dict_pointers[list_of_names[i]] << endl; 
    }*/
    return dict_pointers;
}



vector<string> read_atom_name(ifstream &file, map<string,int>& dict_pointers,int position)
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
                atoms[_j]=strip(line.substr(i, 4));
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


vector<tuple<string,float>> read_charge(ifstream &file, map<string,int>& dict_pointers, vector<string> atom_names,int position)
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

vector<int> read_atoms_per_molecule(ifstream &file, map<string,int>& dict_pointers,int position)
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

    int atom_index=0;
    vector<int> atom_number_in_molec(dict_pointers["NATOM"]);
    
    for(int molecule_index=0; molecule_index<nspm;++molecule_index)
    {
        for(int _k=0; _k<atoms_per_molecule[molecule_index]; ++_k)
        {
            atom_number_in_molec[atom_index]= molecule_index;
            ++atom_index;
        }
    }
    return(atom_number_in_molec);
}

vector<int> read_ati(ifstream &file, map<string,int>& dict_pointers,int position)
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
            //continue;
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
                ati[_j]=stoi(strip(line.substr(i, 8)));
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
        cout << "Atom " << i << ": " << "atom_type_index" <<":"<< _number << endl; 
    }*/
    return(ati);
}

vector<float> read_mass(ifstream &file, map<string,int>& dict_pointers,int position)
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
                mass[_j]=stof(strip(line.substr(i, 16)));
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

map<tuple<int,int>,tuple<float,float>> read_lj(ifstream &file, map<string,int>& dict_pointers,int position)
{
    string line="";
    bool inFlag = false;
    map<tuple<int,int>,tuple<float,float>> lj_coefficient;
    file.clear();
    file.seekg(position);
    int _j=0;
    int _i=1;//atomo 1
    int _k=1;//atomo 2
    
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
                float coef_a=stof(strip(line.substr(i, 16)));
                
                //lj_coefficient[{_i,_k}]={coef_a,0.0};
                get<0>(lj_coefficient[{_i,_k}])=coef_a;
                get<0>(lj_coefficient[{_k,_i}])=coef_a;
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
                float coef_b=stof(strip(line.substr(i, 16)));

                get<1>(lj_coefficient[{_i,_k}])=coef_b;
                get<1>(lj_coefficient[{_k,_i}])=coef_b;

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
                a=get<0>(lj_coefficient[{_i,_k}]);
                b=get<1>(lj_coefficient[{_i,_k}]);
                sigma=pow(a / b, 1.0 / 6.0);
                epsilon=pow(b,2)/(4*a);
                get<0>(lj_coefficient[{_i,_k}])=sigma; 
                get<1>(lj_coefficient[{_i,_k}])=epsilon;
            }
        }
    return(lj_coefficient);
}

map<string, int> flag_position(ifstream &file)
{
    string line="";
    map<string, int> flag;
    bool inFlag=false;
    while (getline(file, line)) {
        
        if (line.find("%FLAG") != string::npos) {
            line=strip(line.substr(5, 80));
            inFlag = true;
            flag[line]=file.tellg();
            flag[line]-=81;
            //cout<<line<<":"<<flag[line]<<endl;
            continue;
        }
    }
    return(flag);
}

int main() {
    clock_t start = clock();

    for(int i= 0; i < 1000; i++) {
        ifstream file("jja.prmtop");
        map<string, int> positions_of_flags=flag_position(file);
        map<string,int> dict_pointers = read_pointers(file,positions_of_flags["POINTERS"]);
        vector<tuple<string,float>> map_atoms = read_charge(file, dict_pointers, read_atom_name(file, dict_pointers,positions_of_flags["ATOM_NAME"]),positions_of_flags["CHARGE"]);
        vector<int> atoms_per_molecule = read_atoms_per_molecule(file, dict_pointers,positions_of_flags["SOLVENT_POINTERS"]); 
        vector<int> atom_type_index=read_ati(file,dict_pointers,positions_of_flags["ATOM_TYPE_INDEX"]); 
        vector<float> mass=read_mass(file,dict_pointers,positions_of_flags["MASS"]);
        map<tuple<int,int>,tuple<float,float>> lj_coefficient=read_lj(file,dict_pointers,positions_of_flags["LENNARD_JONES_ACOEF"]);
    
    }

    clock_t end = clock();
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    cout << "Tiempo de ejecucion: " << elapsed_secs << endl;
    return 0;
}
