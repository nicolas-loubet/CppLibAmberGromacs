#ifndef CONFIGURATIONLIPID_HPP
#define CONFIGURATIONLIPID_HPP

/**
 * Version: July 2025
 * Author: Ezequiel Cuenca
 */

#include "Configuration.hpp"
#include <array>
#include <vector>

/**
 * This class creates a Configuration object, with an array of molecules (specific for lipid systems)
 */
class ConfigurationLipid : public Configuration {
    private:

        inline ToolKit::ArrInt findNearbyAtoms(const int ID_MOLEC, const int ID_CENTER, const float D_MAX_NEI) const {
			int* i_nearby= new int[getMolec(ID_MOLEC).getNAtoms()];
			int counter= 0;

			for(int i= 1; i <= getMolec(ID_MOLEC).getNAtoms(); i++) {
				if(i == ID_CENTER) continue; //If it is the same Atom
				if(getMolec(ID_MOLEC).getAtom(ID_CENTER).distanceTo(getMolec(ID_MOLEC).getAtom(i), bounds) <= D_MAX_NEI)
					i_nearby[counter++]= i;
			}

			ToolKit::ArrInt output;
			output.arr= i_nearby;
			output.size= counter;
			return output;
		}
        
        inline vector<int> findCH3(const int ID_MOLEC) const
            {
            vector<int> list_CH3;
            int n_CH3=0;
            for(int i= 1; i <= getMolec(ID_MOLEC).getNAtoms(); i++) 
                {
                ToolKit::ArrInt NearbyAtoms;
                if(getMolec(ID_MOLEC).getAtom(i).getZ()==6)
                    {
                    NearbyAtoms=findNearbyAtoms(ID_MOLEC, i, 1.65);//1.58
                    if(NearbyAtoms.size==3)
                        {
                        if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[0]).getZ()==7 ||
                            getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[1]).getZ()==7 ||
                            getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[2]).getZ()==7 ||
                            getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[3]).getZ()==7)
                            {
                                continue;  
                            }
                        int n_hidrogen=0;
                        for(int j=0; j<3; j++)
                            {
                            if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==1) {n_hidrogen+=1;}  
                            }
                        if(n_hidrogen==3)
                            {
                                list_CH3[n_CH3]=i; 
                                n_CH3+=1; 
                            }

                        }
                    }
                }
                return list_CH3;
            }

        inline pair<int,vector<map<int,ToolKit::ArrInt>>> analizeChain(const int ID_MOLEC, const int ID_CH3) const
            {
            ToolKit::ArrInt NearbyAtoms;
            vector<map<int,ToolKit::ArrInt>> cadena;
            NearbyAtoms=findNearbyAtoms(ID_MOLEC, ID_CH3, 1.65);
            int anterior=ID_CH3;
            int siguiente=0;
            int largo_cadena=0;
            bool salir=false;
            int SN_hidrogen=0;

            for(int j=0; j<NearbyAtoms.size; j++)
                {
                if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==6) 
                    {
                    if(NearbyAtoms.arr[j]!=anterior)
                        {
                        map<int,ToolKit::ArrInt> atomo_analizado;
                        atomo_analizado[anterior]=NearbyAtoms;
                        siguiente=NearbyAtoms.arr[j];
                        cadena.insert(cadena.begin(),atomo_analizado);
                        //cadena[largo_cadena]=atomo_analizado;
                        largo_cadena+=1;
                        break;
                        }
                    }  
                }
            int failsafe=0;
            while(salir==false)
                {
                NearbyAtoms=findNearbyAtoms(ID_MOLEC, siguiente, 1.65);
                failsafe+=1;
                if(largo_cadena>=20 || failsafe>25){salir=true; break;}

                for(int j=0; j<NearbyAtoms.size; j++)
                    {
                    if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==8) 
                        {
                        ToolKit::ArrInt NearbyAtoms_oxigen=findNearbyAtoms(ID_MOLEC, NearbyAtoms.arr[j], 1.65);
                        for(int k=0; k<NearbyAtoms_oxigen.size; k++)
                            {
                            if(getMolec(ID_MOLEC).getAtom(NearbyAtoms_oxigen.arr[k]).getZ()==6) 
                                {
                                if(NearbyAtoms_oxigen.arr[k]!=anterior)
                                    {
                                    ToolKit::ArrInt NearbyAtoms_SN2=findNearbyAtoms(ID_MOLEC, NearbyAtoms_oxigen.arr[k], 1.65);    
                                    for(int l=0; l<NearbyAtoms_SN2.size; l++)
                                        {
                                        if(getMolec(ID_MOLEC).getAtom(NearbyAtoms_SN2.arr[l]).getZ()==1) {SN_hidrogen+=1;}  
                                        }
                                    salir=true;
                                    break;
                                    }
                                }
                            }
                        if(salir==true){break;}
                        continue;
                        }
                    if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==6) 
                        {
                        if(NearbyAtoms.arr[j]!=anterior)
                            {
                            map<int,ToolKit::ArrInt> atomo_analizado;
                            atomo_analizado[anterior]=NearbyAtoms;
                            anterior=siguiente;
                            siguiente=NearbyAtoms.arr[j];
                            cadena.insert(cadena.begin(),atomo_analizado);
                            //cadena[largo_cadena]=atomo_analizado;
                            largo_cadena+=1;
                            break;
                            }
                        }  
                    }
                }
            return make_pair(SN_hidrogen,cadena);
            }

        inline int Molecule_analizer(const int ID_MOLEC) const
            {
            vector<int> CH3_found = findCH3(ID_MOLEC);
            pair<int,vector<map<int,ToolKit::ArrInt>>> cadena = analizeChain(ID_MOLEC, CH3_found[0]);
            vector<map<int,ToolKit::ArrInt>> cadena_analizada=cadena.second;
            for(int j=0; j<cadena_analizada.size(); j++)
                {
                map<int,ToolKit::ArrInt> atomo_analizado=cadena_analizada[0];
                for(const auto& mol_pair : atomo_analizado)
                    {
                    int carbono_numero=mol_pair.first;
                    cout<< "C"<<carbono_numero<<"[";  
                    ToolKit::ArrInt lista=mol_pair.second;
                    for(int k=0;k<lista.size;k++)
                        {
                        cout<<lista.arr[k]<<",";    
                        }
                    cout<<"]"<<endl;
                    }
                }
             return 1;
            }

	public:
        ConfigurationLipid(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info) :
        Configuration(coord_reader, filename, topol_info) {}
        
        array<vector<float>,2> orderParameter() {
            array<vector<float>,2> op;
            //Ejemplo de insertar al comienzo un valor (10.5), hace el "kick" automático
            op[0].insert(op[0].begin(), 10.5f);
            //TODO
            //throw runtime_error("Method not implemented");
            int test = Molecule_analizer(1);
            return op;
        }
		
};

#endif
