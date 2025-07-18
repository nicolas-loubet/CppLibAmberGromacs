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
        /**
         * This function find the nearest atoms in the same molecule
         * @param ID_MOLEC   ID of the molecule
         * @param ID_CENTER  ID of the atom center to search
         * @param D_MAX_NEI  float that indicate maximum bond distanec. 1.65 is OK.
         * @return A struct array of int s
         */
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

        /**
         * This function find the CH3 groups in molecule. It do not add CH3 of Choline head group.
         * @param ID_MOLEC   ID of the molecule
         * @return a vector of IDs generally 2
         */

        inline vector<int> findCH3(const int ID_MOLEC) const
            {
            vector<int> list_CH3;
            int n_CH3=0;
            for(int i= 1; i <= getMolec(ID_MOLEC).getNAtoms(); i++) //for each atom in molecule
                {
                ToolKit::ArrInt NearbyAtoms;
                if(getMolec(ID_MOLEC).getAtom(i).getZ()==6) //if it is C
                    {
                    NearbyAtoms=findNearbyAtoms(ID_MOLEC, i, 1.65);//1.58 is maximum distance found in simulation
                    bool exit=false;
                    if(NearbyAtoms.size>0)
                        {
                        for(int j=0; j<NearbyAtoms.size; j++)
                            {
                            if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==7) // Discards CH3 bound to N
                                {
                                exit=true;  
                                }
                            }
                        if(exit==true){ exit=false; continue;}

                        int n_hydrogen=0;
                        for(int j=0; j<NearbyAtoms.size; j++)
                            {
                            if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==1) {n_hydrogen+=1;}  // if hydrogen is found
                            }
                        if(n_hydrogen==3)
                            {
                            list_CH3.insert(list_CH3.begin(),i); // add CH3 carbon atom ID to the list
                            n_CH3+=1; 
                            }

                        }
                    }
                }
                return list_CH3;
            }

        /**
         * This function analize the chain from terminal CH3 going up until it finds an oxigen atom. 
         * @param ID_MOLEC   ID of the molecule
         * @param ID_CH3  ID of the CH3 residue to search up the molecule
         * @return returns a pair in which the first term holds if the chain is SN1 or SN1 and
         * a vector of maps in which the first term hold the carbon number and the second term holds the nearby atoms, generally 2 C and 2
         */

        inline pair<int,vector<map<int,ToolKit::ArrInt>>> analizeChain(const int ID_MOLEC, const int ID_CH3) const
            {
            ToolKit::ArrInt NearbyAtoms;
            vector<map<int,ToolKit::ArrInt>> current_chain;
            NearbyAtoms=findNearbyAtoms(ID_MOLEC, ID_CH3, 1.65);
            int previous_c=ID_CH3;
            int next_C=0;
            int chain_lenght=0;
            bool exit=false;
            int SN_hydrogen=0;
            
            map<int,ToolKit::ArrInt> primer_atomo_analizado;
            primer_atomo_analizado[previous_c]=NearbyAtoms;   
            current_chain.insert(current_chain.begin(),primer_atomo_analizado);// Adds CH3 to the chain

            for(int j=0; j<NearbyAtoms.size; j++) //searchs one time for the next carbon
                {
                if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==6) 
                    {
                    if(NearbyAtoms.arr[j]!=previous_c)
                        {
                        next_C=NearbyAtoms.arr[j]; // and sets as the next
                        chain_lenght+=1;
                        break;
                        }
                    }  
                }

            int failsafe=0; //filesafe to avoid infinite while()
            while(exit==false)
                {
                NearbyAtoms=findNearbyAtoms(ID_MOLEC, next_C, 1.65); //find next atom neighbours
                failsafe+=1;
                if(chain_lenght>=22 || failsafe>25){exit=true; break;}

                for(int j=0; j<NearbyAtoms.size; j++)
                    {
                    if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==8) // An Oxigen has been found
                        {
                        ToolKit::ArrInt NearbyAtoms_oxigen=findNearbyAtoms(ID_MOLEC, NearbyAtoms.arr[j], 1.65);
                        for(int k=0; k<NearbyAtoms_oxigen.size; k++)
                            {
                            if(getMolec(ID_MOLEC).getAtom(NearbyAtoms_oxigen.arr[k]).getZ()==6) //Search in the oxigen for the next carbon
                                {
                                if(NearbyAtoms_oxigen.arr[k]!=previous_c) // that is not the previous, so it is SN1 or SN2 Carbon
                                    {
                                    ToolKit::ArrInt NearbyAtoms_SN2=findNearbyAtoms(ID_MOLEC, NearbyAtoms_oxigen.arr[k], 1.65);    
                                    for(int l=0; l<NearbyAtoms_SN2.size; l++)
                                        {
                                        if(getMolec(ID_MOLEC).getAtom(NearbyAtoms_SN2.arr[l]).getZ()==1) {SN_hydrogen+=1;}  //Count the number of Hydrogen to determne if it is SN1 or SN2
                                        }
                                    exit=true;
                                    break;
                                    }
                                }
                            }
                        if(exit==true) // if oxigen found, add C to the list and break from the while loop
                            {
                            map<int,ToolKit::ArrInt> analyzed_atom;
                            analyzed_atom[next_C]=NearbyAtoms;
                            current_chain.insert(current_chain.begin(),analyzed_atom);
                            chain_lenght+=1;
                            break;
                            }
                        continue;
                        }
                    if(getMolec(ID_MOLEC).getAtom(NearbyAtoms.arr[j]).getZ()==6) // if no Oxigen was found and Carbon was found
                        {
                        if(NearbyAtoms.arr[j]!=previous_c) // And it is no the previuos Carbon, thien is the next.
                            {
                            map<int,ToolKit::ArrInt> analyzed_atom;
                            analyzed_atom[next_C]=NearbyAtoms;
                            current_chain.insert(current_chain.begin(),analyzed_atom);//en vez de insert tiene que ser un pushback para considerar cadenas distintas y quelos ch3 coincidan
                            previous_c=next_C; 
                            next_C=NearbyAtoms.arr[j];//moves to the next
                            chain_lenght+=1;//Adds the carbon to the chain
                            break;
                            }
                        }  
                    }
                }
            return make_pair(SN_hydrogen,current_chain);
            }
        
        /**
         * This function analize order parameter of the molecule. 
         * @param ID_MOLEC   ID of the molecule
         * @return returns a list of atoms and the order parameter of the chain.
         */

        inline vector<float> Molecule_analizer(const int ID_MOLEC) const
            {
            vector<float> order_per_carbon(22);  // sets a maximun chain lenght of 22 
            vector<int> CH3_found = findCH3(ID_MOLEC); // finds the terminal CH3
            for(int i=0; i<CH3_found.size();i++) //for each CH3 found
                {
                pair<int,vector<map<int,ToolKit::ArrInt>>> cadena = analizeChain(ID_MOLEC, CH3_found[i]); //Chain is analyzed
                vector<map<int,ToolKit::ArrInt>> cadena_analizada=cadena.second; //Grabs the chain
                
                for(int j=0; j<cadena_analizada.size(); j++)
                    {
                    map<int,ToolKit::ArrInt> analyzed_atom=cadena_analizada[j];
                    float order=0.0;
                    int n_hydrogen=0;
                    for(const auto& mol_pair : analyzed_atom)
                        {
                        int carbono_numero=mol_pair.first; //obteins the carbon
                        ToolKit::ArrInt lista=mol_pair.second;
                        for(int k=0;k<lista.size;k++)
                            {
                            if(getMolec(ID_MOLEC).getAtom(lista.arr[k]).getZ()==1) //finds the Hydrogen
                                {
                                Vector CH=getMolec(ID_MOLEC).getAtom(carbono_numero).getPosition() - getMolec(ID_MOLEC).getAtom(lista.arr[k]).getPosition(); //Gets the unit vector of the CH bond
                                float z = (CH/CH.magnitude()).z; //calculates the cosene (proyection of the unit vector with z)
                                float cos2_tita=z*z;
                                order+=0.5*(3*cos2_tita-1);
                                n_hydrogen+=1;
                                }
                            }

                        }
                        if(n_hydrogen==0){order=0;n_hydrogen=1;}
                        order_per_carbon[j]+=order/(n_hydrogen*CH3_found.size()); //adds order from the hydrogen considered divided by the total number of Hydrogen atoms for that Carbon
                    }
                }
            return order_per_carbon;
            }
	public:
        ConfigurationLipid(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info) :
        Configuration(coord_reader, filename, topol_info) {}
        
        vector<float> orderParameter(const int ID_MOLEC) {
            return Molecule_analizer(ID_MOLEC);
        }
		
};

#endif
