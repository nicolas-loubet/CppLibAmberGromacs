#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

/**
 * Version: May 2025
 * Author: Nicolás Loubet
 */

#include "../Molecules/Water.hpp"
#include "../General/ToolKit.hpp"
#include "../General/Sorter.hpp"
#include "../General/Geometrics.hpp"
#include "../ReaderInterfaces.hpp"

using namespace std;

/**
 * This class creates a Configuration object, with an array of molecules (generic)
 */
class Configuration {
	protected:
		Molecule** molecs; //Array of molecule pointers
		int N_MOLEC; //The number of Molecule objects in the array
		Vector bounds; //The bounds of the system

	public:
		static constexpr int CLASSIFICATION_D_MOLECULE= 0;
		static constexpr int CLASSIFICATION_T0_MOLECULE= 1;
		static constexpr int CLASSIFICATION_T1_MOLECULE= 2;
		static constexpr int CLASSIFICATION_T2_MOLECULE= 3;

		//Getters
		int getNMolec() const { return N_MOLEC; }
		Molecule* getMolec(int id) const { return molecs[id-1]; }
		Vector getBounds() const { return bounds; }

		/**
		 * Constructor. It reads the topology and coordinates files
		 * @param coord_reader Object that reads the coordinates file, given as pointer
		 * @param topol_info Object that reads the topology file
		 */
		Configuration(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info) {
			N_MOLEC= topol_info.num_molecules;
			molecs= new Molecule*[N_MOLEC];
	
			if(!coord_reader->readCoordinates(filename, topol_info, molecs, bounds))
				throw runtime_error("Failed to read coordinates");
		}
	
		/**
		 * Destructor. It destroys the molecule array
		 */
		~Configuration() {
			for(int i= 0; i < N_MOLEC; i++)
				delete(molecs[i]);
			delete(molecs);
		}

		/**
		 * Finds the molecules nearby the specified one
		 * @param ID_CENTER int The ID of the Molecule that is the center of the search
		 * @param D_MAX_NEI float Maximum radium of neighbour search
		 * @return A ToolKit::ArrInt with the IDs of the molecules that are at a distance of D_MAX_NEI Angstrom or less
		 */
		ToolKit::ArrInt findNearby(const int ID_CENTER, const float D_MAX_NEI) const {
			int* i_nearby= new int[N_MOLEC];
			int counter= 0;

			for(int i= 1; i <= N_MOLEC; i++) {
				if(i == ID_CENTER) continue; //If it is the same molecule
				if(molecs[ID_CENTER-1]->distanceTo(*molecs[i-1], bounds) <= D_MAX_NEI)
					i_nearby[counter++]= i;
			}

			ToolKit::ArrInt output;
			output.arr= i_nearby;
			output.size= counter;
			return output;
		}

		/**
		 * It returns the list of all the potentials detected ordered by magnitud
		 * @param ID_CENTER int The ID of the Molecule that is the center of the search
		 * @param MAX_V4 The radii of cut-off to not compare all the molecules
		 * @return The potential number V_index of a sorted list of all potentials
		 */
		float* getVList(const int ID_CENTER, const float MAX_V4= 5.5) {
			int ls_V_i= 0;
			float* ls_V= new float[N_MOLEC];

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID_CENTER) continue; //If they are the same molecule
				if(molecs[ID_CENTER-1]->distanceTo(*molecs[j], bounds) > MAX_V4) continue; //Cutoff use
				Water* w= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
				Water* w2= dynamic_cast<Water*>(molecs[j]);
				if(w==nullptr || w2==nullptr) continue;
				ls_V[ls_V_i++]= w->potentialWith(*w2, bounds);
			}

			Sorter::sort(ls_V, ls_V_i, true);
			return ls_V;
		}

		/**
		 * It calculates the ith major potential of the Water molecule with another
		 * @param ID_CENTER int The ID of the Molecule to calculate
		 * @param V_index The i to return the ith potential, default is 4 (V4)
		 * @return The potential number V_index of a sorted list of all potentials
		 */
		float vI(const int ID_CENTER, const int V_index=4) {
			float *ls_V= getVList(ID_CENTER);
			float v_i= ls_V[V_index-1];
			delete(ls_V);
			return v_i;
		}

		/**
		 * It indicates if the molecule is a D_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4 (V4)
		 * @param threshold The value of potential to which compare the vI return value, default is -12.0
		 * @return If the vI value if higher that the threshold, or false if it is not Water
		 */
		bool isD(const int ID_CENTER, const int V_index=4, const float threshold=-12.0) {
			//If it has been classified, it returns that value
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m==nullptr) return false; //If it is not a Water
			if(m->getClassification() != NOT_CLASSIFIED) return m->getClassification()==CLASSIFICATION_D_MOLECULE;

			float v4= vI(ID_CENTER, V_index);

			if(v4 > threshold) {
				m->setClassification(CLASSIFICATION_D_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * It expands recursively the classification of molecules if the system is already classified in D_MOLECULES and T2_MOLECULES
		 * @param ID_CENTER int The ID of the Water molecule to check
		 */
		void expandClassification(const int ID_CENTER) {
			ToolKit::ArrInt neighbours= findNearby(ID_CENTER,3.2);

			Water* molec= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(molec==nullptr) return;
			for(int i= 0; i < neighbours.size; i++) {
				Water* molec2= dynamic_cast<Water*>(molecs[neighbours.arr[i]-1]);
				if(molec2==nullptr) continue;
				if(molec2->getClassification() <= molec->getClassification()+1) continue;

				molec2->setClassification(molec->getClassification()+1);
				expandClassification(neighbours.arr[i]); //When you don't replace it, it's not called anymore; this is indirectly a cut parameter
			}
			delete(neighbours.arr);
		}

		/**
		 * Finds the water molecules nearby the specified one in terms of potential energy
		 * @param m *Water that is the center of the search
		 * @param pots *float To return the potential
		 * @param identificators *ToolKit::ArrInt (use &) To return the ids of the molecules for each potential
		 * @param potential_matrix **float where to register the potentials to avoid two times search
		 */
		void getNeighboursByPotential(Water* m, vector<float>& pots, vector<int>& identificators, float** potential_matrix) {
			const float MAX_V4= 5.5; //Cutoff to not compare all the molecules

			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue; //If they are the same molecule
				if(m->distanceTo(*molecs[i], bounds) > MAX_V4) continue; //Cutoff use
				Water* w2= dynamic_cast<Water*>(molecs[i]);
				if(w2==nullptr) continue;
				int min= i<m->getID()-1?i:m->getID()-1;
				int max= i>m->getID()-1?i:m->getID()-1;
				if(potential_matrix[max][min] == NOT_CLASSIFIED) {
					float pot= m->potentialWith(*w2, bounds);
					pots.push_back(pot);
					potential_matrix[max][min]= pot;
					identificators.push_back(i+1);
				} else {
					pots.push_back(potential_matrix[max][min]);
					identificators.push_back(i+1);
				}
			}
		}

		/**
		 * It creates a matrix for the getNeighboursByPotential() function
		 * @return A float** to the matrix
		 */
		float** createPotentialMatrix() {
			float** potential_matrix= new float*[N_MOLEC];
			for(int i= 1; i < N_MOLEC; i++) {
				potential_matrix[i]= new float[i];
				for(int j= 0; j < i; j++)
					potential_matrix[i][j]= NOT_CLASSIFIED;
			}
			return potential_matrix;
		}

		/**
		 * It dealocates the matrix generated for the getNeighboursByPotential() function (created with createPotentialMatrix() function)
		 * @param A float** to the matrix
		 */
		void deletePotentialMatrix(float** potential_matrix) {
			for(int i= 1; i < N_MOLEC; i++)
				delete(potential_matrix[i]);
			delete(potential_matrix);
		}

		/**
		 * Define the classification of each molecule in the configuration
		 * The neighbours are searched with the potential
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4
		 * @param threshold The value of potential to which compare the vI return value, default is -12
		 */
		void classifyMolecules(const int V_index=4, const float threshold=-12.0) {
			const int NUMBER_OF_NEIGHBOURS= 4;

			float** potential_matrix= createPotentialMatrix();

			//Firstly, I want to know every D molecule
			ToolKit::ArrInt ids;
			for(int i= 0; i < N_MOLEC; i++) {
				vector<float> pots_neighs_vector;
				vector<int> neighs_vector;

				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w==nullptr) continue;

				getNeighboursByPotential(w, pots_neighs_vector, neighs_vector, potential_matrix);
				Sorter::sort(pots_neighs_vector,true);

				if(pots_neighs_vector[V_index-1] > threshold)
					w->setClassification(CLASSIFICATION_D_MOLECULE);
				else
					w->setClassification(CLASSIFICATION_T2_MOLECULE); //At this point, I want that if it's not D, its assigned T2
			}

			//I have all the D molecules, and the rest are classified as T2
			//Now, I search for the T2 molecules if they have a D molecule as a 4th potential or less -> T0
			//The same for T1, searching T0 neighbours
			for(int t_order= CLASSIFICATION_T0_MOLECULE; t_order <= CLASSIFICATION_T1_MOLECULE; t_order++)
				for(int i= 0; i < N_MOLEC; i++) {
					Water* w= dynamic_cast<Water*>(molecs[i]);
					if(w==nullptr) continue;
					if(w->getClassification() != CLASSIFICATION_T2_MOLECULE) continue;

					vector<float> pots_neighs_vector;
					vector<int> neighs_vector;

					getNeighboursByPotential(w, pots_neighs_vector, neighs_vector, potential_matrix);
					Sorter::cosort(pots_neighs_vector, neighs_vector, true);

					for(int i_v= 0; i_v < V_index; i_v++) {
						Water* neigh= dynamic_cast<Water*>(molecs[neighs_vector[i_v]-1]);
						
						if(neigh->getClassification() == t_order-1) {
							w->setClassification(t_order);
							break;
						}
					}
				}

			deletePotentialMatrix(potential_matrix);
		}

		/**
		 * It returns the interactions per site of the molecule specified by its ID
		 * @param ID The ID of the molecule
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The interactions per site, sorted in descending order
		 */
		vector<float> getInteractionsPerSite(const int ID, float** potential_matrix= nullptr, ToolKit::ArrInt* neighbours= nullptr, const float R_CUT_OFF= 5.) {
			Water* molecule= dynamic_cast<Water*>(molecs[ID-1]);
			if(molecule == nullptr) throw runtime_error("Error: getInteractionsPerSite(ID, float**, ToolKit::ArrInt*, const float R_CUT_OFF) -> molecule is not a water");
			Vector o= molecule->getOxygen().getPosition();
			Vector h1= molecule->getHydrogen_1().getPosition();
			Vector h2= molecule->getHydrogen_2().getPosition();
			Geometrics::TetrahedronVertices t= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds);

			vector<float> sum_per_site(4,0.0f);
			Vector sites[4]= {t.H1, t.H2, t.L1, t.L2};
			
			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1==ID) continue;
				bool consider_ion= false;
				if(molecule->distanceTo(*molecs[j], bounds) > R_CUT_OFF+1.1) {
					if(-1 < molecs[j]->getCharge() && molecs[j]->getCharge() < 1) continue;
					consider_ion= true;
				}
				int i_close= 0;
				float d_close= distancePBC(sites[0],molecs[j]->getPosition(),bounds);
				for(int i= 1; i < 4; i++) {
					float d_new= distancePBC(sites[i],molecs[j]->getPosition(),bounds);
					if(d_close > d_new) {
						i_close= i;
						d_close= d_new;
					}
				}
				if(d_close <= R_CUT_OFF || consider_ion) {
					if(potential_matrix != nullptr) {
						int min= ID-1 < j  ?  ID-1 : j;
						int max= ID-1 > j  ?  ID-1 : j;
						if(potential_matrix[max][min] == NOT_CLASSIFIED)
							potential_matrix[max][min]= molecule->potentialWith(*molecs[j],bounds);
						sum_per_site[i_close]+= potential_matrix[max][min];
						if(neighbours != nullptr)
							neighbours[ID-1].arr[neighbours[ID-1].size++]= j+1;
					} else
						sum_per_site[i_close]+= molecule->potentialWith(*molecs[j],bounds);
				}
			}

			Sorter::sort(sum_per_site,false);
			return sum_per_site;
		}

		/**
		 * It returns the V_4S value of the molecule specified by its ID
		 * @param ID The ID of the molecule
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The V_4S value
		 */
		float v_4S(const int ID, const float R_CUT_OFF= 5.) {
			vector<float> interactions= getInteractionsPerSite(ID, nullptr, nullptr, R_CUT_OFF);
			float v4s= interactions[3];
			return v4s;
		}

		/**
		 * It returns an array with the V_4S values of all molecules
		 * @param i_V The potential number V_index of a sorted list of all potentials, default is 4
		 * @param neighbours An array with the neighbours of each molecule, default is nullptr
		 * @return An array with the V_4S values
		 */
		float* v_4S_arr(const int i_V= 4, ToolKit::ArrInt* neighbours= nullptr) {
			float* output= new float[N_MOLEC];
			float** pm= createPotentialMatrix();
			for(int i= 0; i < N_MOLEC; i++)
				output[i]= getInteractionsPerSite(i+1,pm,neighbours)[i_V-1];
			deletePotentialMatrix(pm);
			return output;
		}

		/**
		 * It calculates the Tanaka ζ value of the molecule specified
		 * @param m The *Water to analyse
		 * @param MAX_D_HB The maximum distance O-O that it could be considereded an HB (near 3.5A)
		 * @param MAX_A_HB The angle O-O-H that it could be considereded an HB (near 30°)
		 * @return The potential number V_index of a sorted list of all potentials
		 */
		float getTanaka(Water* m, const float MAX_D_HB= 3.5, const float MAX_A_HB= 30.) {
			const float MAX_D_ANALYSIS= 6.0; //Maximum distance for analysis
			float dist;
			vector<float> ls_d_HB, ls_d_nHB;

			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue;
				dist= m->distanceTo(*molecs[i], bounds);

				if(dist > MAX_D_ANALYSIS) continue;

				Water* w2= dynamic_cast<Water*>(molecs[i]);
				if(w2==nullptr) continue;

				if(m->isHB(*w2, bounds, MAX_D_HB, MAX_A_HB))
					ls_d_HB.push_back(dist);
				else
					ls_d_nHB.push_back(dist);
			}

			float min_value= *min_element(ls_d_HB.begin(),ls_d_HB.end());
			float max_value= *max_element(ls_d_nHB.begin(),ls_d_nHB.end());
			float output= min_value-max_value;
			
			return output;
		}

		/**
		 * It calculates the Local Structure Index value of the molecule specified (https://doi.org/10.1039/C1CP22076D)
		 * @param id The ID of the molecule (1-N)
		 * @return the LSI value
		 */
		float LSI(const int id) {
			const float R_MAX= 3.7; //A (CUT-OFF)
			vector<float> distances;

			float dist_peripheral= bounds.x+bounds.y+bounds.z;

			for(int j= 1; j <= N_MOLEC; j++) {
				if(j==id) continue;
				float d= molecs[id-1]->distanceTo(*molecs[j-1], bounds);
				if(d > R_MAX) {
					if(d < dist_peripheral)
						dist_peripheral= d;
					continue;
				}
				distances.push_back(d);
			}

			//Add the first molecule outside the sphere R_MAX
			distances.push_back(dist_peripheral);

			Sorter::sort(distances, true);
			int N= distances.size()-1;

			float sum_deltas= 0.;
			float sum_squared_deltas= 0.;
			for(int j= 0; j < N; j++) {
				float delta_j= distances[j+1]-distances[j];
				sum_deltas+= delta_j;
				sum_squared_deltas+= delta_j*delta_j;
			}

			float mean_delta= sum_deltas/N;
			return sum_squared_deltas/N - mean_delta*mean_delta;
		}

};

#endif
