#ifndef CONFIGURATIONBULK_HPP
#define CONFIGURATIONBULK_HPP

/**
 * Version: July 2025
 * Author: Nicolás Loubet
 */

#include "Configuration.hpp"

/**
 * This class creates a Configuration object, with an array of molecules (specific for water bulk systems)
 */
class ConfigurationBulk : public Configuration {
	private:
		// Helper, try to identify a defect
		bool isDefect(const int ID_CENTER, const Real threshold=-12.0, const int V_index=4) {
			return (vI(ID_CENTER, V_index) > threshold);
		}

	public:
		static constexpr int CLASSIFICATION_D_MOLECULE= 0;
		static constexpr int CLASSIFICATION_T0_MOLECULE= 1;
		static constexpr int CLASSIFICATION_T1_MOLECULE= 2;
		static constexpr int CLASSIFICATION_T2_MOLECULE= 3;
		static constexpr int CLASSIFICATION_D3_MOLECULE= 0;
		static constexpr int CLASSIFICATION_D5_MOLECULE= 1;
		static constexpr int CLASSIFICATION_TA_MOLECULE= 2;
		static constexpr int CLASSIFICATION_TB_MOLECULE= 3;

		ConfigurationBulk(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info) :
			Configuration(coord_reader, filename, topol_info) {}

		/**
		 * It indicates if the molecule is a D_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the vI return value, default is -12.0
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4 (V4)
		 * @return If the v4 value is higher that the threshold, or false if it is not Water
		 */
		bool isD(const int ID_CENTER, const Real threshold=-12.0, const int V_index=4) {
			//If it has been classified, it returns that value
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false; //If it is not a Water
			if(m->getClassification() != NOT_CLASSIFIED)
				return (m->getClassification() == CLASSIFICATION_D_MOLECULE || m->getClassification() == CLASSIFICATION_D3_MOLECULE);

			if(isDefect(ID_CENTER, threshold, V_index)) {
				m->setClassification(CLASSIFICATION_D_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * It indicates if the molecule is a D3_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the v4, default is -12.0
		 * @return If the v4 value is higher that the threshold, or false if it is not Water
		 */
		bool isD3(const int ID_CENTER, const Real threshold=-12.0) {
			//If it has been classified, it returns that value
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false; //If it is not a Water
			if(m->getClassification() != NOT_CLASSIFIED)
				return (m->getClassification() == CLASSIFICATION_D_MOLECULE || m->getClassification() == CLASSIFICATION_D3_MOLECULE);

			if(isDefect(ID_CENTER, threshold, 4)) {
				m->setClassification(CLASSIFICATION_D3_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * It indicates if the molecule is a D5_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the v5, default is -12.0
		 * @return If the v5 value is lower that the threshold, or false if it is not Water
		 */
		bool isD5(const int ID_CENTER, const Real threshold=-12.0) {
			//If it has been classified, it returns that value
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false; //If it is not a Water
			if(m->getClassification() != NOT_CLASSIFIED)
				return (m->getClassification() == CLASSIFICATION_D5_MOLECULE);

			if(!isDefect(ID_CENTER, threshold, 5)) {
				m->setClassification(CLASSIFICATION_D5_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * It indicates if the molecule is a D3_MOLECULE or a D5_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the v4 and v5, default is -12.0
		 * @return If the v4 value is higher that the threshold or v5 value is lower, or false if it is not Water
		 */
		bool isDX(const int ID_CENTER, const Real threshold=-12.0) {
			return isD3(ID_CENTER, threshold) || isD5(ID_CENTER, threshold);
		}

		/**
		 * It creates a matrix for the getNeighboursByPotential() function
		 * @return A Real** to the matrix
		 */
		Real** createPotentialMatrix() {
			Real** potential_matrix= new Real*[N_MOLEC];
			for(int i= 1; i < N_MOLEC; i++) {
				potential_matrix[i]= new Real[i];
				for(int j= 0; j < i; j++)
					potential_matrix[i][j]= NOT_CLASSIFIED;
			}
			return potential_matrix;
		}

		/**
		 * It dealocates the matrix generated for the getNeighboursByPotential() function (created with createPotentialMatrix() function)
		 * @param A Real** to the matrix
		 */
		void deletePotentialMatrix(Real** potential_matrix) {
			for(int i= 1; i < N_MOLEC; i++)
				delete[] potential_matrix[i];
			delete[] potential_matrix;
		}

		/**
		 * Define the classification of each molecule in the configuration (D-T2), JCP 2023
		 * The neighbours are searched with the potential
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4
		 * @param threshold The value of potential to which compare the vI return value, default is -12
		 */
		void classifyMolecules(const int V_index=4, const Real threshold=-12.0) {
			const int NUMBER_OF_NEIGHBOURS= 4;

			Real** potential_matrix= createPotentialMatrix();

			//Firstly, I want to know every D molecule
			for(int i= 0; i < N_MOLEC; i++) {
				vector<Real> pots_neighs_vector;
				vector<int> neighs_vector;

				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w == nullptr) continue;

				getNeighboursByPotential(w, pots_neighs_vector, neighs_vector, potential_matrix);
				Sorter::sort(pots_neighs_vector, Sorter::Order::Ascending);

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
					if(w == nullptr) continue;
					if(w->getClassification() != CLASSIFICATION_T2_MOLECULE) continue;

					vector<Real> pots_neighs_vector;
					vector<int> neighs_vector;

					getNeighboursByPotential(w, pots_neighs_vector, neighs_vector, potential_matrix);
					Sorter::cosort(pots_neighs_vector, neighs_vector, Sorter::Order::Ascending);

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
		 * Define the classification of each molecule in the configuration (D3,D5,TA,TB), PRE 2024
		 * The neighbours are searched with the potential
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4
		 * @param threshold The value of potential to which compare the vI return value, default is -12
		 */
		void classifyMolecules_includePentacoordinated(const int V_index=4, const Real threshold=-12.0) {
			const int NUMBER_OF_NEIGHBOURS= 4;
			
			Real** potential_matrix= createPotentialMatrix();

			//Firstly, I want to know every D molecule
			for(int i= 0; i < N_MOLEC; i++) {
				vector<Real> pots_neighs_vector;
				vector<int> neighs_vector;

				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w == nullptr) continue;

				getNeighboursByPotential(w, pots_neighs_vector, neighs_vector, potential_matrix);
				Sorter::cosort(pots_neighs_vector, neighs_vector, Sorter::Order::Ascending);

				if(pots_neighs_vector[V_index-1] > threshold)
					w->setClassification(CLASSIFICATION_D3_MOLECULE);
				else if(pots_neighs_vector[V_index] < threshold)
					w->setClassification(CLASSIFICATION_D5_MOLECULE);
				else
					w->setClassification(CLASSIFICATION_TB_MOLECULE); //At this point, I want that if it's not DX, its assigned TB
			}

			//I have all the DX molecules, and the rest are classified as TB
			//Now, I search for the TB molecules if they have a DX molecule as a 4th potential or less -> TA
			for(int i= 0; i < N_MOLEC; i++) {
				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w == nullptr) continue;
				if(w->getClassification() != CLASSIFICATION_TB_MOLECULE) continue;

				vector<Real> pots_neighs_vector;
				vector<int> neighs_vector;

				getNeighboursByPotential(w, pots_neighs_vector, neighs_vector, potential_matrix);
				Sorter::cosort(pots_neighs_vector, neighs_vector, Sorter::Order::Ascending);

				for(int i_v= 0; i_v < V_index; i_v++) {
					Water* neigh= dynamic_cast<Water*>(molecs[neighs_vector[i_v]-1]);
					
					if(neigh->getClassification() == CLASSIFICATION_D3_MOLECULE ||
					   neigh->getClassification() == CLASSIFICATION_D5_MOLECULE) {

						w->setClassification(CLASSIFICATION_TA_MOLECULE);
						break;
					}
				}
			}
			deletePotentialMatrix(potential_matrix);
		}

		/**
		 * It returns an array with the V_4S values of all molecules
		 * @param inic_value The initial value of the array (first water molecule), default is 1
		 * @param i_V The potential number V_index of a sorted list of all potentials, default is 4
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @param neighbours An array with the neighbours of each molecule, default is nullptr
		 * @return An array with the V_4S values
		 */
		Real* v_4S_arr(const int inic_value= 1, const int i_V= 4, const Real R_CUT_OFF= 5., ToolKit::ArrInt* neighbours= nullptr) {
			Real* output= new Real[N_MOLEC];
			Real** pm= createPotentialMatrix();
			for(int i= inic_value; i < N_MOLEC; i++)
				output[i-inic_value]= getInteractionsPerSite(i,R_CUT_OFF,pm,neighbours)[i_V-1];
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
		Real Tanaka(Water* m, const Real MAX_D_HB= 3.5, const Real MAX_A_HB= 30.) {
			const Real MAX_D_ANALYSIS= 6.0; //Maximum distance for analysis
			Real dist;
			vector<Real> ls_d_HB, ls_d_nHB;

			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue;
				dist= m->distanceTo(getMolec(i+1), bounds);

				if(dist > MAX_D_ANALYSIS) continue;

				Water* w2= dynamic_cast<Water*>(molecs[i]);
				if(w2==nullptr) continue;

				if(m->isHB(*w2, bounds, MAX_D_HB, MAX_A_HB))
					ls_d_HB.push_back(dist);
				else
					ls_d_nHB.push_back(dist);
			}

			Real min_value= *min_element(ls_d_HB.begin(),ls_d_HB.end());
			Real max_value= *max_element(ls_d_nHB.begin(),ls_d_nHB.end());
			Real output= min_value-max_value;
			
			return output;
		}

		/**
		 * It calculates the Local Structure Index value of the molecule specified (https://doi.org/10.1039/C1CP22076D)
		 * @param id The ID of the molecule (1-N)
		 * @return the LSI value
		 */
		Real LSI(const int id) {
			const Real R_MAX= 3.7; //A (CUT-OFF)
			vector<Real> distances;

			Real dist_peripheral= bounds.x+bounds.y+bounds.z;

			for(int j= 1; j <= N_MOLEC; j++) {
				if(j==id) continue;
				Real d= getMolec(id).distanceTo(getMolec(j), bounds);
				if(d > R_MAX) {
					if(d < dist_peripheral)
						dist_peripheral= d;
					continue;
				}
				distances.push_back(d);
			}

			//Add the first molecule outside the sphere R_MAX
			distances.push_back(dist_peripheral);

			Sorter::sort(distances, Sorter::Order::Ascending);
			int N= distances.size()-1;

			Real sum_deltas= 0.;
			Real sum_squared_deltas= 0.;
			for(int j= 0; j < N; j++) {
				Real delta_j= distances[j+1]-distances[j];
				sum_deltas+= delta_j;
				sum_squared_deltas+= delta_j*delta_j;
			}

			Real mean_delta= sum_deltas/N;
			return sum_squared_deltas/N - mean_delta*mean_delta;
		}

};

#endif
