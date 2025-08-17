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
		// Helper function to find the closest site and its distance
		pair<int,Real> findClosestSite(vector<Vector>& sites, const Vector& target, const Vector& bounds) {
			int closest_idx= 0;
			Real closest_dist= distancePBC(sites[0], target, bounds);
			for (int i= 1; i < 4; ++i) {
				Real dist= distancePBC(sites[i], target, bounds);
				if(dist < closest_dist) {
					closest_idx= i;
					closest_dist= dist;
				}
			}
			return {closest_idx, closest_dist};
		}

		// Helper function to process interactions with an ion
		void processIonInteraction(Water& molecule, const Molecule& ion, vector<Vector>& sites, const Vector& bounds, Real R_CUT_OFF, Real V_CUT_OFF,
			                       vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices, vector<Real>& sum_per_site) {
			auto [closest_idx, closest_dist]= findClosestSite(sites, ion.getAtom(1).getPosition(), bounds);
			if(closest_dist > R_CUT_OFF) return;
			Real pot= molecule.potentialWith(ion.getAtom(1), bounds);
			sum_per_site[closest_idx]+= pot;
			if(pot <= V_CUT_OFF) {
				ww_interactions[closest_idx].push_back(pot);
				ww_distances[closest_idx].push_back(molecule.distanceTo(ion, bounds));
				ww_indices[closest_idx].push_back(ion.getID());
			}
		}

		// Helper function to process interactions with another water molecule
		void processWaterInteraction(Water& molecule, Water& other, vector<Vector>& sites, const Vector& bounds, Real R_CUT_OFF, Real V_CUT_OFF,
			                         vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices, vector<Real>& sum_per_site) {
			Real d_ww= molecule.distanceTo(other, bounds);
			if(d_ww > R_CUT_OFF + 1.1) return;

			auto [closest_idx, closest_dist]= findClosestSite(sites, other.getPosition(), bounds);
			if(closest_dist > R_CUT_OFF) return;
			Real pot= molecule.potentialWith(other, bounds);
			sum_per_site[closest_idx]+= pot;
			if(pot <= V_CUT_OFF) {
				ww_interactions[closest_idx].push_back(pot);
				ww_distances[closest_idx].push_back(d_ww);
				ww_indices[closest_idx].push_back(other.getID());
			}
		}

		// Helper function to process interactions with a non-water molecule
		void processSoluteInteraction(Water& molecule, const Molecule& solute, vector<Vector>& sites, const Vector& bounds, Real R_CUT_OFF, Real V_CUT_OFF,
			                          vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices, vector<Real>& sum_per_site) {
			vector<Real> total_potential(4, 0.0);
			vector<Real> min_distance(4, numeric_limits<Real>::infinity());

			for(int a= 1; a <= solute.getNAtoms(); a++) {
				if(molecule.distanceTo(solute.getAtom(a), bounds) >= R_CUT_OFF + 1.1) continue;

				auto [closest_idx, closest_dist]= findClosestSite(sites, solute.getAtom(a).getPosition(), bounds);
				if(closest_dist > R_CUT_OFF) return;
				Real pot= molecule.potentialWith(solute.getAtom(a), bounds);
				sum_per_site[closest_idx]+= pot;
				total_potential[closest_idx]+= pot;
				if(closest_dist < min_distance[closest_idx]) min_distance[closest_idx]= closest_dist;
			}

			for(int i= 0; i < 4; ++i) {
				if(total_potential[i] > V_CUT_OFF) continue;
				ww_interactions[i].push_back(total_potential[i]);
				ww_distances[i].push_back(min_distance[i]);
				ww_indices[i].push_back(solute.getID());
				break;
			}
		}

		// Helper function to initialize site interaction vectors
		void initializeSiteVectors(vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<Real>& sum_per_site) {
			ww_interactions.assign(4, vector<Real>(4, 0.0));
			ww_distances.assign(4, vector<Real>(4, 0.0));
			sum_per_site.assign(4, 0.0);
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
			if(m->getClassification() != NOT_CLASSIFIED) return m->getClassification() == CLASSIFICATION_D_MOLECULE;

			Real v4= vI(ID_CENTER, V_index);

			if(v4 > threshold) {
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
			return isD(ID_CENTER, threshold);
		}

		/**
		 * It indicates if the molecule is a D5_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the v5, default is -12.0
		 * @return If the v5 value is lower that the threshold, or false if it is not Water
		 */
		bool isD5(const int ID_CENTER, const Real threshold=-12.0) {
			return !isD(ID_CENTER, threshold, 5);
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

		struct DJInfo {
			bool is_DJ;
			vector<Real> sum_per_site;
			pair<Real,Real> bifurcated_individual_potentials;
			pair<Real,Real> bifurcated_individual_distances;
			pair<int,int> bifurcated_individual_indices;
			Real bifurcated_site_potential;
			Real lacking_site_potential;
			Vector lacking_site_position;
			Vector bifurcated_site_position;

			DJInfo(): is_DJ(false) {}

			void chargeData(vector<Vector>& sites, vector<Real>& sum_per_site, vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices) {
				sum_per_site= sum_per_site;
				bool hasD3= false, hasD5= false;
				for(int i= 0; i < 4; ++i) {
					if(ww_interactions[i].size() == 0) {
						hasD3 = true;
						lacking_site_position= sites[i];
						lacking_site_potential= sum_per_site[i];
					} else if(ww_interactions[i].size() >= 2) {
						hasD5 = true;
						bifurcated_site_position = sites[i];
						bifurcated_site_potential = sum_per_site[i];
						if(ww_interactions[i][0] < ww_interactions[i][1]) {
							bifurcated_individual_potentials= {ww_interactions[i][0], ww_interactions[i][1]};
							bifurcated_individual_distances= {ww_distances[i][0], ww_distances[i][1]};
							bifurcated_individual_indices= {ww_indices[i][0], ww_indices[i][1]};
						} else {
							bifurcated_individual_potentials= {ww_interactions[i][1], ww_interactions[i][0]};
							bifurcated_individual_distances= {ww_distances[i][1], ww_distances[i][0]};
							bifurcated_individual_indices= {ww_indices[i][1], ww_indices[i][0]};
						}
					}
				}
				is_DJ= hasD3 && hasD5;
			}
		};
		/**
         * It indicates if the molecule is DJ according to per-site counts:
         * DJ if there exists at least one site with 0 neighbours (D3-like)
         * and at least one site with >=2 neighbours (D5-like),
         * where neighbours are waters whose potential with the center is <= V_CUT_OFF.
         *
         * @param ID int The ID of the Molecule to check
         * @param R_CUT_OFF Cutoff radius to consider neighbours (default 5)
         * @param V_CUT_OFF Potential threshold to count neighbours (default -12)
         * @return DJInfo: is_DJ, sum_per_site, bifurcated_individual_potentials, bifurcated_indivisual_distances, bifurcated_site_potential, lacking_site_potential, lacking_site_position, bifurcated_site_position
         */
        DJInfo isDJ(const int ID, const Real R_CUT_OFF = 5.0, const Real V_CUT_OFF = -12.0) {
			DJInfo output;
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");

			Water& molecule= *static_cast<Water*>(molecs[ID - 1]);
			Vector o = molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			Geometrics::TetrahedronVertices t= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds);
			vector<Vector> sites= {t.H1, t.H2, t.L1, t.L2};

			vector<vector<Real>> ww_interactions, ww_distances;
			vector<vector<int>> ww_indices;
			vector<Real> sum_per_site;
			initializeSiteVectors(ww_interactions, ww_distances, sum_per_site);

			for(int j= 0; j < N_MOLEC; ++j) {
				if(j+1 == ID) continue;

				const Molecule& other_molec= getMolec(j + 1);
				if(other_molec.getNAtoms() == 1 && (other_molec.getCharge() >= 1 || other_molec.getCharge() <= -1)) {
					processIonInteraction(molecule, other_molec, sites, bounds, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				} else if(molecs[j]->isWater()) {
					Water& other = *static_cast<Water*>(molecs[j]);
					processWaterInteraction(molecule, other, sites, bounds, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				} else {
					processSoluteInteraction(molecule, other_molec, sites, bounds, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				}
			}

			Sorter::sort(sum_per_site, Sorter::Order::Descending);
			output.chargeData(sites, sum_per_site, ww_interactions, ww_distances, ww_indices);
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
