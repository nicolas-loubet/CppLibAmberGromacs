#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

/**
 * Version: July 2025
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

		/**
		 * Adds the potential of an atom with a water molecule to the sum_per_site vector, used in getInteractionsPerSite
		 */
		void addToSumVector(vector<Vector>& sites, vector<float>& sum_per_site, Water& center_water, Atom& atom, const float R_CUT_OFF) {
			int i_close= 0;
			float d_close= distancePBC(sites[0],atom.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				float d_new= distancePBC(sites[i],atom.getPosition(),bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			if(d_close <= R_CUT_OFF) {
				sum_per_site[i_close]+= center_water.potentialWith(atom,bounds);
			}
		}

		/**
		 * Adds the potential of a water molecule with another water molecule to the sum_per_site vector, used in getInteractionsPerSite
		 */
		void addToSumVector(vector<Vector>& sites, vector<float>& sum_per_site, Water& center_water, Water& other, float** potential_matrix, ToolKit::ArrInt* neighbours, const float R_CUT_OFF) {
			int i_close= 0;
			float d_close= distancePBC(sites[0],other.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				float d_new= distancePBC(sites[i],other.getPosition(),bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			if(d_close > R_CUT_OFF) return;

			float pot;
			if(potential_matrix != nullptr) {
				int min= (center_water.getID() < other.getID()  ?  center_water.getID() : other.getID()) - 1;
				int max= (center_water.getID() > other.getID()  ?  center_water.getID() : other.getID()) - 1;
				if(potential_matrix[max][min] == NOT_CLASSIFIED)
					potential_matrix[max][min]= center_water.potentialWith(other,bounds);
				pot= potential_matrix[max][min];
				if(neighbours != nullptr)
					neighbours[center_water.getID()-1].arr[
						neighbours[center_water.getID()-1].size++
					]= other.getID();
			} else {
				pot= center_water.potentialWith(other,bounds);
			}
			sum_per_site[i_close]+= pot;
		}

		/**
		 * Adds the potential of a water molecule with another water molecule to the sum_per_site vector, used in getInteractionsPerSite flagged with water-water interactions < V_CUT_OFF
		 */
		void addToSumVector(vector<Vector>& sites, vector<float>& sum_per_site, Water& center_water, Water& other, float** potential_matrix, ToolKit::ArrInt* neighbours, const float R_CUT_OFF, bool* ww_interaction, const float V_CUT_OFF) {
			int i_close= 0;
			float d_close= distancePBC(sites[0],other.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				float d_new= distancePBC(sites[i],other.getPosition(),bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			if(d_close > R_CUT_OFF) return;

			float pot;
			if(potential_matrix != nullptr) {
				int min= (center_water.getID() < other.getID()  ?  center_water.getID() : other.getID()) - 1;
				int max= (center_water.getID() > other.getID()  ?  center_water.getID() : other.getID()) - 1;
				if(potential_matrix[max][min] == NOT_CLASSIFIED)
					potential_matrix[max][min]= center_water.potentialWith(other,bounds);
				pot= potential_matrix[max][min];
				if(neighbours != nullptr)
					neighbours[center_water.getID()-1].arr[neighbours[center_water.getID()-1].size++]= other.getID();
			} else {
				pot= center_water.potentialWith(other,bounds);
			}
			sum_per_site[i_close]+= pot;
			ww_interaction[i_close]= ww_interaction[i_close] || (pot <= V_CUT_OFF);
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

		//Getters
		int getNMolec() const { return N_MOLEC; }
		const Molecule& getMolec(int id) const { return (*molecs[id-1]); }
		Molecule* getMolec_ptr(int id) const { return molecs[id-1]; }
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
		 * Copy constructor
		 */
		Configuration(const Configuration& other) : N_MOLEC(other.N_MOLEC), bounds(other.bounds), molecs(nullptr) {
			molecs= new Molecule*[N_MOLEC];
			for(int i = 0; i < N_MOLEC; i++)
				molecs[i]= other.molecs[i] ? new Molecule(*other.molecs[i]) : nullptr;
		}

		/**
		 * Default constructor
		 */
		Configuration() : N_MOLEC(0), bounds(Vector(0,0,0)), molecs(nullptr) {}

		Configuration& operator=(const Configuration& other) {
			if(this == &other) return *this;
			for(int i= 0; i < N_MOLEC; i++)
				delete molecs[i];
			delete[] molecs;

			N_MOLEC= other.N_MOLEC;
			bounds= other.bounds;
			molecs= new Molecule*[N_MOLEC];
			for(int i= 0; i < N_MOLEC; i++)
				molecs[i]= other.molecs[i] ? new Molecule(*other.molecs[i]) : nullptr;
			return *this;
		}
	
		/**
		 * Destructor. It destroys the molecule array
		 */
		~Configuration() {
			for(int i= 0; i < N_MOLEC; i++)
				delete molecs[i];
			delete[] molecs;
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
				if(getMolec(ID_CENTER).distanceTo(getMolec(i), bounds) <= D_MAX_NEI)
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
				if(getMolec(ID_CENTER).distanceTo(getMolec(j+1), bounds) > MAX_V4) continue; //Cutoff use
				Water* w= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
				Water* w2= dynamic_cast<Water*>(molecs[j]);
				if(w==nullptr || w2==nullptr) continue;
				ls_V[ls_V_i++]= w->potentialWith(*w2, bounds);
			}

			Sorter::sort(ls_V, ls_V_i, Sorter::Order::Ascending);
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
		 * @param threshold The value of potential to which compare the vI return value, default is -12.0
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4 (V4)
		 * @return If the v4 value is higher that the threshold, or false if it is not Water
		 */
		bool isD(const int ID_CENTER, const float threshold=-12.0, const int V_index=4) {
			//If it has been classified, it returns that value
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false; //If it is not a Water
			if(m->getClassification() != NOT_CLASSIFIED) return m->getClassification() == CLASSIFICATION_D_MOLECULE;

			float v4= vI(ID_CENTER, V_index);

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
		bool isD3(const int ID_CENTER, const float threshold=-12.0) {
			return isD(ID_CENTER, threshold);
		}

		/**
		 * It indicates if the molecule is a D5_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the v5, default is -12.0
		 * @return If the v5 value is lower that the threshold, or false if it is not Water
		 */
		bool isD5(const int ID_CENTER, const float threshold=-12.0) {
			return !isD(ID_CENTER, threshold, 5);
		}

		/**
		 * It indicates if the molecule is a D3_MOLECULE or a D5_MOLECULE (and also assigns this value to the classification so you don't have to use this function two times)
		 * @param ID_CENTER int The ID of the Molecule to check
		 * @param threshold The value of potential to which compare the v4 and v5, default is -12.0
		 * @return If the v4 value is higher that the threshold or v5 value is lower, or false if it is not Water
		 */
		bool isDX(const int ID_CENTER, const float threshold=-12.0) {
			return isD3(ID_CENTER, threshold) || isD5(ID_CENTER, threshold);
		}

		/**
		 * Finds the water molecules nearby the specified one in terms of potential energy
		 * @param m *Water that is the center of the search
		 * @param pots *float To return the potential
		 * @param identificators *ToolKit::ArrInt (use &) To return the ids of the molecules for each potential
		 * @param potential_matrix **float where to register the potentials to avoid two times search
		 */
		void getNeighboursByPotential(Water* m, vector<float>& pots, vector<int>& identificators, float** potential_matrix) {
			const float MAX_R_V4= 5.5; //Cutoff to not compare all the molecules

			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue;
				if(m->distanceTo(getMolec(i+1), bounds) > MAX_R_V4) continue;
				Water* w2= dynamic_cast<Water*>(molecs[i]);
				if(w2==nullptr) continue;

				int min= i < m->getID()-1 ? i:m->getID()-1;
				int max= i > m->getID()-1 ? i:m->getID()-1;

				if(potential_matrix[max][min] == NOT_CLASSIFIED)
					potential_matrix[max][min]= m->potentialWith(*w2, bounds);
				pots.push_back(potential_matrix[max][min]);
				identificators.push_back(i+1);
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
				delete[] potential_matrix[i];
			delete[] potential_matrix;
		}

		/**
		 * Define the classification of each molecule in the configuration (D-T2), JCP 2023
		 * The neighbours are searched with the potential
		 * @param V_index Same of vI [See vI(m,V_index)], default is 4
		 * @param threshold The value of potential to which compare the vI return value, default is -12
		 */
		void classifyMolecules(const int V_index=4, const float threshold=-12.0) {
			const int NUMBER_OF_NEIGHBOURS= 4;

			float** potential_matrix= createPotentialMatrix();

			//Firstly, I want to know every D molecule
			for(int i= 0; i < N_MOLEC; i++) {
				vector<float> pots_neighs_vector;
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

					vector<float> pots_neighs_vector;
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
		void classifyMolecules_includePentacoordinated(const int V_index=4, const float threshold=-12.0) {
			const int NUMBER_OF_NEIGHBOURS= 4;
			
			float** potential_matrix= createPotentialMatrix();

			//Firstly, I want to know every D molecule
			for(int i= 0; i < N_MOLEC; i++) {
				vector<float> pots_neighs_vector;
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

				vector<float> pots_neighs_vector;
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
		 * It returns the interactions per site of the molecule specified by its ID
		 * @param ID The ID of the molecule
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The interactions per site, sorted in descending order
		 */
		vector<float> getInteractionsPerSite(const int ID, const float R_CUT_OFF= 5., float** potential_matrix= nullptr, ToolKit::ArrInt* neighbours= nullptr) {
			Water& molecule= *dynamic_cast<Water*>(molecs[ID-1]);
			Vector o= molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			Geometrics::TetrahedronVertices t= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds);
			
			vector<float> sum_per_site(4,0.0f);
			vector<Vector> sites= {t.H1, t.H2, t.L1, t.L2};
			
			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue; // Same molecule, skip
				
				if((getMolec(j+1).getNAtoms() == 1) && (getMolec(j+1).getCharge() >= 1 || getMolec(j+1).getCharge() <= -1)) {
					// We found a ion, consider it always
					addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(1), R_CUT_OFF);
					continue;
				}
				
				Water* other= dynamic_cast<Water*>(molecs[j]);
				if(other != nullptr) { // We found a water molecule
					if(molecule.distanceTo(*other,bounds) > R_CUT_OFF+1.1) continue;
					addToSumVector(sites, sum_per_site, molecule, *other, potential_matrix, neighbours, R_CUT_OFF);
				} else { // We found another type of molecule, study atom by atom
					for(int a= 1; a <= getMolec(j+1).getNAtoms(); a++)
						if(molecule.distanceTo(getMolec(j+1).getAtom(a),bounds) < R_CUT_OFF+1.1) {
							addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(a), R_CUT_OFF);
						}
				}
			}

			Sorter::sort(sum_per_site, Sorter::Order::Ascending);
			return sum_per_site;
		}

		/**
		 * It returns the interactions per site of the molecule specified by its ID, with a flag to know if there is a water-water interaction in each site
		 * @param ID The ID of the molecule
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @param V_CUT_OFF The potential cutoff, default is -12
		 * @param flag_ww A flag to know if there is a water-water interaction
		 * @return The interactions per site, sorted in descending order
		 */
		vector<float> getInteractionsPerSite(const int ID, bool& flag_ww, const float V_CUT_OFF= -12.0f, const float R_CUT_OFF= 5., float** potential_matrix= nullptr, ToolKit::ArrInt* neighbours= nullptr) {
			Water& molecule= *dynamic_cast<Water*>(molecs[ID-1]);
			Vector o= molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			Geometrics::TetrahedronVertices t= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds);
			
			vector<float> sum_per_site(4,0.0f);
			vector<Vector> sites= {t.H1, t.H2, t.L1, t.L2};

			bool* water_water_interaction= new bool[4];
			for(int i= 0; i < 4; i++)
				water_water_interaction[i]= false;
			
			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue; // Same molecule, skip
				
				if((getMolec(j+1).getNAtoms() == 1) && (getMolec(j+1).getCharge() >= 1 || getMolec(j+1).getCharge() <= -1)) {
					// We found a ion, consider it always
					addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(1), R_CUT_OFF);
					continue;
				}
				
				Water* other= dynamic_cast<Water*>(molecs[j]);
				if(other != nullptr) { // We found a water molecule
					if(molecule.distanceTo(*other,bounds) > R_CUT_OFF+1.1) continue;
					addToSumVector(sites, sum_per_site, molecule, *other, potential_matrix, neighbours, R_CUT_OFF, water_water_interaction, V_CUT_OFF);
				} else { // We found another type of molecule, study atom by atom
					for(int a= 1; a <= getMolec(j+1).getNAtoms(); a++)
						if(molecule.distanceTo(getMolec(j+1).getAtom(a),bounds) < R_CUT_OFF+1.1) {
							addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(a), R_CUT_OFF);
						}
				}
			}

			flag_ww= true;
			for(int iv= 0; iv < 4; iv++) {
				if(!water_water_interaction[iv]) {
					flag_ww= false;
					break;
				}
			}

			Sorter::sort(sum_per_site, Sorter::Order::Descending);
			return sum_per_site;
		}

		/**
		 * It returns the V_4S value of the molecule specified by its ID
		 * @param ID The ID of the molecule
		 * @param i_V The potential number V_index of a sorted list of all potentials, default is 4
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The V_4S value
		 */
		float v_4S(const int ID, const int i_V= 4, const float R_CUT_OFF= 5.) {
			return getInteractionsPerSite(ID, R_CUT_OFF, nullptr, nullptr)[i_V-1];
		}

		/**
		 * It returns an array with the V_4S values of all molecules
		 * @param inic_value The initial value of the array (first water molecule), default is 1
		 * @param i_V The potential number V_index of a sorted list of all potentials, default is 4
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @param neighbours An array with the neighbours of each molecule, default is nullptr
		 * @return An array with the V_4S values
		 */
		float* v_4S_arr(const int inic_value= 1, const int i_V= 4, const float R_CUT_OFF= 5., ToolKit::ArrInt* neighbours= nullptr) {
			float* output= new float[N_MOLEC];
			float** pm= createPotentialMatrix();
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
		float Tanaka(Water* m, const float MAX_D_HB= 3.5, const float MAX_A_HB= 30.) {
			const float MAX_D_ANALYSIS= 6.0; //Maximum distance for analysis
			float dist;
			vector<float> ls_d_HB, ls_d_nHB;

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
				float d= getMolec(id).distanceTo(getMolec(j), bounds);
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
