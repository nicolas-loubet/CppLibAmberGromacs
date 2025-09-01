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

		// Helper: finds the index of the closest site
		inline bool closestSiteIndex(const vector<Vector>& sites, const Vector& pos, const Vector& bounds, Real R_CUT_OFF, int& i_close) {
			i_close= 0;
			Real d_close= distancePBC(sites[0], pos, bounds);
			for(int i= 1; i < 4; i++) {
				Real d_new= distancePBC(sites[i], pos, bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			return (d_close <= R_CUT_OFF);
		}

		// Helper: check in the potential matrix the stored value
		inline Real checkInPotentialMatrix(Water& center_water, Water& other, Real** potential_matrix, ToolKit::ArrInt* neighbours) {
			int min= (center_water.getID() < other.getID() ? center_water.getID() : other.getID()) - 1;
			int max= (center_water.getID() > other.getID() ? center_water.getID() : other.getID()) - 1;

			if(potential_matrix[max][min] == NOT_CLASSIFIED)
				potential_matrix[max][min]= center_water.potentialWith(other, bounds);

			Real pot= potential_matrix[max][min];
			if(neighbours != nullptr)
				neighbours[center_water.getID()-1].arr[neighbours[center_water.getID()-1].size++]= other.getID();
			return pot;
		}

		// Helper: Adds the potential of an atom with a water molecule to the sum_per_site vector, used in getInteractionsPerSite
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Atom& atom, const Real R_CUT_OFF) {
			int i_close;
			if(closestSiteIndex(sites, atom.getPosition(), bounds, R_CUT_OFF, i_close))
				sum_per_site[i_close]+= center_water.potentialWith(atom, bounds);
		}

		// Helper: Adds the potential of a water molecule with another water molecule to the sum_per_site vector, used in getInteractionsPerSite
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Water& other,
							Real** potential_matrix, ToolKit::ArrInt* neighbours, const Real R_CUT_OFF) {
			int i_close;
			if(!closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close))
				return;

			sum_per_site[i_close]+= (potential_matrix != nullptr) ?
									checkInPotentialMatrix(center_water, other, potential_matrix, neighbours) :
									center_water.potentialWith(other, bounds);
		}


		// Helper: Adds the potential of a water molecule with another water molecule to the sum_per_site vector, used in getInteractionsPerSite
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water,
							Water& other, const Real R_CUT_OFF, vector<Real>& sum_only_water) {
			int i_close;
			if(!closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close))
				return;

			Real pot= center_water.potentialWith(other, bounds);
			sum_per_site[i_close]+= pot;
			sum_only_water[i_close]+= pot;
		}


		// Helper: Adds the potential of a water molecule with another water molecule to the sum_per_site vector,
		// used in getInteractionsPerSite flagged with water-water interactions < V_CUT_OFF
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Water& other, Real** potential_matrix,
							ToolKit::ArrInt* neighbours, const Real R_CUT_OFF, bool* ww_interaction, const Real V_CUT_OFF) {
			int i_close;
			if(!closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close))
				return;

			Real pot= (potential_matrix != nullptr) ?
					  checkInPotentialMatrix(center_water, other, potential_matrix, neighbours) :
					  center_water.potentialWith(other, bounds);

			sum_per_site[i_close]+= pot;
			ww_interaction[i_close]= ww_interaction[i_close] || (pot <= V_CUT_OFF);
		}

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
		void initializeSiteVectors(vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices, vector<Real>& sum_per_site) {
			ww_interactions.assign(4, vector<Real>());
			ww_distances.assign(4, vector<Real>());
			ww_indices.assign(4, vector<int>());
			sum_per_site.assign(4, 0.0);
		}

	public:
		//Getters
		int getNMolec() const { return N_MOLEC; }
		const Molecule& getMolec(int id) const { return (*molecs[id-1]); }
		Molecule* getMolec_ptr(int id) const { return molecs[id-1]; }
		Vector getBounds() const { return bounds; }

		/**
		 * Getters for water molecules, with cast
		 * @param id The ID of the molecule
		 * @return The water molecule
		 */
		const Water& getWater(int id) const {
			if(molecs[id-1]->isWater())
				return *((Water*) molecs[id-1]);
			throw runtime_error("The molecule is not a water");
		}

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
			for(int i= 0; i < N_MOLEC; i++)
				molecs[i]= other.molecs[i] ? new Molecule(*other.molecs[i]) : nullptr;
		}

		/**
		 * Default constructor
		 */
		Configuration() : N_MOLEC(0), bounds(Vector(0,0,0)), molecs(nullptr) {}

		Configuration& operator=(const Configuration& other) {
			if(this == &other) return *this;
			delete(this);

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
			molecs= nullptr;
		}

		/**
		 * Finds the molecules nearby the specified one
		 * @param ID_CENTER int The ID of the Molecule that is the center of the search
		 * @param D_MAX_NEI Real Maximum radium of neighbour search
		 * @return A ToolKit::ArrInt with the IDs of the molecules that are at a distance of D_MAX_NEI Angstrom or less
		 */
		ToolKit::ArrInt findNearby(const int ID_CENTER, const Real D_MAX_NEI) const {
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
		Real* getVList(const int ID_CENTER, const Real MAX_V4= 5.5) {
			if(!getMolec(ID_CENTER).isWater()) return nullptr;
			Water w1= *static_cast<Water*>(molecs[ID_CENTER-1]);

			int ls_V_i= 0;
			Real* ls_V= new Real[N_MOLEC];

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID_CENTER) continue; //If they are the same molecule
				if(!getMolec(j+1).isWater()) continue;
				if(getMolec(ID_CENTER).distanceTo(getMolec(j+1), bounds) > MAX_V4) continue; //Cutoff use
				Water w2= *static_cast<Water*>(molecs[j]);
				ls_V[ls_V_i++]= w1.potentialWith(w2, bounds);
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
		Real vI(const int ID_CENTER, const int V_index=4) {
			Real *ls_V= getVList(ID_CENTER);
			Real v_i= ls_V[V_index-1];
			delete(ls_V);
			return v_i;
		}

		/**
		 * Finds the water molecules nearby the specified one in terms of potential energy
		 * @param m *Water that is the center of the search
		 * @param pots vector<Real> To return the potential
		 * @param identificators vector<int> To return the ids of the molecules for each potential
		 * @param potential_matrix **Real where to register the potentials to avoid two times search
		 */
		void getNeighboursByPotential(Water* m, vector<Real>& pots, vector<int>& identificators, Real** potential_matrix) {
			const Real MAX_R_V4= 5.5; //Cutoff to not compare all the molecules

			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue;
				if(!getMolec(i+1).isWater()) continue;
				if(m->distanceTo(getMolec(i+1), bounds) > MAX_R_V4) continue;
				Water* w2= static_cast<Water*>(molecs[i]);

				int min= i < m->getID()-1 ? i:m->getID()-1;
				int max= i > m->getID()-1 ? i:m->getID()-1;

				if(potential_matrix[max][min] == NOT_CLASSIFIED)
					potential_matrix[max][min]= m->potentialWith(*w2, bounds);
				pots.push_back(potential_matrix[max][min]);
				identificators.push_back(i+1);
			}
		}

		/**
		 * It returns the interactions per site of the molecule specified by its ID
		 * @param ID The ID of the molecule
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @param labels Identification of the sorted list, default is nullptr. If you need it, declare it ass "new int[4]", result would be 0-3
		 * @return The interactions per site, sorted in descending order
		 */
		vector<Real> getInteractionsPerSite(const int ID, const Real R_CUT_OFF= 5.0, Real** potential_matrix= nullptr, ToolKit::ArrInt* neighbours= nullptr, int* labels= nullptr) {
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			Vector o= molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			vector<Vector> sites= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds).toVector();
			
			vector<Real> sum_per_site(4,0.0);
			
			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue; // Same molecule, skip
				
				if((getMolec(j+1).getNAtoms() == 1) && (getMolec(j+1).getCharge() >= 1 || getMolec(j+1).getCharge() <= -1)) {
					// We found a ion, consider it always
					addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(1), R_CUT_OFF);
					continue;
				}
				
				if(molecs[j]->isWater()) { // We found a water molecule
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other,bounds) > R_CUT_OFF+1.1) continue;
					addToSumVector(sites, sum_per_site, molecule, other, potential_matrix, neighbours, R_CUT_OFF);
				} else { // We found another type of molecule, study atom by atom
					for(int a= 1; a <= getMolec(j+1).getNAtoms(); a++)
						if(molecule.distanceTo(getMolec(j+1).getAtom(a),bounds) < R_CUT_OFF+1.1) {
							addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(a), R_CUT_OFF);
						}
				}
			}

			if(labels == nullptr) {
				Sorter::sort(sum_per_site, Sorter::Order::Ascending);
			} else {
				vector<int> vector_labels= {0,1,2,3};
				Sorter::cosort(sum_per_site, vector_labels, Sorter::Order::Ascending);
				for(int i= 0; i < 4; i++) labels[i]= vector_labels[i];
			}
			return sum_per_site;
		}

		/**
		 * It returns the interactions per site of the molecule specified by its ID, with a flag to know if there is a water-water interaction in each site
		 * @param ID The ID of the molecule
		 * @param flag_ww A flag to know if there is a water-water interaction
		 * @param V_CUT_OFF The potential cutoff, default is -12
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @return The interactions per site, sorted in descending order
		 */
		vector<Real> getInteractionsPerSite(const int ID, bool& flag_ww, const Real V_CUT_OFF= -12.0, const Real R_CUT_OFF= 5.0, Real** potential_matrix= nullptr, ToolKit::ArrInt* neighbours= nullptr) {
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			Vector o= molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			vector<Vector> sites= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds).toVector();
			
			vector<Real> sum_per_site(4,0.0);

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
				
				if(molecs[j]->isWater()) { // We found a water molecule
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other,bounds) > R_CUT_OFF+1.1) continue;
					addToSumVector(sites, sum_per_site, molecule, other, potential_matrix, neighbours, R_CUT_OFF, water_water_interaction, V_CUT_OFF);
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
		 * It returns the interactions per site of the molecule specified by its ID, but it only considers water-water interactions
		 * @param ID The ID of the molecule
		 * @param V_CUT_OFF The potential cutoff, default is -12
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The interactions per site, sorted in descending order
		 */
		vector<Real> getInteractionsPerSite_waterOnly(const int ID, const Real V_CUT_OFF= -12.0, const Real R_CUT_OFF= 5.0) {
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			Vector o= molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			vector<Vector> sites= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds).toVector();
			
			vector<Real> sum_per_site(4,0.0);
			vector<Real> sum_only_water(4,0.0);

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue; // Same molecule, skip
				
				if((getMolec(j+1).getNAtoms() == 1) && (getMolec(j+1).getCharge() >= 1 || getMolec(j+1).getCharge() <= -1)) {
					// We found a ion, consider it always
					addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(1), R_CUT_OFF);
					continue;
				}
				
				if(molecs[j]->isWater()) { // We found a water molecule
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other,bounds) > R_CUT_OFF+1.1) continue;
					addToSumVector(sites, sum_per_site, molecule, other, R_CUT_OFF, sum_only_water);
				} else { // We found another type of molecule, study atom by atom
					for(int a= 1; a <= getMolec(j+1).getNAtoms(); a++)
						if(molecule.distanceTo(getMolec(j+1).getAtom(a),bounds) < R_CUT_OFF+1.1) {
							addToSumVector(sites, sum_per_site, molecule, getMolec(j+1).getAtom(a), R_CUT_OFF);
						}
				}
			}

			vector<int> index_seq= {0,1,2,3};
			Sorter::cosort(sum_per_site, index_seq, Sorter::Order::Descending);
			return {sum_only_water[index_seq[0]], sum_only_water[index_seq[1]], sum_only_water[index_seq[2]], sum_only_water[index_seq[3]]};
		}

		/**
		 * It returns the V_4S value of the molecule specified by its ID
		 * @param ID The ID of the molecule
		 * @param i_V The potential number V_index of a sorted list of all potentials, default is 4
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The V_4S value
		 */
		Real v_4S(const int ID, const int i_V= 4, const Real R_CUT_OFF= 5.0) {
			return getInteractionsPerSite(ID, R_CUT_OFF, nullptr, nullptr)[i_V-1];
		}

		struct DefectInfo {
			bool is_DJ;
			bool is_D3;
			bool is_D5;
			int lacking_sites;
			int bifurcated_sites;
			vector<Real> sum_per_site;
			pair<Real,Real> bifurcated_individual_potentials;
			pair<Real,Real> bifurcated_individual_distances;
			pair<int,int> bifurcated_individual_indices;
			Real bifurcated_site_potential;
			Real lacking_site_potential;
			Vector lacking_site_position;
			Vector bifurcated_site_position;

			DefectInfo(): is_DJ(false), is_D3(false), is_D5(false), lacking_sites(0), bifurcated_sites(0) {}

			void chargeData(vector<Vector>& sites, vector<Real>& sum_per_site, vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices) {
				this->sum_per_site= sum_per_site;
				bool hasD3= false, hasD5= false;
				for(int i= 0; i < 4; ++i) {
					if(ww_interactions[i].size() == 0) {
						hasD3= true;
						lacking_sites++;
						lacking_site_position= sites[i];
						lacking_site_potential= sum_per_site[i];
					} else if(ww_interactions[i].size() >= 2) {
						hasD5= true;
						bifurcated_sites++;
						bifurcated_site_position= sites[i];
						bifurcated_site_potential= sum_per_site[i];
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
				if(!is_DJ) {
					is_D3= hasD3;
					is_D5= hasD5;
				}
			}
		};
		/**
		 * It indicates if the molecule is DJ according to per-site counts:
		 * DJ if there exists at least one site with 0 neighbours (D3-like) and at least one site with >=2 neighbours (D5-like),
		 * where neighbours are waters whose potential with the center is <= V_CUT_OFF.
		 * @param ID int The ID of the Molecule to check
		 * @param R_CUT_OFF Cutoff radius to consider neighbours (default 5)
		 * @param V_CUT_OFF Potential threshold to count neighbours (default -12)
		 * @return DefectInfo: Information about the molecule, see struct DefectInfo
		 */
		DefectInfo classifyDefect(const int ID, const Real R_CUT_OFF= 5.0, const Real V_CUT_OFF= -12.0) {
			DefectInfo output;
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");

			Water& molecule= *static_cast<Water*>(molecs[ID - 1]);
			Vector o= molecule.getOxygen().getPosition();
			Vector h1= molecule.getHydrogen_1().getPosition();
			Vector h2= molecule.getHydrogen_2().getPosition();
			Geometrics::TetrahedronVertices t= Geometrics::getPerfectTetrahedron(o, h1, h2, bounds);
			vector<Vector> sites= {t.H1, t.H2, t.L1, t.L2};

			vector<vector<Real>> ww_interactions, ww_distances;
			vector<vector<int>> ww_indices;
			vector<Real> sum_per_site;
			initializeSiteVectors(ww_interactions, ww_distances, ww_indices, sum_per_site);

			for(int j= 0; j < N_MOLEC; ++j) {
				if(j+1 == ID) continue;

				const Molecule& other_molec= getMolec(j + 1);
				if(other_molec.getNAtoms() == 1 && (other_molec.getCharge() >= 1 || other_molec.getCharge() <= -1)) {
					processIonInteraction(molecule, other_molec, sites, bounds, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				} else if(molecs[j]->isWater()) {
					Water& other= *static_cast<Water*>(molecs[j]);
					processWaterInteraction(molecule, other, sites, bounds, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				} else {
					processSoluteInteraction(molecule, other_molec, sites, bounds, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				}
			}

			Sorter::sort(sum_per_site, Sorter::Order::Descending);
			output.chargeData(sites, sum_per_site, ww_interactions, ww_distances, ww_indices);
			return output;
		}

};

#endif
