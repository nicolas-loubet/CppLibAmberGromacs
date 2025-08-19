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
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Atom& atom, const Real R_CUT_OFF) {
			int i_close= 0;
			Real d_close= distancePBC(sites[0],atom.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				Real d_new= distancePBC(sites[i],atom.getPosition(),bounds);
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
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Water& other, Real** potential_matrix, ToolKit::ArrInt* neighbours, const Real R_CUT_OFF) {
			int i_close= 0;
			Real d_close= distancePBC(sites[0],other.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				Real d_new= distancePBC(sites[i],other.getPosition(),bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			if(d_close > R_CUT_OFF) return;

			Real pot;
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
		 * Adds the potential of a water molecule with another water molecule to the sum_per_site vector, used in getInteractionsPerSite
		 */
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Water& other, const Real R_CUT_OFF, vector<Real>& sum_only_water) {
			int i_close= 0;
			Real d_close= distancePBC(sites[0],other.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				Real d_new= distancePBC(sites[i],other.getPosition(),bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			if(d_close > R_CUT_OFF) return;

			Real pot= center_water.potentialWith(other,bounds);
			sum_per_site[i_close]+= pot;
			sum_only_water[i_close]+= pot;
		}

		/**
		 * Adds the potential of a water molecule with another water molecule to the sum_per_site vector, used in getInteractionsPerSite flagged with water-water interactions < V_CUT_OFF
		 */
		void addToSumVector(vector<Vector>& sites, vector<Real>& sum_per_site, Water& center_water, Water& other, Real** potential_matrix, ToolKit::ArrInt* neighbours, const Real R_CUT_OFF, bool* ww_interaction, const Real V_CUT_OFF) {
			int i_close= 0;
			Real d_close= distancePBC(sites[0],other.getPosition(),bounds);
			for(int i= 1; i < 4; i++) {
				Real d_new= distancePBC(sites[i],other.getPosition(),bounds);
				if(d_close > d_new) {
					i_close= i;
					d_close= d_new;
				}
			}
			if(d_close > R_CUT_OFF) return;

			Real pot;
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
			for(int i = 0; i < N_MOLEC; i++)
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
		 * @param pots *Real To return the potential
		 * @param identificators *ToolKit::ArrInt (use &) To return the ids of the molecules for each potential
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
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @return The interactions per site, sorted in descending order
		 */
		vector<Real> getInteractionsPerSite(const int ID, const Real R_CUT_OFF= 5.0, Real** potential_matrix= nullptr, ToolKit::ArrInt* neighbours= nullptr) {
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
		 * It returns the interactions per site of the molecule specified by its ID, with a flag to know if there is a water-water interaction in each site
		 * @param ID The ID of the molecule
		 * @param potential_matrix The matrix with the potential values, default is nullptr
		 * @param neighbours The neighbours of the molecule, default is nullptr
		 * @param R_CUT_OFF The cutoff radius, default is 5.
		 * @param V_CUT_OFF The potential cutoff, default is -12
		 * @param flag_ww A flag to know if there is a water-water interaction
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

			Sorter::cosort(sum_per_site, sum_only_water, Sorter::Order::Descending);
			return sum_only_water;
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

};

#endif
