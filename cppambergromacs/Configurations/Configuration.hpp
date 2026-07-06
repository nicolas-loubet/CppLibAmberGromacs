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
#include <array>

using namespace std;

/**
 * This class creates a Configuration object with an array of molecule pointers.
 * Extended by ConfigurationDefects and ConfigurationDiscriminated (included from CppLibAmberGromacs.hpp).
 */
class Configuration {
	protected:
		Molecule** molecs;
		int        N_MOLEC;
		Vector     bounds;
		map<pair<string,string>, pair<Real,Real>> special_interactions;

		// =====================================================================
		//  HOT-PATH HELPERS
		//  Called inside the inner loop over all molecules.
		// =====================================================================

		/**
		 * Returns the index (0-3) of the tetrahedral site closest to `pos`.
		 * Returns true if that distance is within R_CUT_OFF.
		 */
		inline bool closestSiteIndex(const vector<Vector>& sites, const Vector& pos, const Vector& bnd, Real R_CUT_OFF, int& i_close) const {
			i_close= 0;
			Real d_close= distancePBC(sites[0], pos, bnd);
			for(int i= 1; i < 4; i++) {
				Real d= distancePBC(sites[i], pos, bnd);
				if(d < d_close) { i_close= i; d_close= d; }
			}
			return (d_close <= R_CUT_OFF);
		}

		/**
		 * Returns both the closest-site index and its distance. Used in defect analysis.
		 */
		inline pair<int,Real> closestSiteIndexAndDistance(const vector<Vector>& sites, const Vector& pos, const Vector& bnd) const {
			int idx= 0;
			Real d_close= distancePBC(sites[0], pos, bnd);
			for(int i= 1; i < 4; i++) {
				Real d= distancePBC(sites[i], pos, bnd);
				if(d < d_close) { idx= i; d_close= d; }
			}
			return {idx, d_close};
		}

		/** Accumulates the potential of `atom` (ion or solute atom) onto the closest site. */
		inline void accumulateSiteFromAtom(const vector<Vector>& sites, vector<Real>& sum_per_site, Water& center, const Atom& atom, Real R_CUT_OFF) {
			int i_close;
			if(closestSiteIndex(sites, atom.getPosition(), bounds, R_CUT_OFF, i_close))
				sum_per_site[i_close]+= center.potentialWith(atom, bounds, special_interactions);
		}

		/** Accumulates the potential of `other` water onto the closest site. */
		inline void accumulateSiteFromWater(const vector<Vector>& sites, vector<Real>& sum_per_site, Water& center, Water& other, Real R_CUT_OFF) {
			int i_close;
			if(closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close))
				sum_per_site[i_close]+= center.potentialWith(other, bounds);
		}

		/** Accumulates onto both `sum_per_site` and `sum_water_only` (water-only variant). */
		inline void accumulateSiteFromWater_waterOnly(const vector<Vector>& sites, vector<Real>& sum_per_site, vector<Real>& sum_water_only,
		                                              Water& center, Water& other, Real R_CUT_OFF) {
			int i_close;
			if(!closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close)) return;
			Real pot= center.potentialWith(other, bounds);
			sum_per_site[i_close]  += pot;
			sum_water_only[i_close]+= pot;
		}

		/** Accumulates potential, and flags the site if pot <= V_CUT_OFF. */
		inline void accumulateSiteFromWater_flagged(const vector<Vector>& sites, vector<Real>& sum_per_site, Water& center, Water& other,
		                                            Real R_CUT_OFF, array<bool,4>& ww_flag, Real V_CUT_OFF) {
			int i_close;
			if(!closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close)) return;
			Real pot= center.potentialWith(other, bounds);
			sum_per_site[i_close]+= pot;
			ww_flag[i_close]= ww_flag[i_close] || (pot <= V_CUT_OFF);
		}

		/** Returns the tetrahedral sites for a water molecule. */
		inline vector<Vector> tetrahedralSites(const Water& w) const {
			return Geometrics::getPerfectTetrahedron(
				w.getOxygen().getPosition(),
				w.getHydrogen_1().getPosition(),
				w.getHydrogen_2().getPosition(),
				bounds
			).toVector();
		}

		/** Returns true if molecule at index j is an ion (single atom, |charge| >= 1). */
		inline bool isIon(int j) const {
			return molecs[j]->getNAtoms() == 1 && (molecs[j]->getCharge() >= 1.0 || molecs[j]->getCharge() <= -1.0);
		}

	public:
		int             getNMolec()          const { return N_MOLEC; }
		const Molecule& getMolec(int id)     const { return *molecs[id-1]; }
		Molecule*       getMolec_ptr(int id) const { return molecs[id-1]; }
		Vector          getBounds()          const { return bounds; }

		const Water& getWater(int id) const {
			if(molecs[id-1]->isWater()) return *static_cast<Water*>(molecs[id-1]);
			throw runtime_error("The molecule is not a water");
		}

		// =====================================================================
		//  CONSTRUCTION / DESTRUCTION
		// =====================================================================

		Configuration(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info) {
			N_MOLEC= topol_info.num_molecules;
			molecs=  new Molecule*[N_MOLEC];
			if(!coord_reader->readCoordinates(filename, topol_info, molecs, bounds))
				throw runtime_error("Failed to read coordinates: " + filename);
			special_interactions= topol_info.special_interaction;
		}

		Configuration(const Configuration& other): N_MOLEC(other.N_MOLEC), bounds(other.bounds), molecs(nullptr), special_interactions(other.special_interactions) {
			molecs= new Molecule*[N_MOLEC];
			for(int i= 0; i < N_MOLEC; i++)
				molecs[i]= other.molecs[i] ? new Molecule(*other.molecs[i]) : nullptr;
		}

		Configuration(): N_MOLEC(0), bounds(Vector(0,0,0)), molecs(nullptr) {}

		Configuration& operator=(const Configuration& other) {
			if(this == &other) return *this;
			for(int i= 0; i < N_MOLEC; i++) delete molecs[i];
			delete[] molecs;
			N_MOLEC= other.N_MOLEC;
			bounds= other.bounds;
			special_interactions= other.special_interactions;
			molecs= new Molecule*[N_MOLEC];
			for(int i= 0; i < N_MOLEC; i++)
				molecs[i]= other.molecs[i] ? new Molecule(*other.molecs[i]) : nullptr;
			return *this;
		}

		~Configuration() {
			for(int i= 0; i < N_MOLEC; i++) delete molecs[i];
			delete[] molecs;
			molecs= nullptr;
		}

		// =====================================================================
		//  NEIGHBOUR / POTENTIAL UTILITIES
		// =====================================================================

		/**
		 * Returns IDs of all molecules within D_MAX_NEI of ID_CENTER.
		 */
		vector<int> findNearby(const int ID_CENTER, const Real D_MAX_NEI) const {
			vector<int> nearby;
			nearby.reserve(64);
			for(int i= 1; i <= N_MOLEC; i++) {
				if(i == ID_CENTER) continue;
				if(getMolec(ID_CENTER).distanceTo(getMolec(i), bounds) <= D_MAX_NEI)
					nearby.push_back(i);
			}
			return nearby;
		}

		/**
		 * Fills `pots` and `ids` with the potential and ID of every water within MAX_R_V4 of `m`.
		 */
		void getNeighboursByPotential(Water* m, vector<Real>& pots, vector<int>& ids, Real** potential_matrix) {
			const Real MAX_R_V4= 5.5;
			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue;
				if(!molecs[i]->isWater()) continue;
				if(m->distanceTo(getMolec(i+1), bounds) > MAX_R_V4) continue;

				Water* w2= static_cast<Water*>(molecs[i]);
				int lo= i < m->getID()-1 ? i : m->getID()-1;
				int hi= i > m->getID()-1 ? i : m->getID()-1;
				if(potential_matrix[hi][lo] == NOT_CLASSIFIED)
					potential_matrix[hi][lo]= m->potentialWith(*w2, bounds);
				pots.push_back(potential_matrix[hi][lo]);
				ids.push_back(i+1);
			}
		}

		/**
		 * Returns all water-water potentials within MAX_V4, sorted ascending.
		 */
		vector<Real> getVList(const int ID_CENTER, const Real MAX_V4= 5.5) {
			if(!getMolec(ID_CENTER).isWater()) return {};
			Water& w1= *static_cast<Water*>(molecs[ID_CENTER-1]);
			vector<Real> ls;
			ls.reserve(64);
			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID_CENTER || !molecs[j]->isWater()) continue;
				if(getMolec(ID_CENTER).distanceTo(getMolec(j+1), bounds) > MAX_V4) continue;
				ls.push_back(w1.potentialWith(*static_cast<Water*>(molecs[j]), bounds));
			}
			Sorter::sort(ls, Sorter::Order::Ascending);
			return ls;
		}

		/**
		 * Returns the V_index-th more negative water-water potential for ID_CENTER (1-based).
		 */
		Real vI(const int ID_CENTER, const int V_index= 4) {
			return getVList(ID_CENTER)[V_index-1];
		}

		// =====================================================================
		//  INTERACTIONS PER SITE
		// =====================================================================

		/**
		 * Returns the 4 site-resolved interaction energies for molecule ID, sorted ascending.
		 * @param R_CUT_OFF Cutoff radius (default 5.0 Å)
		 * @param labels    Optional: receives original site indices if non-null (new int[4]).
		 */
		vector<Real> getInteractionsPerSite(const int ID, const Real R_CUT_OFF= 5.0, int* labels= nullptr) {
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			vector<Vector> sites= tetrahedralSites(molecule);
			vector<Real> sum_per_site(4, 0.0);

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue;
				if(isIon(j)) {
					accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(1), R_CUT_OFF);
					continue;
				}
				if(molecs[j]->isWater()) {
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other, bounds) > R_CUT_OFF + 1.1) continue;
					accumulateSiteFromWater(sites, sum_per_site, molecule, other, R_CUT_OFF);
				} else {
					for(int a= 1; a <= molecs[j]->getNAtoms(); a++)
						if(molecule.distanceTo(molecs[j]->getAtom(a), bounds) < R_CUT_OFF + 1.1)
							accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(a), R_CUT_OFF);
				}
			}

			if(labels == nullptr) {
				Sorter::sort(sum_per_site, Sorter::Order::Ascending);
			} else {
				vector<int> idx= {0,1,2,3};
				Sorter::cosort(sum_per_site, idx, Sorter::Order::Ascending);
				for(int i= 0; i < 4; i++) labels[i]= idx[i];
			}
			return sum_per_site;
		}

		/**
		 * Like getInteractionsPerSite, but also sets flag_ww=true if every site has at least one water-water interaction with pot <= V_CUT_OFF.
		 */
		vector<Real> getInteractionsPerSite(const int ID, bool& flag_ww, const Real V_CUT_OFF= -12.0, const Real R_CUT_OFF= 5.0) {
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			vector<Vector> sites= tetrahedralSites(molecule);
			vector<Real> sum_per_site(4, 0.0);
			array<bool,4> ww_flag= {false, false, false, false};

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue;
				if(isIon(j)) {
					accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(1), R_CUT_OFF);
					continue;
				}
				if(molecs[j]->isWater()) {
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other, bounds) > R_CUT_OFF + 1.1) continue;
					accumulateSiteFromWater_flagged(sites, sum_per_site, molecule, other, R_CUT_OFF, ww_flag, V_CUT_OFF);
				} else {
					for(int a= 1; a <= molecs[j]->getNAtoms(); a++)
						if(molecule.distanceTo(molecs[j]->getAtom(a), bounds) < R_CUT_OFF + 1.1)
							accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(a), R_CUT_OFF);
				}
			}

			flag_ww= ww_flag[0] && ww_flag[1] && ww_flag[2] && ww_flag[3];
			Sorter::sort(sum_per_site, Sorter::Order::Ascending);
			return sum_per_site;
		}

		/**
		 * Returns both the full site interactions and the water-only contribution, co-sorted so water-only follows the same site ordering as the full result.
		 */
		pair<vector<Real>,vector<Real>> getInteractionsPerSite_waterOnly(
			const int ID, const Real V_CUT_OFF= -12.0, const Real R_CUT_OFF= 5.0) {
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			vector<Vector> sites= tetrahedralSites(molecule);
			vector<Real> sum_per_site(4, 0.0);
			vector<Real> sum_water_only(4, 0.0);

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue;
				if(isIon(j)) {
					accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(1), R_CUT_OFF);
					continue;
				}
				if(molecs[j]->isWater()) {
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other, bounds) > R_CUT_OFF + 1.1) continue;
					accumulateSiteFromWater_waterOnly(sites, sum_per_site, sum_water_only, molecule, other, R_CUT_OFF);
				} else {
					for(int a= 1; a <= molecs[j]->getNAtoms(); a++)
						if(molecule.distanceTo(molecs[j]->getAtom(a), bounds) < R_CUT_OFF + 1.1)
							accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(a), R_CUT_OFF);
				}
			}

			vector<int> idx= {0,1,2,3};
			Sorter::cosort(sum_per_site, idx, Sorter::Order::Ascending);
			return {   sum_per_site, {sum_water_only[idx[0]], sum_water_only[idx[1]], sum_water_only[idx[2]], sum_water_only[idx[3]]}   };
		}

		/**
		 * Returns the i_V-th site interaction energy for molecule ID (1-based).
		 */
		Real v_4S(const int ID, const int i_V= 4, const Real R_CUT_OFF= 5.0) {
			return getInteractionsPerSite(ID, R_CUT_OFF)[i_V-1];
		}
};

#endif // CONFIGURATION_HPP
