#ifndef CONFIGURATION_DEFECTS_HPP
#define CONFIGURATION_DEFECTS_HPP

/**
 * Version: July 2025
 * Author: Nicolás Loubet
 */

#include "Configuration.hpp"

/**
 * Result of classifyDefect(). Describes the tetrahedral coordination state of one water molecule.
 * A molecule is DJ if it has simultaneously at least one lacking site (D3-like) and one bifurcated site (D5-like). Otherwise it is purely D3 or D5.
 */
struct DefectInfo {
	bool is_DJ= false;
	bool is_D3= false;
	bool is_D5= false;
	int  lacking_sites=     0;
	int  bifurcated_sites=  0;

	vector<Real>     sum_per_site;
	pair<Real,Real>  bifurcated_individual_potentials;
	pair<Real,Real>  bifurcated_individual_distances;
	pair<int,int>    bifurcated_individual_indices;
	Real             bifurcated_site_potential= 0.0;
	Real             lacking_site_potential=    0.0;
	Vector           lacking_site_position;
	Vector           bifurcated_site_position;

	/**
	 * Fills this struct from per-site interaction data.
	 * `ww_interactions[i]` holds the potentials of all water/ion neighbours assigned to site i that passed the V_CUT_OFF threshold.
	 */
	void chargeData(const vector<Vector>& sites, const vector<Real>& sps, const vector<vector<Real>>& ww_interactions,
	                const vector<vector<Real>>& ww_distances, const vector<vector<int>>&  ww_indices) {
		sum_per_site= sps;
		bool hasD3= false, hasD5= false;
		for(int i= 0; i < 4; i++) {
			if(ww_interactions[i].empty()) {
				hasD3= true;
				lacking_sites++;
				lacking_site_position= sites[i];
				lacking_site_potential= sps[i];
			} else if(ww_interactions[i].size() >= 2) {
				hasD5= true;
				bifurcated_sites++;
				bifurcated_site_position= sites[i];
				bifurcated_site_potential= sps[i];
				int a= 0, b= 1;
				if(ww_interactions[i][0] > ww_interactions[i][1]) swap(a, b);
				bifurcated_individual_potentials= {ww_interactions[i][a], ww_interactions[i][b]};
				bifurcated_individual_distances=  {ww_distances[i][a],    ww_distances[i][b]};
				bifurcated_individual_indices=    {ww_indices[i][a],      ww_indices[i][b]};
			}
		}
		is_DJ= hasD3 && hasD5;
		is_D3= !is_DJ && hasD3;
		is_D5= !is_DJ && hasD5;
	}
};

class ConfigurationDefects : public virtual Configuration {
	protected:
		// -----------------------------------------------------------------
		//  Per-site accumulation helpers for defect analysis.
		//  These track not just the potential sum but also a per-site list
		//  of neighbours that pass the V_CUT_OFF threshold.
		// -----------------------------------------------------------------
		void accumulateDefectSiteFromIon(Water& molecule, const Molecule& ion, const vector<Vector>& sites, Real R_CUT_OFF, Real V_CUT_OFF,
		                                 vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices,
										 vector<Real>& sum_per_site) {
			auto [idx, dist]= closestSiteIndexAndDistance(sites, ion.getPosition(), bounds);
			if(dist > R_CUT_OFF) return;
			Real pot= molecule.potentialWith(ion, bounds, special_interactions);
			sum_per_site[idx]+= pot;
			if(pot <= V_CUT_OFF) {
				ww_interactions[idx].push_back(pot);
				ww_distances[idx].push_back(molecule.distanceTo(ion, bounds));
				ww_indices[idx].push_back(ion.getID());
			}
		}

		void accumulateDefectSiteFromWater(Water& molecule, Water& other, const vector<Vector>& sites, Real R_CUT_OFF, Real V_CUT_OFF,
		                                   vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices,
										   vector<Real>& sum_per_site) {
			if(molecule.distanceTo(other, bounds) > R_CUT_OFF + 1.1) return;
			auto [idx, dist]= closestSiteIndexAndDistance(sites, other.getPosition(), bounds);
			if(dist > R_CUT_OFF) return;
			Real pot= molecule.potentialWith(other, bounds);
			sum_per_site[idx]+= pot;
			if(pot <= V_CUT_OFF) {
				ww_interactions[idx].push_back(pot);
				ww_distances[idx].push_back(molecule.distanceTo(other, bounds));
				ww_indices[idx].push_back(other.getID());
			}
		}

		void accumulateDefectSiteFromSolute(Water& molecule, const Molecule& solute, const vector<Vector>& sites, Real R_CUT_OFF, Real V_CUT_OFF,
		                                    vector<vector<Real>>& ww_interactions, vector<vector<Real>>& ww_distances, vector<vector<int>>& ww_indices,
											vector<Real>& sum_per_site) {
			vector<Real> site_pot(4, 0.0);
			vector<Real> site_min_dist(4, numeric_limits<Real>::infinity());

			for(int a= 1; a <= solute.getNAtoms(); a++) {
				if(molecule.distanceTo(solute.getAtom(a), bounds) >= R_CUT_OFF + 1.1) continue;
				auto [idx, dist]= closestSiteIndexAndDistance(sites, solute.getAtom(a).getPosition(), bounds);
				if(dist > R_CUT_OFF) continue;
				Real pot= molecule.potentialWith(solute.getAtom(a), bounds, special_interactions);
				sum_per_site[idx]+= pot;
				site_pot[idx]    += pot;
				if(dist < site_min_dist[idx]) site_min_dist[idx]= dist;
			}

			for(int i= 0; i < 4; i++) {
				if(site_pot[i] > V_CUT_OFF) continue;
				ww_interactions[i].push_back(site_pot[i]);
				ww_distances[i].push_back(site_min_dist[i]);
				ww_indices[i].push_back(solute.getID());
				break;
			}
		}

	public:
		using Configuration::Configuration;

		/**
		 * Classifies the tetrahedral coordination state of molecule ID.
		 * D3: at least one site with no qualified neighbour (V <= V_CUT_OFF).
		 * D5: at least one site with two or more qualified neighbours.
		 * DJ: simultaneously D3 and D5.
		 */
		DefectInfo classifyDefect(const int ID, const Real R_CUT_OFF= 5.0, const Real V_CUT_OFF= -12.0) {
			DefectInfo output;
			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");

			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			auto tv= Geometrics::getPerfectTetrahedron(
				molecule.getOxygen().getPosition(),
				molecule.getHydrogen_1().getPosition(),
				molecule.getHydrogen_2().getPosition(),
				bounds
			);
			vector<Vector> sites= {tv.H1, tv.H2, tv.L1, tv.L2};

			vector<vector<Real>> ww_interactions(4), ww_distances(4);
			vector<vector<int>>  ww_indices(4);
			vector<Real>         sum_per_site(4, 0.0);

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue;
				const Molecule& other_m= *molecs[j];
				if(isIon(j)) {
					accumulateDefectSiteFromIon(molecule, other_m, sites, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				} else if(molecs[j]->isWater()) {
					Water& other= *static_cast<Water*>(molecs[j]);
					accumulateDefectSiteFromWater(molecule, other, sites, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				} else {
					accumulateDefectSiteFromSolute(molecule, other_m, sites, R_CUT_OFF, V_CUT_OFF, ww_interactions, ww_distances, ww_indices, sum_per_site);
				}
			}

			Sorter::sort(sum_per_site, Sorter::Order::Ascending);
			output.chargeData(sites, sum_per_site, ww_interactions, ww_distances, ww_indices);
			return output;
		}
};

#endif // CONFIGURATION_DEFECTS_HPP
