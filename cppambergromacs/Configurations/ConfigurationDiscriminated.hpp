#ifndef CONFIGURATION_DISCRIMINATED_HPP
#define CONFIGURATION_DISCRIMINATED_HPP

/**
 * Version: July 2025
 * Author: Nicolás Loubet
 */

#include "Configuration.hpp"

using SiteQuantity= tuple<int, Real, map<string, vector<vector<int>>>>;
using SiteEnergy  = tuple<int, Real, map<string, vector<vector<Real>>>>;

class ConfigurationDiscriminated : public virtual Configuration {
	private:
		/**
		 * Re-orders the inner vectors of a per-site map according to a site-label permutation.
		 * The permutation `labels` maps output slot i -> original site labels[i].
		 * Works for both int and Real inner vectors via template.
		 */
		template<typename T>
		static map<string, vector<vector<T>>> remapByLabels(
			const map<string, vector<vector<T>>>& src,
			const vector<int>& labels) {
			map<string, vector<vector<T>>> dst;
			for(const auto& [key, rows] : src) {
				if(rows.empty()) continue;
				for(const auto& row : rows) {
					vector<T> remapped(4);
					for(int i= 0; i < 4; i++) remapped[i]= row[labels[i]];
					dst[key].push_back(move(remapped));
				}
			}
			return dst;
		}

		// -----------------------------------------------------------------
		//  Discriminated accumulation helpers.
		//  These call the potentialWith_discriminated overloads on Water,
		//  which break the result into Coulomb and LJ and record it by name.
		// -----------------------------------------------------------------

		void accumulateDiscriminatedFromWater(const vector<Vector>& sites, vector<Real>& sum_per_site, Water& center, Water& other, Real R_CUT_OFF,
		                                      SiteQuantity& qty, SiteEnergy& coulomb, SiteEnergy& lj, string& atom_name, const int ID) {
			int i_close;
			if(!closestSiteIndex(sites, other.getPosition(), bounds, R_CUT_OFF, i_close)) return;
			sum_per_site[i_close]+= center.potentialWith_discriminated(other, bounds, qty, coulomb, lj, atom_name, i_close, ID);
		}

		void accumulateDiscriminatedFromAtom(const vector<Vector>& sites, vector<Real>& sum_per_site, Water& center, const Atom& atom, Real R_CUT_OFF,
		                                     SiteQuantity& qty, SiteEnergy& coulomb, SiteEnergy& lj, string& atom_name, const int ID) {
			int i_close;
			if(!closestSiteIndex(sites, atom.getPosition(), bounds, R_CUT_OFF, i_close)) return;
			sum_per_site[i_close]+= center.potentialWith_discriminated(atom, bounds, qty, coulomb, lj, atom_name, i_close, ID);
		}

	public:
		using Configuration::Configuration;

		/**
		 * Like getInteractionsPerSite, but also accumulates per-atom-name breakdowns of quantity, Coulomb, and LJ into the provided output vectors.
		 * Only molecules whose V4S falls within [MIN_BINS, MAX_BINS] are appended.
		 *
		 * @param ID                               Molecule ID (1-based)
		 * @param ti                               Topology info (for atom name lookup in solutes)
		 * @param atom_mapnames_quantity_complete  Output: per-atom-name neighbour counts by site
		 * @param atom_mapnames_Coulomb_complete   Output: per-atom-name Coulomb energies by site
		 * @param atom_mapnames_VLJ_complete       Output: per-atom-name LJ energies by site
		 * @param V4S_preferencia                  Output: (ID, V4S, preferred-site mask) tuples
		 * @param MIN_BINS                         Range filter for V4S, minimum
		 * @param MAX_BINS                         Range filter for V4S, maximum
		 * @param N_BINS                           Range filter for V4S, number of bins
		 */
		vector<Real> getInteractionsPerSite_discriminated(const int ID, const TopolInfo& ti, vector<SiteQuantity>& atom_mapnames_quantity_complete,
														  vector<SiteEnergy>& atom_mapnames_Coulomb_complete, vector<SiteEnergy>& atom_mapnames_VLJ_complete,
														  vector<tuple<int,Real,vector<int>>>& V4S_preferencia, const Real MIN_BINS, const Real MAX_BINS, const int N_BINS) {
			const Real R_CUT_OFF= 5.0;

			if(!getMolec(ID).isWater()) throw invalid_argument("The molecule is not a water molecule.");
			Water& molecule= *static_cast<Water*>(molecs[ID-1]);
			vector<Vector> sites= tetrahedralSites(molecule);

			SiteQuantity qty;
			SiteEnergy   coulomb_data;
			SiteEnergy   lj_data;
			vector<Real> sum_per_site(4, 0.0);

			for(int j= 0; j < N_MOLEC; j++) {
				if(j+1 == ID) continue;
				if(isIon(j)) {
					// Ions: use non-discriminated accumulator (no per-name breakdown needed)
					accumulateSiteFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(1), R_CUT_OFF);
					continue;
				}
				if(molecs[j]->isWater()) {
					Water& other= *static_cast<Water*>(molecs[j]);
					if(molecule.distanceTo(other, bounds) > R_CUT_OFF + 1.1) continue;
					string name= "OW";
					accumulateDiscriminatedFromWater(sites, sum_per_site, molecule, other,
					                                R_CUT_OFF, qty, coulomb_data, lj_data, name, ID);
				} else {
					for(int a= 1; a <= molecs[j]->getNAtoms(); a++) {
						if(molecule.distanceTo(molecs[j]->getAtom(a), bounds) >= R_CUT_OFF + 1.1) continue;
						string name= get<1>(ti.atom_type_name_charge_mass[molecs[j]->getID()-1].at(a-1));
						accumulateDiscriminatedFromAtom(sites, sum_per_site, molecule, molecs[j]->getAtom(a), R_CUT_OFF, qty, coulomb_data, lj_data, name, ID);
					}
				}
			}

			vector<int> labels= {0,1,2,3};
			Sorter::cosort(sum_per_site, labels, Sorter::Order::Ascending);

			if(sum_per_site[3] >= MIN_BINS && sum_per_site[3] <= MAX_BINS) {
				// Record which original site was V4S
				vector<int> v4s_mask(4, 0);
				v4s_mask[labels[3]]= 1;
				V4S_preferencia.emplace_back(ID, sum_per_site[3], v4s_mask);

				// Remap all per-site data to match the sorted site order
				atom_mapnames_quantity_complete.emplace_back(ID, sum_per_site[3], remapByLabels(get<2>(qty), labels));
				atom_mapnames_Coulomb_complete.emplace_back(ID, sum_per_site[3], remapByLabels(get<2>(coulomb_data), labels));
				atom_mapnames_VLJ_complete.emplace_back(ID, sum_per_site[3], remapByLabels(get<2>(lj_data), labels));
			}

			return sum_per_site;
		}
};

#endif // CONFIGURATION_DISCRIMINATED_HPP
