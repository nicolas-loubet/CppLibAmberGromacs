#ifndef CONFIGURATIONLIPID_HPP
#define CONFIGURATIONLIPID_HPP

/**
 * Version: July 2025
 * Author: Ezequiel Cuenca
 */

#include "Configuration.hpp"

/**
 * Configuration subclass for lipid membrane systems. Provides acyl-chain order parameter calculation.
 */
class ConfigurationLipid: public Configuration {
	private:
		/**
		 * Finds the terminal CH3 carbon IDs in molecule ID_MOLEC.
		 * Excludes CH3 groups bound to nitrogen (choline head group).
		 */
		inline vector<int> findCH3(const int ID_MOLEC) const {
			vector<int> list_CH3;
			for(int i= 1; i <= getMolec(ID_MOLEC).getNAtoms(); i++) {
				if(getMolec(ID_MOLEC).getAtom(i).getZ() != 6) continue;

				vector<int> nearby= getMolec(ID_MOLEC).findNearbyAtoms(i, 1.65, bounds);
				if(nearby.empty()) continue;

				// Skip CH3 bound to N
				bool bound_to_N= false;
				for(int nb: nearby)
					if(getMolec(ID_MOLEC).getAtom(nb).getZ() == 7) { bound_to_N= true; break; }
				if(bound_to_N) continue;

				int n_H= 0;
				for(int nb: nearby)
					if(getMolec(ID_MOLEC).getAtom(nb).getZ() == 1) n_H++;
				if(n_H == 3) list_CH3.insert(list_CH3.begin(), i);
			}
			return list_CH3;
		}

		/**
		 * Traces the acyl chain from terminal CH3 upward until an oxygen is found.
		 * Returns {SN_hydrogen_count, chain} where chain is ordered from ester toward CH3, {carbon_id -> nearby_atom_ids}.
		 * SN_hydrogen_count determines whether the chain is SN1 or SN2.
		 */
		inline pair<int, vector<map<int,vector<int>>>> analyzeChain(const int ID_MOLEC, const int ID_CH3) const {
			vector<map<int,vector<int>>> chain;
			int previous_C=  ID_CH3;
			int next_C=      0;
			int SN_hydrogen= 0;

			vector<int> nearby= getMolec(ID_MOLEC).findNearbyAtoms(ID_CH3, 1.65, bounds);
			chain.insert(chain.begin(), {{previous_C, nearby}});

			// Find the first C neighbor of CH3
			for(int nb: nearby)
				if(getMolec(ID_MOLEC).getAtom(nb).getZ() == 6 && nb != previous_C) { next_C= nb; break; }

			int failsafe= 0;
			while(true) {
				nearby= getMolec(ID_MOLEC).findNearbyAtoms(next_C, 1.65, bounds);
				if(++failsafe > 25 || (int)chain.size() >= 22) break;

				bool found_O= false;
				for(int nb: nearby) {
					if(getMolec(ID_MOLEC).getAtom(nb).getZ() != 8) continue;

					// Check its neighbours for the ester carbon
					vector<int> O_nearby= getMolec(ID_MOLEC).findNearbyAtoms(nb, 1.65, bounds);
					if(O_nearby.size() < 2) continue;

					for(int ok: O_nearby) {
						if(getMolec(ID_MOLEC).getAtom(ok).getZ() != 6 || ok == next_C) continue;
						// Count H on SN carbon to distinguish SN1/SN2
						for(int snk: getMolec(ID_MOLEC).findNearbyAtoms(ok, 1.65, bounds))
							if(getMolec(ID_MOLEC).getAtom(snk).getZ() == 1) SN_hydrogen++;
						found_O= true;
						break;
					}
					if(found_O) {
						chain.insert(chain.begin(), {{next_C, nearby}});
						break;
					}
				}
				if(found_O) break;

				int next_next_C= 0;
				for(int nb: nearby)
					if(getMolec(ID_MOLEC).getAtom(nb).getZ() == 6 && nb != previous_C)
						{ next_next_C= nb; break; }

				chain.insert(chain.begin(), {{next_C, nearby}});
				previous_C= next_C;
				next_C=     next_next_C;
				if(next_C == 0) break;
			}

			return {SN_hydrogen, chain};
		}

	public:
		ConfigurationLipid(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info): Configuration(coord_reader, filename, topol_info) {}

		/**
		 * Calculates the acyl-chain order parameter for molecule ID_MOLEC.
		 * Returns a 3×22 matrix: first dimension = chain (SN1/SN2), second = per-carbon order.
		 */
		vector<vector<Real>> orderParameter(const int ID_MOLEC) {
			vector<int>          CH3_list= findCH3(ID_MOLEC);
			vector<vector<Real>> order_per_chain(3, vector<Real>(22, 0.0));

			for(int i= 0; i < (int)CH3_list.size(); i++) {
				auto [sn, chain]= analyzeChain(ID_MOLEC, CH3_list[i]);
				vector<Real> order_per_C(22, 0.0);

				for(int j= 0; j < (int)chain.size(); j++) {
					const auto& [carbon_id, atom_ids]= *chain[j].begin();
					Real order= 0.0;
					int  n_H=   0;
					for(int k: atom_ids) {
						if(getMolec(ID_MOLEC).getAtom(k).getZ() != 1) continue;
						Vector CH= getMolec(ID_MOLEC).getAtom(carbon_id).getPosition() - getMolec(ID_MOLEC).getAtom(k).getPosition();
						Real z= (CH / CH.magnitude()).z;
						order+= 0.5 * (3*z*z - 1);
						n_H++;
					}
					if(n_H == 0) n_H= 1;
					order_per_C[j]= order / n_H;
				}
				order_per_chain[sn]= order_per_C;
			}
			return order_per_chain;
		}
};

#endif // CONFIGURATIONLIPID_HPP
