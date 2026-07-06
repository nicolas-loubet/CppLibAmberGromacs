#ifndef CONFIGURATIONBULK_HPP
#define CONFIGURATIONBULK_HPP

/**
 * Version: July 2025
 * Author: Nicolás Loubet
 */

#include "ConfigurationDefects.hpp"

/**
 * Configuration subclass for water bulk systems.
 */
class ConfigurationBulk : public ConfigurationDefects {
	private:
		bool isDefect(const int ID_CENTER, const Real threshold= -12.0, const int V_index= 4) {
			return (vI(ID_CENTER, V_index) > threshold);
		}

	public:
		static constexpr int CLASSIFICATION_D_MOLECULE=  0;
		static constexpr int CLASSIFICATION_T0_MOLECULE= 1;
		static constexpr int CLASSIFICATION_T1_MOLECULE= 2;
		static constexpr int CLASSIFICATION_T2_MOLECULE= 3;
		static constexpr int CLASSIFICATION_D3_MOLECULE= 0;
		static constexpr int CLASSIFICATION_D5_MOLECULE= 1;
		static constexpr int CLASSIFICATION_TA_MOLECULE= 2;
		static constexpr int CLASSIFICATION_TB_MOLECULE= 3;

		ConfigurationBulk(CoordinateReader* coord_reader, const string& filename, TopolInfo& topol_info): Configuration(coord_reader, filename, topol_info) {}

		/**
		 * Returns true if molecule ID_CENTER is a D-type defect (vI > threshold).
		 */
		bool isD(const int ID_CENTER, const Real threshold= -12.0, const int V_index= 4) {
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false;
			if(m->getClassification() != NOT_CLASSIFIED)
				return (m->getClassification() == CLASSIFICATION_D_MOLECULE || m->getClassification() == CLASSIFICATION_D3_MOLECULE);
			if(isDefect(ID_CENTER, threshold, V_index)) {
				m->setClassification(CLASSIFICATION_D_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * Returns true if molecule ID_CENTER is D3 (under-coordinated defect).
		 */
		bool isD3(const int ID_CENTER, const Real threshold= -12.0) {
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false;
			if(m->getClassification() != NOT_CLASSIFIED)
				return (m->getClassification() == CLASSIFICATION_D_MOLECULE ||
				        m->getClassification() == CLASSIFICATION_D3_MOLECULE);
			if(isDefect(ID_CENTER, threshold, 4)) {
				m->setClassification(CLASSIFICATION_D3_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * Returns true if molecule ID_CENTER is D5 (over-coordinated defect).
		 */
		bool isD5(const int ID_CENTER, const Real threshold= -12.0) {
			Water* m= dynamic_cast<Water*>(molecs[ID_CENTER-1]);
			if(m == nullptr) return false;
			if(m->getClassification() != NOT_CLASSIFIED)
				return (m->getClassification() == CLASSIFICATION_D5_MOLECULE);
			if(!isDefect(ID_CENTER, threshold, 5)) {
				m->setClassification(CLASSIFICATION_D5_MOLECULE);
				return true;
			}
			return false;
		}

		/**
		 * Returns true if molecule ID_CENTER is a DX-type defect (D3 or D5).
		 */
		bool isDX(const int ID_CENTER, const Real threshold= -12.0) {
			return isD3(ID_CENTER, threshold) || isD5(ID_CENTER, threshold);
		}

		/**
		 * Allocates the lower-triangle potential matrix used by getNeighboursByPotential.
		 * Caller is responsible for freeing with deletePotentialMatrix.
		 */
		Real** createPotentialMatrix() {
			Real** pm= new Real*[N_MOLEC];
			for(int i= 1; i < N_MOLEC; i++) {
				pm[i]= new Real[i];
				for(int j= 0; j < i; j++) pm[i][j]= NOT_CLASSIFIED;
			}
			return pm;
		}

		void deletePotentialMatrix(Real** pm) {
			for(int i= 1; i < N_MOLEC; i++) delete[] pm[i];
			delete[] pm;
		}

		/**
		 * Classifies all water molecules as D/T0/T1/T2 (JCP 2023 scheme).
		 * T0 = tetrahedral neighbour of a D; T1 = tetrahedral neighbour of a T0; T2 = bulk-like.
		 */
		void classifyMolecules(const int V_index= 4, const Real threshold= -12.0) {
			Real** pm= createPotentialMatrix();

			for(int i= 0; i < N_MOLEC; i++) {
				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w == nullptr) continue;
				vector<Real> pots; vector<int> ids;
				getNeighboursByPotential(w, pots, ids, pm);
				Sorter::sort(pots, Sorter::Order::Ascending);
				w->setClassification(pots[V_index-1] > threshold
					? CLASSIFICATION_D_MOLECULE : CLASSIFICATION_T2_MOLECULE);
			}

			for(int t= CLASSIFICATION_T0_MOLECULE; t <= CLASSIFICATION_T1_MOLECULE; t++)
				for(int i= 0; i < N_MOLEC; i++) {
					Water* w= dynamic_cast<Water*>(molecs[i]);
					if(w == nullptr || w->getClassification() != CLASSIFICATION_T2_MOLECULE) continue;
					vector<Real> pots; vector<int> ids;
					getNeighboursByPotential(w, pots, ids, pm);
					Sorter::cosort(pots, ids, Sorter::Order::Ascending);
					for(int iv= 0; iv < V_index; iv++) {
						Water* neigh= dynamic_cast<Water*>(molecs[ids[iv]-1]);
						if(neigh->getClassification() == t-1) {
							w->setClassification(t);
							break;
						}
					}
				}

			deletePotentialMatrix(pm);
		}

		/**
		 * Classifies all water molecules as D3/D5/TA/TB (PRE 2024 scheme).
		 * TA = tetrahedral neighbour of a D3 or D5; TB = bulk-like.
		 */
		void classifyMolecules_includePentacoordinated(const int V_index= 4, const Real threshold= -12.0) {
			Real** pm= createPotentialMatrix();

			for(int i= 0; i < N_MOLEC; i++) {
				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w == nullptr) continue;
				vector<Real> pots; vector<int> ids;
				getNeighboursByPotential(w, pots, ids, pm);
				Sorter::cosort(pots, ids, Sorter::Order::Ascending);
				if(pots[V_index-1] > threshold)
					w->setClassification(CLASSIFICATION_D3_MOLECULE);
				else if(pots[V_index] < threshold)
					w->setClassification(CLASSIFICATION_D5_MOLECULE);
				else
					w->setClassification(CLASSIFICATION_TB_MOLECULE);
			}

			for(int i= 0; i < N_MOLEC; i++) {
				Water* w= dynamic_cast<Water*>(molecs[i]);
				if(w == nullptr || w->getClassification() != CLASSIFICATION_TB_MOLECULE) continue;
				vector<Real> pots; vector<int> ids;
				getNeighboursByPotential(w, pots, ids, pm);
				Sorter::cosort(pots, ids, Sorter::Order::Ascending);
				for(int iv= 0; iv < V_index; iv++) {
					Water* neigh= dynamic_cast<Water*>(molecs[ids[iv]-1]);
					int cls= neigh->getClassification();
					if(cls == CLASSIFICATION_D3_MOLECULE || cls == CLASSIFICATION_D5_MOLECULE) {
						w->setClassification(CLASSIFICATION_TA_MOLECULE);
						break;
					}
				}
			}

			deletePotentialMatrix(pm);
		}

		/**
		 * Calculates the Tanaka ζ parameter for molecule m.
		 * ζ = min(d_HB) - max(d_nHB), where HB neighbours satisfy the H-bond criteria.
		 */
		Real Tanaka(Water* m, const Real MAX_D_HB= 3.5, const Real MAX_A_HB= 30.0) {
			const Real MAX_D_ANALYSIS= 6.0;
			vector<Real> ls_HB, ls_nHB;

			for(int i= 0; i < N_MOLEC; i++) {
				if(i+1 == m->getID()) continue;
				Real d= m->distanceTo(getMolec(i+1), bounds);
				if(d > MAX_D_ANALYSIS) continue;
				Water* w2= dynamic_cast<Water*>(molecs[i]);
				if(w2 == nullptr) continue;
				(m->isHB(*w2, bounds, MAX_D_HB, MAX_A_HB) ? ls_HB : ls_nHB).push_back(d);
			}

			if(ls_HB.empty())  throw runtime_error("Tanaka: no HB neighbours for Water "  + to_string(m->getID()));
			if(ls_nHB.empty()) throw runtime_error("Tanaka: no nHB neighbours for Water " + to_string(m->getID()));
			return *min_element(ls_HB.begin(),  ls_HB.end()) - *max_element(ls_nHB.begin(), ls_nHB.end());
		}

		/**
		 * Calculates the Local Structure Index for molecule id. (https://doi.org/10.1039/C1CP22076D)
		 */
		Real LSI(const int id) {
			const Real R_MAX= 3.7;
			vector<Real> distances;
			Real dist_peripheral= bounds.x + bounds.y + bounds.z;

			for(int j= 1; j <= N_MOLEC; j++) {
				if(j == id) continue;
				Real d= getMolec(id).distanceTo(getMolec(j), bounds);
				if(d > R_MAX) {
					if(d < dist_peripheral) dist_peripheral= d;
					continue;
				}
				distances.push_back(d);
			}
			distances.push_back(dist_peripheral); // first molecule outside R_MAX

			Sorter::sort(distances, Sorter::Order::Ascending);
			int N= (int)distances.size() - 1;

			Real sum_delta= 0.0, sum_sq= 0.0;
			for(int j= 0; j < N; j++) {
				Real dj= distances[j+1] - distances[j];
				sum_delta+= dj;
				sum_sq   += dj * dj;
			}
			Real mean= sum_delta / N;
			return sum_sq / N - mean * mean;
		}
};

#endif // CONFIGURATIONBULK_HPP
