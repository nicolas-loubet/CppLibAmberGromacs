#ifndef ATOM_HPP
#define ATOM_HPP

/**
 * Version: May 2025
 * Author: Ezequiel Cuenca
 */

#include "Particle.hpp"

class Atom : public Particle{
	protected:
		static constexpr float K_COULOMB= 1389.35458; //kJ/mol to e-2
		bool is_HAtom;
		float s_AO;
		float e_AO;
		int Z;
	public:
		void setis_Hatom(bool Hatom){is_HAtom=Hatom;}
		bool getis_Hatom() const {return(is_HAtom);}
		
		void setSigma_AO(float s_AOxigen){s_AO=s_AOxigen;}
		float getSigma_AO() const {return(s_AO);}
		void setEpsilon_AO(float e_AOxigen){e_AO=e_AOxigen;}
		float getEpsilon_AO() const {return(e_AO);}
		void setZ(int atom_type){Z=atom_type;is_HAtom= Z==1;}
		int getZ() const {return(Z);}

		Atom(): Particle(), is_HAtom(false), Z(0), e_AO(0), s_AO(0) {}

		Atom(Vector pos, int id): Particle(pos, id, 0) {}
		
		Atom(Vector pos, int id, float mass, float charge, float e_AO, float s_AO, int Z):
			Particle(pos, id, mass, charge), Z(Z), e_AO(e_AO), s_AO(s_AO) {
				is_HAtom= Z==1;
			}

		/**
		 * Calculates the electrostatic potential energy between two atoms
		 * @param a Atom the other atom
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return The potential energy in kJ/mol between the two atoms
		 */
		float getCoulombPotential(Atom a, Vector bounds) const {
			return (K_COULOMB*charge*a.getCharge()) / distanceTo(a,bounds);
		}

};
#endif // ATOM_HPP