#ifndef WATER_HPP
#define WATER_HPP

/**
 * Version: June 2025
 * Author: Nicolás Loubet
 */

#include "Molecule.hpp"

/**
 * This class creates a Water molecule object, with the 3 atoms: H2O
 */
class Water : public Molecule {
	public:
		// Getters
		Atom getOxygen() const { return atoms[0]; }
		Atom getHydrogen_1() const { return atoms[1]; }
		Atom getHydrogen_2() const { return atoms[2]; }

		/**
		 * Basic constructor
		 * @param id Number of ID to identify the molecule in the configuration
		 * @param atoms *Atom for the oxygen and the two hydrogens
		 */
		Water(const int id, Atom* atoms, const int n_atoms): Molecule(id, atoms, n_atoms, atoms[0].getPosition()) {
			is_water= true;
		}

		Water(): Molecule() {
			is_water= true;
		}

		Water(const Water& other): Molecule(other) {
			is_water= true;
		}

		Water& operator=(const Water& other) {
			if(this == &other) return *this;
			Molecule::operator=(other);
			classif= other.classif;
			return *this;
		}
		/**
		 * Destructor
		 */
		~Water()= default;

		/**
		 * Determines if this molecule is linked to the one specified with an hydrogen bond (HB)
		 * @param m *Water to the other molecule
		 * @param MAX_D_HB The maximum distance O-O that it could be considereded an HB (near 3.5A)
		 * @param MAX_A_HB The angle O-O-H that it could be considereded an HB (near 30°)
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return true if the two molecules form an HB
		 */
		bool isHB(const Water& m, const Vector& bounds, const Real MAX_D_HB= 3.5, const Real MAX_A_HB= 30) const {
			if(distanceTo(m, bounds) > MAX_D_HB) return false;

			const Real MAX_A_HB_RAD= MAX_A_HB*Vector::PI/180.;
			Atom* atoms_other= m.getAtoms();
			if(getAngle(atoms[0].getPosition(), atoms_other[0].getPosition(), atoms_other[1].getPosition(), bounds) < MAX_A_HB_RAD) return true;
			if(getAngle(atoms[0].getPosition(), atoms_other[0].getPosition(), atoms_other[2].getPosition(), bounds) < MAX_A_HB_RAD) return true;
			if(getAngle(atoms[1].getPosition(), atoms[0].getPosition(), atoms_other[0].getPosition(), bounds) < MAX_A_HB_RAD) return true;
			if(getAngle(atoms[2].getPosition(), atoms[0].getPosition(), atoms_other[0].getPosition(), bounds) < MAX_A_HB_RAD) return true;
			return false;
		}

		/**
		 * Calculates the Lennard-Jones potential energy between two atoms
		 * @param m *Particle to the other molecule to calculate the LJPotential with
		 * @param s Parameter of the Lennard-Jones potential ; s (Angstrom)
		 * @param e Parameter of the Lennard-Jones potential ; e (kJ/mol)
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return The L-J potential energy in kJ/mol between the two molecules
		 */
		Real getLJPotential(const Particle& m, const Real s, const Real e, const Vector& bounds) const {
			const Real s_over_R= s/distanceTo(m,bounds);
			return 4*e*(pow(s_over_R,12)-pow(s_over_R,6));
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m Water the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		Real potentialWith(const Water& m, const Vector& bounds) const {
			Real Vtot= getLJPotential(m, atoms[0].getSigma(), atoms[0].getEpsilon(), bounds);

			Atom* arr_other= m.getAtoms();
			for(int i= 0; i < n_atoms; i++)
				for(int j= 0; j < m.getNAtoms(); j++)
					Vtot+= atoms[i].getCoulombPotential(arr_other[j],bounds);
					
			return Vtot;
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m Molecule the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		Real potentialWith(const Molecule& m, const Vector& bounds) const {
			Atom* arr_other= m.getAtoms();
			//Combination rules: Lorentz-Berthelot
			Real s= .5*(atoms[0].getSigma() + arr_other[0].getSigma());
			Real e= sqrt(atoms[0].getEpsilon() * arr_other[0].getEpsilon());

			Real Vtot= getLJPotential(m, s, e, bounds);

			for(int i= 0; i < n_atoms; i++)
				for(int j= 0; j < m.getNAtoms(); j++)
					Vtot+= atoms[i].getCoulombPotential(arr_other[j],bounds);
					
			return Vtot;
		}

		/**
		 * Calculates the potential of this molecule with an atom
		 * @param atom Atom interacting with this molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		Real potentialWith(const Atom& atom, const Vector& bounds) {
			Real s= .5*(atoms[0].getSigma() + atom.getSigma());
			Real e= sqrt(atoms[0].getEpsilon() * atom.getEpsilon());

			Real Vtot= getLJPotential(atom, s, e, bounds);

			for(int i= 0; i < n_atoms; i++)
				Vtot+= atoms[i].getCoulombPotential(atom,bounds);
					
			return Vtot;
		}

};

#endif
