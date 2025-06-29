#ifndef WATER_HPP
#define WATER_HPP

/**
 * Version: June 2025
 * Author: Nicolás Loubet
 */

#include "Molecule.hpp"

#define NOT_CLASSIFIED -314159265

/**
 * This class creates a Water molecule object, with the 3 atoms: H2O
 */
class Water : public Molecule {
	protected:
		int classif;

	public:
		//Setters and getters
		void setClassification(const int c) { classif= c; }
		void removeClassification() { classif= NOT_CLASSIFIED; }
		int getClassification() const { return classif; }
		bool isClassified() const { return classif != NOT_CLASSIFIED; }

		Atom getOxygen() const { return atoms[0]; }
		Atom getHydrogen_1() const { return atoms[1]; }
		Atom getHydrogen_2() const { return atoms[2]; }

		/**
		 * Basic constructor
		 * @param id Number of ID to identify the molecule in the configuration
		 * @param atoms *Atom for the oxygen and the two hydrogens
		 */
		Water(const int id, Atom* atoms, const int n_atoms): Molecule(id, atoms, n_atoms, atoms[0].getPosition()), classif(NOT_CLASSIFIED) {}

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
		bool isHB(const Water& m, const Vector& bounds, const float MAX_D_HB= 3.5, const float MAX_A_HB= 30) const {
			if(distanceTo(m, bounds) > MAX_D_HB) return false;

			const float MAX_A_HB_RAD= MAX_A_HB*Vector::PI/180.;
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
		float getLJPotential(const Particle& m, const float s, const float e, const Vector& bounds) const {
			const float s_over_R= s/distanceTo(m,bounds);
			return 4*e*(pow(s_over_R,12)-pow(s_over_R,6));
		}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m Water the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(const Water& m, const Vector& bounds) const {
			float Vtot= getLJPotential(m, atoms[0].getSigma(), atoms[0].getEpsilon(), bounds);

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
		float potentialWith(const Molecule& m, const Vector& bounds) const {
			Atom* arr_other= m.getAtoms();
			//Combination rules: Lorentz-Berthelot
			float s= .5*(atoms[0].getSigma() + arr_other[0].getSigma());
			float e= sqrt(atoms[0].getEpsilon() * arr_other[0].getEpsilon());

			float Vtot= getLJPotential(m, s, e, bounds);

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
		float potentialWith(const Atom& atom, const Vector& bounds) {
			float s= .5*(atoms[0].getSigma() + atom.getSigma());
			float e= sqrt(atoms[0].getEpsilon() * atom.getEpsilon());

			float Vtot= getLJPotential(atom, s, e, bounds);

			for(int i= 0; i < n_atoms; i++)
				Vtot+= atoms[i].getCoulombPotential(atom,bounds);
					
			return Vtot;
		}

};

#endif
