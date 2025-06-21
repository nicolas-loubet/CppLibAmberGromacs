#ifndef WATER_HPP
#define WATER_HPP

/**
 * Version: April 2025
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
		static constexpr float M_PI= 3.14159f;
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
		Water(const int id, Atom* atoms): Molecule(id, atoms, 3, atoms[0].getPosition()), classif(NOT_CLASSIFIED) {}

		/**
		 * Virtual destructor
		 */
		virtual ~Water()= default;

		/**
		 * Calculates the potential of this molecule with another specified
		 * Virtual, must be redefined explicitly
		 * @param m Water the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		virtual float potentialWith(const Water& m, const Vector& bounds) const= 0;

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

			const float MAX_A_HB_RAD= MAX_A_HB*M_PI/180.;
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
		 * @param m Molecule the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(const Molecule& m, const Vector& bounds) const {
			Atom* arr_other= m.getAtoms();
			//Combination rules: Lorentz-Berthelot
			float s= .5*(atoms[0].getSigma_AO() + arr_other[0].getSigma_AO());
			float e= sqrt(atoms[0].getEpsilon_AO() * arr_other[0].getEpsilon_AO());

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
			float s= .5*(atoms[0].getSigma_AO() + atom.getSigma_AO());
			float e= sqrt(atoms[0].getEpsilon_AO() * atom.getEpsilon_AO());

			float Vtot= getLJPotential(atom, s, e, bounds);

			for(int i= 0; i < n_atoms; i++)
				Vtot+= atoms[i].getCoulombPotential(atom,bounds);
					
			return Vtot;
		}

};


/**
 * This class creates a Water molecule object, with 3 atoms: H2O (SPC/E model https://doi.org/10.1021/j100308a038)
 */
class SPCEWater : public Water {
	public:
		static constexpr float CHARGE_O= -0.8476;
		static constexpr float CHARGE_H= -CHARGE_O/2;
		static constexpr float EPSILON= 0.650;
		static constexpr float SIGMA= 3.166;

		/**
		 * Basic constructor. You must use delete(SPCEWater) after using it
		 * @param id Number of ID to identify the molecule in the configuration
		 * @param atoms *Atom for oxygen, hydrogen 1 and hydrogen 2. You can use the basic constructor Atom(Vector pos, int id)
		 */
		SPCEWater(const int id, Atom* atoms): Water(id, atoms) {
			atoms[0].setMass(15.9994);
			atoms[0].setCharge(CHARGE_O);
			atoms[0].setEpsilon_AO(EPSILON);
			atoms[0].setSigma_AO(SIGMA);

			atoms[1].setMass(1.008);
			atoms[1].setCharge(CHARGE_H);
			atoms[1].setEpsilon_AO(0);
			atoms[1].setSigma_AO(0);

			atoms[2].setMass(1.008);
			atoms[2].setCharge(CHARGE_H);
			atoms[2].setEpsilon_AO(0);
			atoms[2].setSigma_AO(0);

			n_atoms= 3;
		}

		/**
		 * Destructor. It destroys the atoms and the position of the molecule
		 */
		~SPCEWater() {/*It only calls Molecule destructor*/}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m SPCEWater the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(const Water& m, const Vector& bounds) const override {
			float Vtot= getLJPotential(m, SIGMA, EPSILON, bounds);

			Atom* arr_other= m.getAtoms();
			for(int i= 0; i < n_atoms; i++)
				for(int j= 0; j < m.getNAtoms(); j++)
					Vtot+= atoms[i].getCoulombPotential(arr_other[j],bounds);
					
			return Vtot;
		}
};

/**
 * This class creates a Water molecule object, with 3 atoms: H2O (TIP3P https://doi.org/10.1063/1.445869 )
 */
class TIP3PWater : public Water {
	public:
		static constexpr float CHARGE_O= -0.834;
		static constexpr float CHARGE_H= -CHARGE_O/2;
		static constexpr float EPSILON= 0.6364;
		static constexpr float SIGMA= 3.15061;

		/**
		 * Basic constructor. You must use delete(TIP3PWater) after using it
		 * @param id Number of ID to identify the molecule in the configuration
		 * @param atoms *Atom for oxygen, hydrogen 1 and hydrogen 2. You can use the basic constructor Atom(Vector pos, int id)
		 */
		TIP3PWater(const int id, Atom* atoms): Water(id, atoms) {
			atoms[0].setMass(15.9994);
			atoms[0].setCharge(CHARGE_O);
			atoms[0].setEpsilon_AO(EPSILON);
			atoms[0].setSigma_AO(SIGMA);

			atoms[1].setMass(1.008);
			atoms[1].setCharge(CHARGE_H);
			atoms[1].setEpsilon_AO(0);
			atoms[1].setSigma_AO(0);

			atoms[2].setMass(1.008);
			atoms[2].setCharge(CHARGE_H);
			atoms[2].setEpsilon_AO(0);
			atoms[2].setSigma_AO(0);

			n_atoms= 3;
		}

		/**
		 * Destructor. It destroys the atoms and the position of the molecule
		 */
		~TIP3PWater() {/*It only calls Molecule destructor*/}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m SPCEWater the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(const Water& m, const Vector& bounds) const override {
			float Vtot= getLJPotential(m, SIGMA, EPSILON, bounds);

			Atom* arr_other= m.getAtoms();
			for(int i= 0; i < n_atoms; i++)
				for(int j= 0; j < m.getNAtoms(); j++)
					Vtot+= atoms[i].getCoulombPotential(arr_other[j],bounds);
					
			return Vtot;
		}
};

/**
 * This class creates a Water molecule object, with 4 atoms: H2O+M (TIP4P/2005 model https://doi.org/10.1063/1.2121687)
 */
class TIP4PWater : public Water {
	public:
		static constexpr float CHARGE_MW= -1.1128;
		static constexpr float CHARGE_H= -CHARGE_MW/2;
		static constexpr float EPSILON= 0.7749;
		static constexpr float SIGMA= 3.1589;

		/**
		 * Basic constructor. You must use delete(TIP4PWater) after using it
		 * @param id Number of ID to identify the molecule in the configuration
		 * @param atoms *Atom for oxygen, hydrogen 1, hydrogen 2 and dummy atom (electrons). You can use the basic constructor Atom(Vector pos, int id)
		 */
		TIP4PWater(const int id, Atom* atoms): Water(id, atoms) {
			atoms[0].setMass(15.9994);
			atoms[0].setCharge(0);
			atoms[0].setEpsilon_AO(EPSILON);
			atoms[0].setSigma_AO(SIGMA);

			atoms[1].setMass(1.008);
			atoms[1].setCharge(CHARGE_H);
			atoms[1].setEpsilon_AO(0);
			atoms[1].setSigma_AO(0);

			atoms[2].setMass(1.008);
			atoms[2].setCharge(CHARGE_H);
			atoms[2].setEpsilon_AO(0);
			atoms[2].setSigma_AO(0);

			atoms[3].setMass(0);
			atoms[3].setCharge(CHARGE_MW);
			atoms[3].setEpsilon_AO(0);
			atoms[3].setSigma_AO(0);

			n_atoms= 4;
		}

		/**
		 * Destructor. It destroys the atoms and the position of the molecule
		 */
		~TIP4PWater() {/*It only calls Molecule destructor*/}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m SPCEWater the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(const Water& m, const Vector& bounds) const override {
			float Vtot= getLJPotential(m, SIGMA, EPSILON, bounds);

			Atom* arr_other= m.getAtoms();
			for(int i= 0; i < n_atoms; i++)
				for(int j= 0; j < m.getNAtoms(); j++)
					Vtot+= atoms[i].getCoulombPotential(arr_other[j],bounds);
					
			return Vtot;
		}

		Atom getElectrons() const { return atoms[3]; }
};

/**
 * This class creates a Water molecule object, with 5 atoms: H2O+L2 (TIP5P model https://doi.org/10.1063/1.481505) [There is a 2018 update]
 */
class TIP5PWater : public Water {
	public:
		static constexpr float CHARGE_O= -0.641114;
		static constexpr float CHARGE_H= 0.394137;
		static constexpr float CHARGE_L= -0.07358;
		static constexpr float EPSILON= 0.790;
		static constexpr float SIGMA= 3.145;

		/**
		 * Basic constructor. You must use delete(TIP5PWater) after using it
		 * @param id Number of ID to identify the molecule in the configuration
		 * @param atoms *Atom for oxygen, hydrogen 1, hydrogen 2 and 2 dummys (electrons). You can use the basic constructor Atom(Vector pos, int id)
		 */
		TIP5PWater(const int id, Atom* atoms): Water(id, atoms) {
			atoms[0].setMass(15.9994);
			atoms[0].setCharge(CHARGE_O);
			atoms[0].setEpsilon_AO(EPSILON);
			atoms[0].setSigma_AO(SIGMA);

			atoms[1].setMass(1.008);
			atoms[1].setCharge(CHARGE_H);
			atoms[1].setEpsilon_AO(0);
			atoms[1].setSigma_AO(0);

			atoms[2].setMass(1.008);
			atoms[2].setCharge(CHARGE_H);
			atoms[2].setEpsilon_AO(0);
			atoms[2].setSigma_AO(0);

			atoms[3].setMass(0);
			atoms[3].setCharge(CHARGE_L);
			atoms[3].setEpsilon_AO(0);
			atoms[3].setSigma_AO(0);

			atoms[4].setMass(0);
			atoms[4].setCharge(CHARGE_L);
			atoms[4].setEpsilon_AO(0);
			atoms[4].setSigma_AO(0);

			n_atoms= 5;
		}

		/**
		 * Destructor. It destroys the atoms and the position of the molecule
		 */
		~TIP5PWater() {/*It only calls Molecule destructor*/}

		/**
		 * Calculates the potential of this molecule with another specified
		 * @param m SPCEWater the other molecule
		 * @param bounds The coordinate of the last point, so the components are the width, height and length
		 * @return the potential energy in kJ/mol of the interaction between this two molecules
		 */
		float potentialWith(const Water& m, const Vector& bounds) const override {
			float Vtot= getLJPotential(m, SIGMA, EPSILON, bounds);

			Atom* arr_other= m.getAtoms();
			for(int i= 0; i < n_atoms; i++)
				for(int j= 0; j < m.getNAtoms(); j++)
					Vtot+= atoms[i].getCoulombPotential(arr_other[j],bounds);
					
			return Vtot;
		}

		Atom getElectrons_1() const { return atoms[3]; }
		Atom getElectrons_2() const { return atoms[4]; }
};


#endif
