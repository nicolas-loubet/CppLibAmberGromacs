#ifndef MOLECULE_HPP
#define MOLECULE_HPP

/**
 * Version: May 2025
 * Author: Ezequiel Cuenca
 */

#include "Atom.hpp"
#include <vector>

#define NOT_CLASSIFIED -314159265

class Molecule : public Particle{
	protected:
		Atom* atoms;
		int n_atoms;
		int classif;
		bool is_water;

		/**
		 * Calculates the total mass of the molecule
		 * @param Atoms Array of atoms
		 * @param nAtoms Number of atoms
		 * @return The total mass
		 */
		static Real totalMass(Atom* Atoms, int nAtoms){
			Real total_mass=0.0;
			for(int _i=0;_i<nAtoms;++_i){
				total_mass+=Atoms[_i].getMass();
			}
			return(total_mass);
		}

		/**
		 * Calculates the total charge of the molecule
		 * @param Atoms Array of atoms
		 * @param nAtoms Number of atoms
		 * @return The total charge
		 */
		static Real totalCharge(Atom* Atoms, int nAtoms){
			Real total_charge=0.0;
			for(int _i=0;_i<nAtoms;++_i){
				total_charge+=Atoms[_i].getCharge();
			}
			return(total_charge);
		}

	public:
		Molecule(int id, Atom* Atoms, int nAtoms, Vector pos_molecule):
			Particle(pos_molecule, id, totalMass(Atoms,nAtoms), totalCharge(Atoms,nAtoms)),
			atoms(Atoms), n_atoms(nAtoms), classif(NOT_CLASSIFIED), is_water(false) {}
		Molecule(int id, Atom* Atoms, int nAtoms) : Molecule(id, Atoms, nAtoms,
			centerOfMass(Atoms,nAtoms)) {}
		Molecule() : Molecule(0,nullptr,0) {}

		/**
		 * Copy constructor
		 */
		Molecule(const Molecule& other): Particle(other), n_atoms(other.n_atoms), atoms(nullptr) {
			if(n_atoms <= 0) return;
			atoms= new Atom[n_atoms];
			for (int i= 0; i < n_atoms; i++)
				atoms[i]= other.atoms[i];
		}

		Molecule& operator=(const Molecule& other) {
			if(this == &other) return *this;
			Particle::operator=(other);
			delete[] atoms;
			n_atoms= other.n_atoms;
			atoms= nullptr;
			if(n_atoms <= 0) return *this;
			atoms = new Atom[n_atoms];
			for(int i= 0; i < n_atoms; i++)
				atoms[i]= other.atoms[i];
			return *this;
		}
		
		/**
		 * Calculates the center of mass of the molecule
		 * @param Atoms Array of atoms
		 * @param nAtoms Number of atoms
		 * @return The center of mass
		 */
		static Vector centerOfMass(Atom* Atoms, int nAtoms){
			Vector com ={0.0,0.0,0.0};
			for(int _i=0;_i<nAtoms;++_i){
				com+=Atoms[_i].getPosition();
			}
			return(com/static_cast<Real>(nAtoms));
		}

		/**
		 * Destructor
		 */
		virtual ~Molecule() { //I need to set it virtual so I can use dynamic_cast with Water
			if(atoms) delete[] atoms;
		}

		//Getters
		Atom* getAtoms() const {return(atoms);}
		Atom& getAtom(int id) const { return(atoms[id-1]); }
		int getNAtoms() const {return(n_atoms);}
		bool isWater() const { return is_water; }

		void setClassification(const int c) { classif= c; }
		void removeClassification() { classif= NOT_CLASSIFIED; }
		int getClassification() const { return classif; }
		bool isClassified() const { return classif != NOT_CLASSIFIED; }
		void setIsWater(bool is_water) { this->is_water= is_water; }

		/**
		 * This function find the nearest atoms to the parameter, being par of the same molecule
		 * @param ID_CENTER ID of the atom center to search
		 * @param D_MAX_NEI Real that indicate maximum bond distanece to search
		 * @param bounds Vector with the size of the box for PBC
		 * @return A vector<int> with the IDs of the atoms
		 */
		inline std::vector<int> findNearbyAtoms(const int ID_CENTER, const Real D_MAX_NEI, const Vector& bounds) const {
			std::vector<int> i_nearby;
			int counter= 0;

			for(int i= 1; i <= n_atoms; i++) {
				if(i == ID_CENTER) continue;
				if(atoms[ID_CENTER-1].distanceTo(atoms[i-1], bounds) <= D_MAX_NEI)
					i_nearby.push_back(i);
			}
			return i_nearby;
		}

};

#endif
