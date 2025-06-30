#ifndef MOLECULE_HPP
#define MOLECULE_HPP

/**
 * Version: May 2025
 * Author: Ezequiel Cuenca
 */

#include "Atom.hpp"

class Molecule : public Particle{
	protected:
		Atom* atoms;
		int n_atoms;

		/**
		 * Calculates the total mass of the molecule
		 * @param Atoms Array of atoms
		 * @param nAtoms Number of atoms
		 * @return The total mass
		 */
		static float totalMass(Atom* Atoms, int nAtoms){
			float total_mass=0.0f;
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
		static float totalCharge(Atom* Atoms, int nAtoms){
			float total_charge=0.0f;
			for(int _i=0;_i<nAtoms;++_i){
				total_charge+=Atoms[_i].getCharge();
			}
			return(total_charge);
		}

	public:
		Molecule(int id, Atom* Atoms, int nAtoms, Vector pos_molecule) : Particle(pos_molecule, id, totalMass(Atoms,nAtoms), totalCharge(Atoms,nAtoms)), atoms(Atoms), n_atoms(nAtoms){}
		Molecule(int id, Atom* Atoms, int nAtoms) : Molecule(id, Atoms, nAtoms, centerOfMass(Atoms,nAtoms)) {}
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
			Vector com ={0.0f,0.0f,0.0f};
			for(int _i=0;_i<nAtoms;++_i){
				com+=Atoms[_i].getPosition();
			}
			return(com/static_cast<float>(nAtoms));
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
};

#endif
