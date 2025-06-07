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
	static float totalMass(Atom* Atoms, int nAtoms){
		float total_mass=0.0f;
		for(int _i=0;_i<nAtoms;++_i){
			total_mass+=Atoms[_i].getMass();

		}
		return(total_mass);
	}
	public:

	Molecule(int id, Atom* Atoms, int nAtoms) : Particle(centerOfMass(Atoms,nAtoms), id, totalMass(Atoms, nAtoms)), atoms(Atoms),n_atoms(nAtoms){}

	static Vector centerOfMass(Atom* Atoms, int nAtoms){
		Vector com ={0.0f,0.0f,0.0f};
		for(int _i=0;_i<nAtoms;++_i){
			com+=Atoms[_i].getPosition();
		}
		return(com/static_cast<float>(nAtoms));
	}

	virtual ~Molecule() { //I need to set it virtual so I can use dynamic_cast with Water
		delete[] atoms;
	}
	Atom* getAtoms() const {return(atoms);}
	int getNAtoms() const {return(n_atoms);}

};

#endif
