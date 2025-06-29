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
		float sigma;
		float epsilon;
		int Z;
	public:
		void setis_Hatom(bool Hatom){is_HAtom=Hatom;}
		bool getis_Hatom() const {return(is_HAtom);}
		
		void setSigma(float s){sigma=s;}
		float getSigma() const {return(sigma);}
		void setEpsilon(float e){epsilon=e;}
		float getEpsilon() const {return(epsilon);}
		void setZ(int atom_type){Z=atom_type;is_HAtom= Z==1;}
		int getZ() const {return(Z);}

		Atom(): Particle(), is_HAtom(false), Z(0), epsilon(0), sigma(0) {}

		Atom(Vector pos, int id): Particle(pos, id, 0) {}
		
		Atom(Vector pos, int id, float mass, float charge, float e, float s, int Z):
			Particle(pos, id, mass, charge), Z(Z), epsilon(e), sigma(s) {
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

		/**
		 * Returns the atom in PDB format
		 * @param i The index of the atom
		 * @param id_molec The ID of the molecule
		 * @param atom_name The name of the atom, defaults to "UNK"
		 * @param residue_name The name of the residue, defaults to "UNK"
		 * @param chain_id The chain ID, defaults to " "
		 * @param occupancy The occupancy, defaults to 1.00f
		 * @param temp_factor The temperature factor, defaults to 0.00f
		 * @return The atom in PDB format
		 */
		std::string toPDBFormat(const int i= 0, const int id_molec= 0, const std::string& atom_name= "UNK",
			               const std::string& residue_name= "UNK", const std::string& chain_id= " ",
						   const float occupancy= 1.00f, const float temp_factor= 0.00f) const {
			const std::map<int,std::string> atomic_numbers= {
				{1,"H"},{5,"B"},{6,"C"},{7,"N"},{8,"O"},
				{9,"F"},{11,"Na"},{12,"Mg"},{13,"Al"},
				{14,"Si"},{15,"P"},{16,"S"},{17,"Cl"},
				{19,"K"},{35,"Br"},{53,"I"}
			};
			std::ostringstream oss;
			oss << "ATOM  " << std::right << std::setw(5) << i << " " << std::left << std::setw(4) << atom_name.substr(0, 4) << " ";
			oss << std::left << std::setw(3) << residue_name.substr(0, 3) << " " << chain_id.substr(0, 1);
			oss << std::right << std::setw(4) << id_molec << "    ";
			oss << std::fixed << std::setprecision(3) << std::right << std::setw(8) << pos.x;
			oss << std::fixed << std::setprecision(3) << std::right << std::setw(8) << pos.y;
			oss << std::fixed << std::setprecision(3) << std::right << std::setw(8) << pos.z;
			oss << std::fixed << std::setprecision(2) << std::right << std::setw(6) << occupancy;
			oss << std::fixed << std::setprecision(2) << std::right << std::setw(6) << temp_factor;
			oss << "          " << std::right << std::setw(2) << atomic_numbers.at(Z).substr(0,2) << "  ";
			return oss.str();
		}

};
#endif // ATOM_HPP