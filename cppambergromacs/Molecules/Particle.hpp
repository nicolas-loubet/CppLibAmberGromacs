#ifndef PARTICLE_HPP
#define PARTICLE_HPP

/**
 * Version: May 2025
 * Author: Ezequiel Cuenca
 */

#include "../General/Vector.hpp"

class Particle{
	protected:
		Vector pos;
		int ID;
		Real mass;
		Real charge;
	public:
		void setPosition(Vector position){pos=position;}
		Vector getPosition() const {return(pos);}

		void setID(int ID_variable){ID=ID_variable;}
		int getID() const {return(ID);}

		void setMass(Real mass_variable){mass=mass_variable;}
		Real getMass() const {return(mass);}

		void setCharge(Real Charge){charge=Charge;}
		Real getCharge() const {return(charge);}
		
		Particle():pos({0.0,0.0,0.0}), ID(0), mass(0),charge(0) {}
		Particle (Vector position) : pos(position), ID(0), mass(0),charge(0) {}
		Particle (Vector position, int id) : pos(position), ID(id), mass(0),charge(0) {}
		Particle (Vector position, int id, Real Mass) : pos(position), ID(id), mass(Mass),charge(0) {}
		Particle (Vector position, int id, Real Mass, Real Charge) : pos(position), ID(id), mass(Mass), charge(Charge) {}

		Particle(const Particle& other) : pos(other.pos), ID(other.ID), mass(other.mass), charge(other.charge) {}

		Particle& operator=(const Particle& other) {
			if(this == &other) return *this;
			pos= other.pos;
			ID= other.ID;
			mass= other.mass;
			charge= other.charge;
			return *this;
		}

		/**
		 * It returns the distance between two particles in periodic boundary conditions
		 * @param other The other particle
		 * @param box The box size
		 * @return The distance between the two particles
		 */
		Real distanceTo(const Particle& other, const Vector box) const 
			{
			return distancePBC(pos, other.getPosition(), box);
			}
};

#endif // PARTICLE_HPP