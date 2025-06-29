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
		float mass;
		float charge;
	public:
		void setPosition(Vector position){pos=position;}
		Vector getPosition() const {return(pos);}

		void setID(int ID_variable){ID=ID_variable;}
		int getID() const {return(ID);}

		void setMass(float mass_variable){mass=mass_variable;}
		float getMass() const {return(mass);}

		void setCharge(float Charge){charge=Charge;}
		float getCharge() const {return(charge);}
		
		Particle():pos({0.0f,0.0f,0.0f}), ID(0), mass(0),charge(0) {}
		Particle (Vector position) : pos(position), ID(0), mass(0),charge(0) {}
		Particle (Vector position, int id) : pos(position), ID(id), mass(0),charge(0) {}
		Particle (Vector position, int id, float Mass) : pos(position), ID(id), mass(Mass),charge(0) {}
		Particle (Vector position, int id, float Mass, float Charge) : pos(position), ID(id), mass(Mass), charge(Charge) {}

		/**
		 * It returns the distance between two particles in periodic boundary conditions
		 * @param other The other particle
		 * @param box The box size
		 * @return The distance between the two particles
		 */
		float distanceTo(const Particle& other, const Vector box) const 
			{
			return distancePBC(pos, other.getPosition(), box);
			}
};

#endif // PARTICLE_H