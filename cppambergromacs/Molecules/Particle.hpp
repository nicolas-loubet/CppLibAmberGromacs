#ifndef PARTICLE_HPP
#define PARTICLE_HPP

/**
 * Version: May 2025
 * Author: Ezequiel Cuenca
 */

#include "../General/Vector.hpp"

class Particle{
	private:
		Vector pos;
		int ID;
		float mass;
	public:
		void setPosition(Vector position){pos=position;}
		Vector getPosition() const {return(pos);}

		void setID(int ID_variable){ID=ID_variable;}
		int getID() const {return(ID);}

		void setMass(float mass_variable){mass=mass_variable;}
		float getMass() const {return(mass);}
		
		Particle():pos({0.0f,0.0f,0.0f}), ID(0), mass(0){}
		Particle (Vector position) : pos(position), ID(0), mass(0) {}
		Particle (Vector position, int id) : pos(position), ID(id), mass(0) {}
		Particle (Vector position, int id, float Mass) : pos(position), ID(id), mass(Mass) {}

		float distanceTo(Particle other, const Vector box) const 
			{
				return(distancePBC(pos,other.getPosition(),box));
			}
};

#endif // PARTICLE_H