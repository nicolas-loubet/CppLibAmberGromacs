#ifndef VECTOR_HPP
#define VECTOR_HPP

/**
 * Version: April 2025
 * Author: Nicolás Loubet
 */

#include <ostream>
#include <cmath>
#include <string>

/**
 * This class creates a mathematical Vector, with components x, y and z (3D system)
 */
class Vector {
	public:
		static constexpr double PI = 3.14159265358979323846;
		static constexpr float EPSILON= 1e-5f;
		static constexpr float RAD2DEG= 180.0f / static_cast<float>(PI);
	
		float x, y, z;

		/**
		 * Basic constructor
		 * @param v_x Is the x component of the vector (default=0)
		 * @param v_y Is the y component of the vector (default=0)
		 * @param v_z Is the z component of the vector (default=0)
		 */
		Vector(float v_x= 0.0f, float v_y= 0.0f, float v_z= 0.0f): x(v_x), y(v_y), z(v_z) {}

		/**
		 * Constructor for assigment
		 * @param v Other vector already defined
		 */
		Vector(const Vector& v): x(v.x), y(v.y), z(v.z) {}

		/**
		 * Compute the module of the vector
		 * @return module of the vector
		 */
		float magnitude() const {
			return sqrt(x*x+y*y+z*z);
		}

		/**
		 * To compare two vectors in magnitude. Otherwise, use ==
		 * @param v Other vector to compare
		 * @return true if both vectors have the same magnitude (tolerance of e-5)
		 */
		bool hasEqualMagnitude(const Vector& v) const {
			return std::fabs(magnitude() - v.magnitude()) < EPSILON;
		}

		/**
		 * OVERLOAD of operator + for two vectors
		 * @param v Other vector to sum
		 * @return The sum of each component in a new Vector object
		 */
		Vector operator +(const Vector& v) const {
			return Vector(x+v.x, y+v.y, z+v.z);
		}

		/**
		 * OVERLOAD of operator - for two vectors
		 * @param v Other vector to substract
		 * @return The substraction of each component in a new Vector object
		 */
		Vector operator -(const Vector& v) const {
			return Vector(x-v.x, y-v.y, z-v.z);
		}

		/**
		 * OVERLOAD of operator * for two vectors
		 * @param v Other vector to operate
		 * @return The scalar product of both vectors
		 */
		float operator *(const Vector& v) const {
			return x*v.x + y*v.y + z*v.z;
		}

		/**
		 * OVERLOAD of operator * for one vector and a scalar
		 * @param k A constant (scalar)
		 * @return The product of each component with the scalar k in a new Vector object
		 */
		Vector operator *(const float k) const {
			return Vector(x*k, y*k, z*k);
		}

		/**
		 * OVERLOAD of operator / for one vector and a scalar
		 * @param k A constant (scalar)
		 * @return The division of each component with the scalar k in a new Vector object
		 */
		Vector operator /(const float k) const {
			return Vector(x/k, y/k, z/k);
		}


		/**
		 * OVERLOAD of operator % for two vectors
		 * @param v Other vector to operate
		 * @return The cross (vector) product of both vectors, a new Vector object
		 */
		Vector operator %(const Vector& v) const {
			return Vector(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
		}

		/**
		 * Dot product for two vectors (this one and the parameter)
		 * @param v Other vector to operate
		 * @return The scalar product of both vectors
		 */
		float dot(const Vector& v) const {
			return x*v.x + y*v.y + z*v.z;
		}

		/**
		 * Cross product for two vectors (this one and the parameter)
		 * @param v Other vector to operate
		 * @return The cross (vector) product of both vectors, a new Vector object
		 */
		Vector cross(const Vector& v) const {
			return Vector(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
		}

		/**
		 * Normalize the vector (unit vector)
		 * @return The normalized vector in a new Vector object
		 */
		void normalize() {
			float m= magnitude();
			if(m > EPSILON) {
				x/= m;
				y/= m;
				z/= m;
			}
		}

		/**
		 * Normalize the vector (unit vector)
		 * @return The normalized vector in a new Vector object
		 */
		Vector getNormalized() const {
			float m= magnitude();
			if(m > EPSILON) return Vector(x/m, y/m, z/m);
			return Vector();
		}

		/**
		 * OVERLOAD of operator += for two vectors
		 * @param v Other vector to operate
		 * @return The sum of each component, replacing the actual object
		 */
		Vector& operator+=(const Vector& v) {
			x+= v.x; y+= v.y; z+= v.z;
			return *this;
		}

		/**
		 * OVERLOAD of operator -= for two vectors
		 * @param v Other vector to operate
		 * @return The substraction of each component, replacing the actual object
		 */
		Vector& operator-=(const Vector& v) {
			x-= v.x; y-= v.y; z-= v.z;
			return *this;
		}

		/**
		 * OVERLOAD of operator *= for a vector (this) and a scalar
		 * @param k A constant (scalar)
		 * @return The product of each component with the scalar k, replacing the actual object
		 */
		Vector& operator*=(const float k) {
			x*= k; y*= k; z*= k;
			return *this;
		}

		/**
		 * OVERLOAD of operator /= for a vector (this) and a scalar
		 * @param v Other vector to operate
		 * @return The division of each component with the scalar k, replacing the actual object
		 */
		Vector& operator/=(const float k) {
			x/= k; y/= k; z/= k;
			return *this;
		}

		/**
		 * OVERLOAD of operator %= for two vectors
		 * @param v Other vector to operate
		 * @return The cross (vector) product of both vectors, replacing the actual object
		 */
		Vector& operator%=(const Vector& v) {
			float new_x= y*v.z - z*v.y;
			float new_y= z*v.x - x*v.z;
			float new_z= x*v.y - y*v.x;
			x= new_x; y= new_y; z= new_z;
			return *this;
		}

		/**
		 * OVERLOAD of operator == for two vectors
		 * To compare two vectors in component. For comparing magnitud use equalsMagnitud()
		 * @param v Other vector to compare
		 * @return true if both vectors have the same components in x, y, and z (tolerance of e-5)
		 */
		bool operator ==(const Vector& v) const {
			return std::fabs(x - v.x) < EPSILON && std::fabs(y - v.y) < EPSILON && std::fabs(z - v.z) < EPSILON;
		}

		/**
		 * OVERLOAD of operator != for two vectors
		 * To compare two vectors in component. For comparing magnitud use !hasEqualsMagnitud()
		 * @param v Other vector to compare
		 * @return true if both vectors DO NOT have the same components in x, y, and z
		 */
		bool operator !=(const Vector& v) const {
			return !(*this==v);
		}

		/**
		 * OVERLOAD of operator < for two vectors
		 * To compare two vectors in magnitud
		 * @param v Other vector to compare
		 * @return true if this vector has a minor magnitud than v
		 */
		bool operator <(const Vector& v) const {
			return magnitude()<v.magnitude();
		}

		/**
		 * OVERLOAD of operator > for two vectors
		 * To compare two vectors in magnitud
		 * @param v Other vector to compare
		 * @return true if this vector has a major magnitud than v
		 */
		bool operator >(const Vector& v) const {
			return magnitude()>v.magnitude();
		}

		/**
		 * OVERLOAD of operator <= for two vectors
		 * To compare two vectors in magnitud
		 * @param v Other vector to compare
		 * @return true if this vector has a minor or equal magnitud than v
		 */
		bool operator <=(const Vector& v) const {
			return magnitude()<=v.magnitude();
		}

		/**
		 * OVERLOAD of operator >= for two vectors
		 * To compare two vectors in magnitud
		 * @param v Other vector to compare
		 * @return true if this vector has a major or equal magnitud than v
		 */
		bool operator >=(const Vector& v) const {
			return magnitude()>=v.magnitude();
		}

		/**
		 * Rotate this vector around a given axis by a given angle (Rodrigues' rotation formula)
		 * @param axis Unit vector representing the axis of rotation
		 * @param angle Angle of rotation in radians
		 * @return Rotated vector
		 */
		Vector rotatedAround(const Vector& axis, float angle) const {
			Vector k= axis;
			k.normalize();
			float cos_theta= std::cos(angle);
			float sin_theta= std::sin(angle);
			return (*this)*cos_theta + (k%(*this))*sin_theta + k*(k*(*this)) * (1-cos_theta);
		}

		/**
		 * Print a VMD-like representation of the vector
		 * @return A string representation of the vector in the format: {x y z}
		 */
		std::string toVMD() const {
			return "{" + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "}";
		}

		/**
		 * Mixed product (triple scalar product) between three vectors
		 * @param a First vector
		 * @param b Second vector
		 * @param c Third vector
		 * @return The scalar value of the mixed product, volume of the parallelepiped formed by the three vectors
		 */
		static float mixedProduct(const Vector& a, const Vector& b, const Vector& c) {
			return a.dot(b % c);
		}

};

/**
 * Scalar multiplication from left side. Same as operator* but with the scalar in the left side
 * @param k A constant (scalar)
 * @param v Vector to operate
 * @return The product of each component with the scalar k in a new Vector object
 */
inline Vector operator*(float k, const Vector& v) {
	return v * k;
}

/**
 * OVERLOAD of operator << for a Vector
 * Writes the components of the vector in the format: {x y z}
 * @param o Ostream to write in
 * @param v Vector to write
 * @return Ostream with the writed string
 */
inline std::ostream& operator <<(std::ostream& o, const Vector& v) {
	o << "{" << v.x << " " << v.y << " " << v.z << "}";
	return o;
}

/**
 * Compute the distance between two vectors without using Periodic Boundary Conditions (PBC)
 * @param a First vector (position)
 * @param b Second vector (position)
 * @return The scalar distance between the two points
 */
inline float distanceWithoutPBC(const Vector& a, const Vector& b) {
	float dx= a.x - b.x;
	float dy= a.y - b.y;
	float dz= a.z - b.z;
	return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Compute the minimum-image distance between two vectors using Periodic Boundary Conditions (PBC)
 * @param a First vector (position)
 * @param b Second vector (position)
 * @param box Box dimensions in x, y, z (assumed orthorhombic)
 * @return The scalar distance between the two points using minimum image convention
 */
inline float distancePBC(const Vector& a, const Vector& b, const Vector& box) {
	float dx= a.x - b.x;
	float dy= a.y - b.y;
	float dz= a.z - b.z;

	dx-= box.x * floorf(dx / box.x + 0.5f);
	dy-= box.y * floorf(dy / box.y + 0.5f);
	dz-= box.z * floorf(dz / box.z + 0.5f);

	return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Compute the 2D distance between two vectors without using PBC
 * @param a First vector (position)
 * @param b Second vector (position)
 * @return The scalar distance between the two points, in 2D (ignoring z component)
 */
inline float distance2D(const Vector& a, const Vector& b) {
	float dx= a.x - b.x;
	float dy= a.y - b.y;
	return std::sqrt(dx*dx + dy*dy);
}

/**
 * Compute the minimum image displacement vector between two positions under PBC
 * @param a First position vector
 * @param b Second position vector
 * @param box Dimensions of the periodic box (assumed orthorhombic)
 * @return Minimum image displacement vector from b to a
 */
inline Vector displacementPBC(const Vector& a, const Vector& b, const Vector& box) {
	float dx= a.x - b.x;
	float dy= a.y - b.y;
	float dz= a.z - b.z;

	dx-= box.x * floorf(dx / box.x + 0.5f);
	dy-= box.y * floorf(dy / box.y + 0.5f);
	dz-= box.z * floorf(dz / box.z + 0.5f);

	return Vector(dx,dy,dz);
}

/**
 * Compute the squared minimum image distance between two positions under PBC
 * @param a First position vector
 * @param b Second position vector
 * @param box Dimensions of the periodic box (assumed orthorhombic)
 * @return Squared minimum image distance
 */
inline float squaredDistancePBC(const Vector& a, const Vector& b, const Vector& box) {
	float dx= a.x - b.x;
	float dy= a.y - b.y;
	float dz= a.z - b.z;

	dx-= box.x * floorf(dx / box.x + 0.5f);
	dy-= box.y * floorf(dy / box.y + 0.5f);
	dz-= box.z * floorf(dz / box.z + 0.5f);

	return dx*dx + dy*dy + dz*dz;
}

/**
 * Compute the angle (in radians) between two vectors
 * @param a First vector
 * @param b Second vector
 * @return Angle between a and b in radians (0 <= angle <= PI)
 */
inline float angleBetweenRadians(const Vector& a, const Vector& b) {
	float dotProduct= a*b;
	float magA= a.magnitude();
	float magB= b.magnitude();
	
	if(magA < Vector::EPSILON || magB < Vector::EPSILON) return 0.0f;

	float cosTheta= dotProduct / (magA * magB);

	// Clamp value to [-1, 1] to avoid NaNs due to floating point errors
	cosTheta= std::fmax(-1.0f, std::fmin(1.0f, cosTheta));

	return std::acos(cosTheta);
}

/**
 * Compute the angle (in degrees) between two vectors
 * @param a First vector
 * @param b Second vector
 * @return Angle between a and b in degrees (0 <= angle <= 180)
 */
inline float angleBetweenDegrees(const Vector& a, const Vector& b) {
	return angleBetweenRadians(a,b) * Vector::RAD2DEG;
}

/**
 * Calculates the angle that 3 points form in a 3D system
 * @param c1 One of the point in the edges
 * @param c2 The center point
 * @param c3 The other point in an edge
 * @param bounds The coordinate of the last point, so the components are the width, height and length
 * @return The angle in radians formed by the 3 points
 */
inline float getAngle(Vector c1, Vector c2, Vector c3, Vector bounds) {
	float a= distancePBC(c1,c3,bounds); //Opposite to the angle
	float b= distancePBC(c1,c2,bounds);
	float c= distancePBC(c2,c3,bounds);
	//Law of cosines
	return abs(acos((pow(b,2)+pow(c,2)-pow(a,2))/(2*b*c)));
}

/**
 * Print a VMD-like representation of the line between two points
 * @param a Origin point
 * @param b Destination point
 * @return A string representation of the line in the format: draw line {x1 y1 z1} {x2 y2 z2}
 */
inline std::string lineVMD(const Vector& a, const Vector& b) {
	return "draw line " + a.toVMD() + " " + b.toVMD();
}

#endif
