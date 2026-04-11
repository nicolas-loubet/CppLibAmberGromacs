#ifndef VECTOR_HPP
#define VECTOR_HPP

/**
 * Version: April 2026
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
		static constexpr Real EPSILON= 1e-5;
		static constexpr Real RAD2DEG= 180.0 / static_cast<Real>(PI);
	
		Real x, y, z;

		/**
		 * Basic constructor
		 * @param v_x Is the x component of the vector (default=0)
		 * @param v_y Is the y component of the vector (default=0)
		 * @param v_z Is the z component of the vector (default=0)
		 */
		Vector(Real v_x= 0.0, Real v_y= 0.0, Real v_z= 0.0): x(v_x), y(v_y), z(v_z) {}

		/**
		 * Constructor for assigment
		 * @param v Other vector already defined
		 */
		Vector(const Vector& v): x(v.x), y(v.y), z(v.z) {}

		/**
		 * Compute the module of the vector
		 * @return module of the vector
		 */
		Real magnitude() const {
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
		Real operator *(const Vector& v) const {
			return x*v.x + y*v.y + z*v.z;
		}

		/**
		 * OVERLOAD of operator * for one vector and a scalar
		 * @param k A constant (scalar)
		 * @return The product of each component with the scalar k in a new Vector object
		 */
		Vector operator *(const Real k) const {
			return Vector(x*k, y*k, z*k);
		}

		/**
		 * OVERLOAD of operator / for one vector and a scalar
		 * @param k A constant (scalar)
		 * @return The division of each component with the scalar k in a new Vector object
		 */
		Vector operator /(const Real k) const {
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
		Real dot(const Vector& v) const {
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
			Real m= magnitude();
			if(m > EPSILON) { x/= m; y/= m; z/= m; }
		}

		/**
		 * Normalize the vector (unit vector)
		 * @return The normalized vector in a new Vector object
		 */
		Vector getNormalized() const {
			Real m= magnitude();
			if(m > EPSILON) return Vector(x/m, y/m, z/m);
			return Vector();
		}

		/**
		 * Compute the volume of the parallelepiped formed by the three vectors
		 * @return The volume of the parallelepiped
		 */
		Real volumeBox() const {
			return x*y*z;
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
		Vector& operator*=(const Real k) {
			x*= k; y*= k; z*= k;
			return *this;
		}

		/**
		 * OVERLOAD of operator /= for a vector (this) and a scalar
		 * @param v Other vector to operate
		 * @return The division of each component with the scalar k, replacing the actual object
		 */
		Vector& operator/=(const Real k) {
			x/= k; y/= k; z/= k;
			return *this;
		}

		/**
		 * OVERLOAD of operator %= for two vectors
		 * @param v Other vector to operate
		 * @return The cross (vector) product of both vectors, replacing the actual object
		 */
		Vector& operator%=(const Vector& v) {
			Real new_x= y*v.z - z*v.y;
			Real new_y= z*v.x - x*v.z;
			Real new_z= x*v.y - y*v.x;
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
		Vector rotatedAround(const Vector& axis, Real angle) const {
			if(angle==0.0 || axis.magnitude()==0.0) return *this;
			Vector k= axis;
			k.normalize();
			Real cos_theta= std::cos(angle);
			Real sin_theta= std::sin(angle);
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
		static Real mixedProduct(const Vector& a, const Vector& b, const Vector& c) {
			return a.dot(b % c);
		}

};

/**
 * Scalar multiplication from left side. Same as operator* but with the scalar in the left side
 * @param k A constant (scalar)
 * @param v Vector to operate
 * @return The product of each component with the scalar k in a new Vector object
 */
inline Vector operator*(Real k, const Vector& v) {
	return v * k;
}

/**
 * OVERLOAD of operator << for a Vector
 * Writes the components of the vector in the format: {x y z}
 * @param o Ostream to write in
 * @param v Vector to write
 * @return Ostream with the writed string
 */
inline std::ostream& operator<<(std::ostream& o, const Vector& v) {
	o << "{" << v.x << " " << v.y << " " << v.z << "}";
	return o;
}

/**
 * Compute the distance between two vectors without using Periodic Boundary Conditions (PBC)
 * @param a First vector (position)
 * @param b Second vector (position)
 * @return The scalar distance between the two points
 */
inline Real distanceWithoutPBC(const Vector& a, const Vector& b) {
	Real dx= a.x - b.x;
	Real dy= a.y - b.y;
	Real dz= a.z - b.z;
	return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Compute the 2D distance between two vectors without using PBC
 * @param a First vector (position)
 * @param b Second vector (position)
 * @return The scalar distance between the two points, in 2D (ignoring z component)
 */
inline Real distance2D(const Vector& a, const Vector& b) {
	Real dx= a.x - b.x;
	Real dy= a.y - b.y;
	return std::sqrt(dx*dx + dy*dy);
}

/**
 * Compute the angle (in radians) between two vectors
 * @param a First vector
 * @param b Second vector
 * @return Angle between a and b in radians (0 <= angle <= PI)
 */
inline Real angleBetweenRadians(const Vector& a, const Vector& b) {
	Real dp= a*b, magA= a.magnitude(), magB= b.magnitude();
	if(magA < Vector::EPSILON || magB < Vector::EPSILON) return 0.0;

	Real cosTheta= dp / (magA * magB);

	// Clamp value to [-1, 1] to avoid NaNs due to floating point errors
	cosTheta= std::max(static_cast<Real>(-1.0), std::min(static_cast<Real>(1.0), cosTheta));

	return std::acos(cosTheta);
}

/**
 * Compute the angle (in degrees) between two vectors
 * @param a First vector
 * @param b Second vector
 * @return Angle between a and b in degrees (0 <= angle <= 180)
 */
inline Real angleBetweenDegrees(const Vector& a, const Vector& b) {
	return angleBetweenRadians(a,b) * Vector::RAD2DEG;
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




// ============================================================
//  Bounds box type selection
// ============================================================

#ifndef USE_TRUNCATED_OCTAHEDRON
    // BRANCH A - Orthorhombic (default)
    using BoundsType= Vector;
#else
    // BRANCH B - Truncated Octahedron / Triclinic
    /**
	 * BoxTO: periodic box for a truncated octahedron (or any triclinic cell).
	 *
	 * Constructed from CRYST1 data: lx, ly, lz, alpha, beta, gamma (degrees).
	 * The 3x3 cell matrix H and its inverse are computed once at construction.
	 *
	 * Layout of H (column = cell vector):
	 *       | ax  bx  cx |
	 *   H = | 0   by  cy |
	 *       | 0   0   cz |
	 * where ax=lx, bx=ly*cos(gamma), by=ly*sin(gamma), etc.
	 */
	struct BoxTO {
		Real lx, ly, lz; // cell lengths (Angstrom)
		Real h[3][3];    // cell matrix H
		Real hinv[3][3]; // H^{-1}, precomputed

		/**
		* Construct from CRYST1 parameters.
		* @param a  Cell length a (Angstrom)
		* @param b  Cell length b (Angstrom)
		* @param c  Cell length c (Angstrom)
		* @param alpha Angle alpha (degrees) — between b and c
		* @param beta  Angle beta  (degrees) — between a and c
		* @param gamma Angle gamma (degrees) — between a and b
		*/
		BoxTO(Real a, Real b, Real c, Real alpha_deg, Real beta_deg, Real gamma_deg): lx(a), ly(b), lz(c) {
			const Real alpha= alpha_deg / RAD2DEG;
			const Real beta = beta_deg  / RAD2DEG;
			const Real gamma= gamma_deg / RAD2DEG;

			const Real cos_a= std::cos(alpha);
			const Real cos_b= std::cos(beta);
			const Real cos_g= std::cos(gamma);
			const Real sin_g= std::sin(gamma);

			// Cell matrix H (upper-triangular convention, column = cell vector)
			// Column 0 (a vector): along x
			h[0][0]= a;
			h[1][0]= 0.0;
			h[2][0]= 0.0;

			// Column 1 (b vector): in xy plane
			h[0][1]= b * cos_g;
			h[1][1]= b * sin_g;
			h[2][1]= 0.0;

			// Column 2 (c vector): general
			h[0][2]= c * cos_b;
			h[1][2]= c * (cos_a - cos_b*cos_g) / sin_g;
			h[2][2]= std::sqrt(c*c - h[0][2]*h[0][2] - h[1][2]*h[1][2]);

			// Analytical inverse of upper-triangular 3x3
			// | a  b  c |^(-1)   =   | 1/a   -b/(a*e)   (b*f-c*e)/(a*e*i) |
			// | 0  e  f |            |  0      1/e         -f/(e*i)       |
			// | 0  0  i |            |  0       0              1/i        |
			const Real inv_a= 1.0 / h[0][0];
			const Real inv_e= 1.0 / h[1][1];
			const Real inv_i= 1.0 / h[2][2];

			hinv[0][0]=  inv_a;
			hinv[0][1]= -h[0][1] * inv_a * inv_e;
			hinv[0][2]=  (h[0][1]*h[1][2] - h[0][2]*h[1][1]) * inv_a * inv_e * inv_i;
			hinv[1][0]=  0.0;
			hinv[1][1]=  inv_e;
			hinv[1][2]= -h[1][2] * inv_e * inv_i;
			hinv[2][0]=  0.0;
			hinv[2][1]=  0.0;
			hinv[2][2]=  inv_i;
		}

		/**
		* Convenience constructor from an orthorhombic Vector (alpha=beta=gamma=90°).
		* Useful for code that receives a BoundsType but started as orthorhombic.
		*/
		explicit BoxTO(const Vector& v): BoxTO(v.x, v.y, v.z, 90.0, 90.0, 90.0) {}
	};
    using BoundsType= BoxTO;
#endif

// ============================================================
//  Generic functions (independent of the box type)
// ============================================================

/**
 * Compute the displacement between two positions under PBC
 * @param a First position vector
 * @param b Second position vector
 * @param box Box dimensions (BoundsType could be orthorhombic or triclinic)
 * @return Displacement vector
 */
inline Vector displacementPBC(const Vector& a, const Vector& b, const BoundsType& box);

/**
 * Compute the squared minimum image distance between two positions under PBC
 * @param a First position vector
 * @param b Second position vector
 * @param box Box dimensions (BoundsType could be orthorhombic or triclinic)
 * @return Squared minimum image distance
 */
inline Real squaredDistancePBC(const Vector& a, const Vector& b, const BoundsType& box) {
    const Vector d= displacementPBC(a, b, box);
    return d.x*d.x + d.y*d.y + d.z*d.z;
}

/**
 * Compute the minimum-image distance between two vectors using Periodic Boundary Conditions (PBC)
 * @param a First vector (position)
 * @param b Second vector (position)
 * @param box Box dimensions (BoundsType could be orthorhombic or triclinic)
 * @return The scalar distance between the two points using minimum image convention
 */
inline Real distancePBC(const Vector& a, const Vector& b, const BoundsType& box) {
    return std::sqrt(squaredDistancePBC(a, b, box));
}

/**
 * Calculates the angle that 3 points form in a 3D system
 * @param c1 One of the point in the edges
 * @param c2 The center point
 * @param c3 The other point in an edge
 * @param box Box dimensions (BoundsType could be orthorhombic or triclinic)
 * @return The angle in radians formed by the 3 points
 */
inline Real getAngle(Vector c1, Vector c2, Vector c3, BoundsType box) {
    Real a= distancePBC(c1, c3, box);
    Real b= distancePBC(c1, c2, box);
    Real c= distancePBC(c2, c3, box);
    Real cos_angle= (b*b + c*c - a*a) / (2*b*c);
    cos_angle= std::max(static_cast<Real>(-1.0), std::min(static_cast<Real>(1.0), cos_angle));
    return std::fabs(std::acos(cos_angle));
}

// ============================================================
//  Specialized functions for each BoundsType
// ============================================================

#ifndef USE_TRUNCATED_OCTAHEDRON
// --------------------- Orthorhombic (Vector) ---------------------
inline Vector displacementPBC(const Vector& a, const Vector& b, const BoundsType& box) {
    Real dx= a.x - b.x;
    Real dy= a.y - b.y;
    Real dz= a.z - b.z;

    dx-= box.x * std::floor(dx / box.x + 0.5f);
    dy-= box.y * std::floor(dy / box.y + 0.5f);
    dz-= box.z * std::floor(dz / box.z + 0.5f);
    return Vector(dx, dy, dz);
}

#else // USE_TRUNCATED_OCTAHEDRON

// --------------------- Triclinic / Truncated Octahedron ---------------------
inline Vector displacementPBC(const Vector& a, const Vector& b, const BoundsType& box) {
    Vector dr= a - b;

    // fractional coordinates
    Real sx= box.hinv[0][0]*dr.x + box.hinv[0][1]*dr.y + box.hinv[0][2]*dr.z;
    Real sy=                       box.hinv[1][1]*dr.y + box.hinv[1][2]*dr.z;
    Real sz=                                             box.hinv[2][2]*dr.z;

    sx-= std::floor(sx + 0.5);
    sy-= std::floor(sy + 0.5);
    sz-= std::floor(sz + 0.5);

    // back to Cartesian
    const Real rx= box.h[0][0]*sx + box.h[0][1]*sy + box.h[0][2]*sz;
    const Real ry=                  box.h[1][1]*sy + box.h[1][2]*sz;
    const Real rz=                                   box.h[2][2]*sz;

    return Vector(rx, ry, rz);
}

#endif // USE_TRUNCATED_OCTAHEDRON

#endif // VECTOR_HPP
