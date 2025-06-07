#ifndef GEOMETRICS_HPP
#define GEOMETRICS_HPP

/**
 * Version: April 2025
 * Author: Nicolás Loubet
 */

#include "Vector.hpp"
#include <vector>

namespace Geometrics {
	const float R_V4S_O_HL= 1; //Anstrongs between O and perfect vertices
	const float PERFECT_ANGLE_FOR_TETRAHEDRON= std::acos(-1./3.)/2; //Perfect angle for tetrahedron
	
	struct TetrahedronVertices {
		Vector H1;
		Vector H2;
		Vector L1;
		Vector L2;

		Vector* toArray() {
			Vector* output= new Vector[4];
			output[0]= H1;
			output[1]= H2;
			output[2]= L1;
			output[3]= L2;
			return output;
		}

		std::vector<Vector> toVector() {
			std::vector<Vector> output;
			output.push_back(H1);
			output.push_back(H2);
			output.push_back(L1);
			output.push_back(L2);
			return output;
		}
	};
	
	/**
	 * Checks if all the molecule is in the same image of the box, if not, it moves all the atoms to the same image
	 * @param O Position of the oxygen
	 * @param H1 Position of the first hydrogen
	 * @param H2 Position of the second hydrogen
	 * @param bounds Periodic box dimensions
	 * @return An array with the positions of the oxygen and the two hydrogens, in that order
	 */
	Vector* checkPBC(Vector O, Vector H1, Vector H2, Vector bounds) {
		Vector* output= new Vector[3];

		output[0]= Vector(O.x,O.y,O.z);
		Vector H_i[2]= {H1, H2};
		for(int i= 0; i < 2; i++)
			output[i+1]= O + displacementPBC(H_i[i],O,bounds);
		return output;
	}


	/**
	 * Calculates the determinant of a 3x3 matrix with the Sarrus Rule
	 * @param matrix A 3x3 float array
	 * @return The matrix's determinant
	 */
	float determinant_3x3(float matrix[3][3]) {
		return matrix[0][0]*matrix[1][1]*matrix[2][2] + matrix[0][1]*matrix[1][2]*matrix[2][0] + matrix[0][2]*matrix[1][0]*matrix[2][1]
			- matrix[0][2]*matrix[1][1]*matrix[2][0] - matrix[0][0]*matrix[1][2]*matrix[2][1] - matrix[0][1]*matrix[1][0]*matrix[2][2];
	}

	/**
	 * Calculates the x,y,z values (a Vector) obteined from the Cramer's Rule
	 * @param matrix A 3x3 float array of the base matrix
	 * @param m_independents A float array of 3 values corresponding to the independt values
	 * @return The Vector with the x, y and z values obtained. (The Vector must be removed with "delete()")
	 */
	Vector CramersRule(float matrix[3][3], float m_independents[3]) {
		float A= determinant_3x3(matrix);
		float Ax_matrix[3][3]= {
			{ m_independents[0] , matrix[0][1] , matrix[0][2] },
			{ m_independents[1] , matrix[1][1] , matrix[1][2] },
			{ m_independents[2] , matrix[2][1] , matrix[2][2] }
		};
		float Ax= determinant_3x3(Ax_matrix);
		float Ay_matrix[3][3]= {
			{ matrix[0][0] , m_independents[0] , matrix[0][2] },
			{ matrix[1][0] , m_independents[1] , matrix[1][2] },
			{ matrix[2][0] , m_independents[2] , matrix[2][2] }
		};
		float Ay= determinant_3x3(Ay_matrix);
		float Az_matrix[3][3]= {
			{ matrix[0][0] , matrix[0][1] , m_independents[0] },
			{ matrix[1][0] , matrix[1][1] , m_independents[1] },
			{ matrix[2][0] , matrix[2][1] , m_independents[2] }
		};
		float Az= determinant_3x3(Az_matrix);
		return Vector(Ax/A, Ay/A, Az/A);
	}

	/**
	 * Calculates the perfect tetrahedron given the position of the oxygen and two hydrogens
	 * @param O_real Position of the oxygen
	 * @param H1_real Position of the first hydrogen
	 * @param H2_real Position of the second hydrogen
	 * @param bounds Periodic box dimensions
	 */
	TetrahedronVertices getPerfectTetrahedron(Vector O_real, Vector H1_real, Vector H2_real, Vector bounds) {
		Vector* pbc_vectors= checkPBC(O_real, H1_real, H2_real, bounds);
		Vector O= pbc_vectors[0];
		Vector H1= pbc_vectors[1];
		Vector H2= pbc_vectors[2];
		delete[] pbc_vectors;

		const float phy= angleBetweenRadians(H1_real-O_real,H2_real-O_real) / 2; //The /2 is because I need the angle between OH and b

		//I step
		Vector OH1= H1-O;
		Vector OH2= H2-O;

		//II step
		Vector b= OH1.getNormalized() + OH2.getNormalized();

		//III step
		Vector nu_OH= (OH1 % OH2);
		nu_OH.normalize();

		//IV step
		Vector hydrogens[2]= {OH1,OH2};
		Vector h[2];
		for(int i= 0; i < 2; i++) {
			Vector OHi= hydrogens[i];
			float k1= R_V4S_O_HL*b.magnitude()*std::cos(PERFECT_ANGLE_FOR_TETRAHEDRON);
			float k2= R_V4S_O_HL*OHi.magnitude()*std::cos(PERFECT_ANGLE_FOR_TETRAHEDRON-phy);

			float A_matriz[3][3]= {
				{   b.x   ,   b.y   ,   b.z   },
				{  OHi.x  ,  OHi.y  ,  OHi.z  },
				{ nu_OH.x , nu_OH.y , nu_OH.z }
			};
			float m_scalars[3]= { k1, k2, 0. };

			Vector oh= CramersRule(A_matriz, m_scalars);
			h[i]= oh + O;
		}

		//V step
		Vector m_H= (h[0]+h[1])*0.5f;
		Vector mH_H= h[0] - m_H;
		float delta= mH_H.magnitude();

		//VI step
		Vector O_2= O * 2.;
		Vector m_L= O_2 - m_H;

		//VII step
		TetrahedronVertices output;
		output.H1= h[0];
		output.H2= h[1];
		output.L1= (nu_OH * delta) + m_L;
		output.L2= (nu_OH * (-delta)) + m_L;

		return output;
	}

	
	/**
	 * Struct representing a 3D plane in space.
	 * The plane is defined by a normal vector and an independent term.
	 * The equation of the plane can be represented as:
	 *   normal.x * x + normal.y * y + normal.z * z + independent_term = 0
	 * normal_vector: (A ; B ; C)
	 * independent_term: D = - A x_0 - B y_0 - C z_0 = - (normal_vector * point_in_plane)
	 * Plane equation: Ax + By + Cz + D = 0
	 */
	struct Plane {
		Vector normal;
		float independent_term;

		Plane(): normal(0.0f, 0.0f, 0.0f), independent_term(0.0f) {}
		Plane(const Vector& n, float d): normal(n), independent_term(d) {}
	};

	/**
	 * Struct representing a 3D line in space.
	 * The line is defined by a direction vector and a point on the line.
	 * The equation of the line can be represented as:
	 *   point + t * dir
	 * where t is a scalar parameter.
	 */
	struct Line {
		Vector dir;
		Vector point;

		Line(): dir(0.0f, 0.0f, 0.0f), point(0.0f, 0.0f, 0.0f) {}
		Line(const Vector& d, const Vector& p): dir(d), point(p) {}
	};

	/**
	 * Struct representing a 3D box in space.
	 * The box is defined by two points: min_pos and max_pos
	 */
	struct Box {
		Vector min_pos;
		Vector max_pos;

		Box(): min_pos(0.0f, 0.0f, 0.0f), max_pos(0.0f, 0.0f, 0.0f) {}
		Box(const Vector& min, const Vector& max): min_pos(min), max_pos(max) {}
		Box(const float min_x, const float min_y, const float min_z, const float max_x, const float max_y, const float max_z):
			min_pos(min_x, min_y, min_z), max_pos(max_x, max_y, max_z) {}
	};

	/**
	 * 
	 */
	struct SphereList {
		int size;
		Vector* centers;

		SphereList(const int N_MAX) {size= 0; centers= new Vector[N_MAX];}
		SphereList(int size, Vector* centers): size(size), centers(centers) {}
		~SphereList() {delete[] centers;}
	};

	/**
	 * Calculates the Plane parameters from 3 points that are contained
	 * @param a One of the three points contained
	 * @param b Other of the three points contained
	 * @param c The last point contained
	 * @param bounds The coordinate of the last point, so the components are the width, height and length
	 * @return A Plane object
	 */
	Plane getPlaneFromPoints(Vector a, Vector b, Vector c, Vector bounds) {
		Vector normal= (b - a) % (c - a);
		normal.normalize();
		return Plane(normal, - (normal * a));
	}

	/**
	 * Calculates the distance of a point to a plane
	 * @param plane The plane object
	 * @param P The point
	 * @return float The distance perpendicular to the Plane up to the point P
	 */
	float distanceToPlane(Plane plane, Vector P) {
		return std::fabs((plane.normal * P) + plane.independent_term);
	}

	/**
	 * Calculates the distance of a point to a line
	 * @param line The line object
	 * @param P The point
	 * @param bounds The coordinate of the last point, so the components are the width, height and length
	 * @return float The distance perpendicular to the line up to the point P
	 */
	float distanceToLine(Line line, Vector P, Vector bounds) {
		float dist_P_P0= distancePBC(P,line.point,bounds);
		return std::sin( getAngle(P, line.point, line.point+line.dir, bounds) )*dist_P_P0;
	}

	/**
	 * Checks if a point is inside a box
	 * @param pos The point to check
	 * @param box The box to check
	 * @return true if the point is inside the box, false otherwise
	 */
	bool isInBox(Vector pos, Box box) {
		return pos.x >= box.min_pos.x && pos.x <= box.max_pos.x &&
			   pos.y >= box.min_pos.y && pos.y <= box.max_pos.y &&
			   pos.z >= box.min_pos.z && pos.z <= box.max_pos.z;
	}

}

#endif
