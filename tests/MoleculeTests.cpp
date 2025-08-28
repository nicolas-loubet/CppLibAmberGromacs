#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <sstream>
#include <cmath>
#include <limits>

TEST_CASE("Particle - Constructor and getters/setters", "[Particle]") {
	SECTION("Should initialize with default values") {
		Particle p;
		REQUIRE(p.getPosition().x == Approx(0.0));
		REQUIRE(p.getPosition().y == Approx(0.0));
		REQUIRE(p.getPosition().z == Approx(0.0));
		REQUIRE(p.getID() == 0);
		REQUIRE(p.getMass() == Approx(0.0));
		REQUIRE(p.getCharge() == Approx(0.0));
	}

	SECTION("Should initialize with position, ID, mass and charge") {
		Vector pos(1.0, 2.0, 3.0);
		Particle p(pos, 42, 18.0, 1.0);
		REQUIRE(p.getPosition().x == Approx(1.0));
		REQUIRE(p.getPosition().y == Approx(2.0));
		REQUIRE(p.getPosition().z == Approx(3.0));
		REQUIRE(p.getID() == 42);
		REQUIRE(p.getMass() == Approx(18.0));
		REQUIRE(p.getCharge() == Approx(1.0));
	}

	SECTION("Should copy correctly") {
		Vector pos(1.0, 2.0, 3.0);
		Particle p1(pos, 42, 18.0, 1.0);
		Particle p2(p1);
		REQUIRE(p2.getPosition().x == Approx(1.0));
		REQUIRE(p2.getID() == 42);
		REQUIRE(p2.getMass() == Approx(18.0));
		REQUIRE(p2.getCharge() == Approx(1.0));
	}

	SECTION("Should assign correctly") {
		Vector pos(1.0, 2.0, 3.0);
		Particle p1(pos, 42, 18.0, 1.0);
		Particle p2;
		p2 = p1;
		REQUIRE(p2.getPosition().x == Approx(1.0));
		REQUIRE(p2.getID() == 42);
		REQUIRE(p2.getMass() == Approx(18.0));
		REQUIRE(p2.getCharge() == Approx(1.0));
	}

	SECTION("Should calculate distance with PBC") {
		Vector pos1(0.0, 0.0, 0.0), pos2(1.0, 1.0, 1.0);
		Vector box(10.0, 10.0, 10.0);
		Particle p1(pos1, 1), p2(pos2, 2);
		Real dist = p1.distanceTo(p2, box);
		REQUIRE(dist == Approx(std::sqrt(3.0))); // Mock distancePBC uses Euclidean distance
	}

	SECTION("Should initialize with extreme values") {
		Vector pos(std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), std::numeric_limits<float>::min());
		Particle p(pos, 9999, std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
		REQUIRE(p.getPosition().x == Approx(std::numeric_limits<float>::max()));
		REQUIRE(p.getPosition().y == Approx(-std::numeric_limits<float>::max()));
		REQUIRE(p.getPosition().z == Approx(std::numeric_limits<float>::min()));
		REQUIRE(p.getID() == 9999);
		REQUIRE(p.getMass() == Approx(std::numeric_limits<float>::max()));
		REQUIRE(p.getCharge() == Approx(std::numeric_limits<float>::max()));
	}

	SECTION("Should initialize with values close to zero") {
		Vector pos(Vector::EPSILON, -Vector::EPSILON, 0.0);
		Particle p(pos, 0, Vector::EPSILON, Vector::EPSILON);
		REQUIRE(p.getPosition().x == Approx(Vector::EPSILON));
		REQUIRE(p.getPosition().y == Approx(-Vector::EPSILON));
		REQUIRE(p.getPosition().z == Approx(0.0));
		REQUIRE(p.getID() == 0);
		REQUIRE(p.getMass() == Approx(Vector::EPSILON));
		REQUIRE(p.getCharge() == Approx(Vector::EPSILON));
	}

	SECTION("Should calculate distance with PBC in multiple dimensions") {
		Vector pos1(0.0, 0.0, 0.0), pos2(11.0, 12.0, 13.0);
		Vector box(10.0, 10.0, 10.0);
		Particle p1(pos1, 1), p2(pos2, 2);
		Real dist = p1.distanceTo(p2, box);
		REQUIRE(dist == Approx(std::sqrt(1+4+9))); // 1.0, 2.0, 3.0 after PBC
	}

	SECTION("Should handle self-assignment") {
		Vector pos(1.0, 2.0, 3.0);
		Particle p(pos, 42, 18.0, 1.0);
		p = p;
		REQUIRE(p.getPosition().x == Approx(1.0));
		REQUIRE(p.getID() == 42);
		REQUIRE(p.getMass() == Approx(18.0));
		REQUIRE(p.getCharge() == Approx(1.0));
	}
}

TEST_CASE("Atom - Constructor and getters/setters", "[Atom]") {
	SECTION("Should initialize with default values") {
		Atom a;
		REQUIRE(a.getPosition().x == Approx(0.0));
		REQUIRE(a.getID() == 0);
		REQUIRE(a.getMass() == Approx(0.0));
		REQUIRE(a.getCharge() == Approx(0.0));
		REQUIRE(a.getis_Hatom() == false);
		REQUIRE(a.getSigma() == Approx(0.0));
		REQUIRE(a.getEpsilon() == Approx(0.0));
		REQUIRE(a.getZ() == 0);
	}

	SECTION("Should initialize with full parameters") {
		Vector pos(1.0, 2.0, 3.0);
		Atom a(pos, 42, 16.0, -0.5, 0.1, 3.5, 8);
		REQUIRE(a.getPosition().x == Approx(1.0));
		REQUIRE(a.getID() == 42);
		REQUIRE(a.getMass() == Approx(16.0));
		REQUIRE(a.getCharge() == Approx(-0.5));
		REQUIRE(a.getis_Hatom() == false); // Z=8 (Oxygen), not 1
		REQUIRE(a.getSigma() == Approx(3.5));
		REQUIRE(a.getEpsilon() == Approx(0.1));
		REQUIRE(a.getZ() == 8);
	}

	SECTION("Should set is_HAtom correctly with Z=1") {
		Atom a;
		a.setZ(1);
		REQUIRE(a.getis_Hatom() == true);
		a.setZ(6);
		REQUIRE(a.getis_Hatom() == false);
	}

	SECTION("Should copy correctly") {
		Vector pos(1.0, 2.0, 3.0);
		Atom a1(pos, 42, 16.0, -0.5, 0.1, 3.5, 8);
		Atom a2(a1);
		REQUIRE(a2.getPosition().x == Approx(1.0));
		REQUIRE(a2.getID() == 42);
		REQUIRE(a2.getMass() == Approx(16.0));
		REQUIRE(a2.getCharge() == Approx(-0.5));
		REQUIRE(a2.getis_Hatom() == false);
		REQUIRE(a2.getSigma() == Approx(3.5));
		REQUIRE(a2.getEpsilon() == Approx(0.1));
		REQUIRE(a2.getZ() == 8);
	}

	SECTION("Should assign correctly") {
		Vector pos(1.0, 2.0, 3.0);
		Atom a1(pos, 42, 16.0, -0.5, 0.1, 3.5, 8);
		Atom a2;
		a2 = a1;
		REQUIRE(a2.getPosition().x == Approx(1.0));
		REQUIRE(a2.getID() == 42);
		REQUIRE(a2.getMass() == Approx(16.0));
		REQUIRE(a2.getCharge() == Approx(-0.5));
		REQUIRE(a2.getis_Hatom() == false);
		REQUIRE(a2.getSigma() == Approx(3.5));
		REQUIRE(a2.getEpsilon() == Approx(0.1));
		REQUIRE(a2.getZ() == 8);
	}

	SECTION("Should initialize with negative values") {
		Vector pos(-1.0, -2.0, -3.0);
		Atom a(pos, -1, -16.0, -1.0, -0.1, -3.5, 8);
		REQUIRE(a.getPosition().x == Approx(-1.0));
		REQUIRE(a.getID() == -1);
		REQUIRE(a.getMass() == Approx(-16.0));
		REQUIRE(a.getCharge() == Approx(-1.0));
		REQUIRE(a.getis_Hatom() == false);
		REQUIRE(a.getSigma() == Approx(-3.5));
		REQUIRE(a.getEpsilon() == Approx(-0.1));
		REQUIRE(a.getZ() == 8);
	}

	SECTION("Should initialize with Z=0") {
		Atom a(Vector(0.0, 0.0, 0.0), 0, 0.0, 0.0, 0.0, 0.0, 0);
		REQUIRE(a.getZ() == 0);
		REQUIRE(a.getis_Hatom() == false);
	}

	SECTION("Should handle extreme values for sigma and epsilon") {
		Atom a(Vector(0.0, 0.0, 0.0), 1, 1.0, 1.0, std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 1);
		REQUIRE(a.getSigma() == Approx(std::numeric_limits<float>::max()));
		REQUIRE(a.getEpsilon() == Approx(std::numeric_limits<float>::max()));
		REQUIRE(a.getis_Hatom() == true);
	}

	SECTION("Should handle self-assignment") {
		Vector pos(1.0, 2.0, 3.0);
		Atom a(pos, 42, 16.0, -0.5, 0.1, 3.5, 8);
		a = a;
		REQUIRE(a.getPosition().x == Approx(1.0));
		REQUIRE(a.getID() == 42);
		REQUIRE(a.getMass() == Approx(16.0));
		REQUIRE(a.getCharge() == Approx(-0.5));
		REQUIRE(a.getis_Hatom() == false);
		REQUIRE(a.getSigma() == Approx(3.5));
		REQUIRE(a.getEpsilon() == Approx(0.1));
		REQUIRE(a.getZ() == 8);
	}

	SECTION("Should set position correctly") {
		Atom a;
		Vector pos(4.0, 5.0, 6.0);
		a.setPosition(pos);
		REQUIRE(a.getPosition() == pos);
	}
}

TEST_CASE("Atom - Coulomb Potential", "[Atom]") {
	Vector box(10.0, 10.0, 10.0);
	Atom a1(Vector(0.0, 0.0, 0.0), 1, 16.0, 1.0, 0.1, 3.5, 8);
	Atom a2(Vector(1.0, 0.0, 0.0), 2, 1.0, -1.0, 0.2, 2.5, 1);

	SECTION("Should compute Coulomb potential correctly") {
		Real potential = a1.getCoulombPotential(a2, box);
		Real expected = (1389.35458 * 1.0 * -1.0) / 1.0; // K_COULOMB * q1 * q2 / distance
		REQUIRE(potential == Approx(expected));
	}

	SECTION("Should handle zero charge") {
		a1.setCharge(0.0);
		Real potential = a1.getCoulombPotential(a2, box);
		REQUIRE(potential == Approx(0.0));
	}

	SECTION("Should handle zero distance with PBC") {
		Atom a2_same(Vector(10.0, 0.0, 0.0), 2, 1.0, -1.0, 0.2, 2.5, 1);
		Real potential = a1.getCoulombPotential(a2_same, box);
		REQUIRE((std::isinf(potential) || std::isnan(potential))); // Distance 0 after PBC
	}

	SECTION("Should handle extreme charges") {
		a1.setCharge(std::numeric_limits<float>::max());
		a2.setCharge(-std::numeric_limits<float>::max());
		Real potential = a1.getCoulombPotential(a2, box);
		REQUIRE((std::isinf(potential) || std::isnan(potential)));
	}

	SECTION("Should handle box with negative dimensions") {
		Vector box_neg(-10.0, -10.0, -10.0);
		Real potential = a1.getCoulombPotential(a2, box_neg);
		REQUIRE(potential == Approx((1389.35458 * 1.0 * -1.0) / 1.0)); // Falls back to non-PBC
	}
}

TEST_CASE("Atom - PDB Format", "[Atom]") {
	Atom a(Vector(1.234, 2.345, 3.456), 42, 16.0, -0.5, 0.1, 3.5, 8);

	SECTION("Should generate PDB format with default values") {
		std::string pdb = a.toPDBFormat();
		REQUIRE(pdb.find("ATOM      0 UNK  UNK     0    ") != std::string::npos);
		REQUIRE(pdb.find("   1.234   2.345   3.456  1.00  0.00           O ") != std::string::npos);
	}

	SECTION("Should generate PDB format with custom parameters") {
		std::string pdb = a.toPDBFormat(1, 100, "O1", "WAT", "A", 0.95, 1.23);
		REQUIRE(pdb.find("ATOM      1 O1   WAT A 100    ") != std::string::npos);
		REQUIRE(pdb.find("   1.234   2.345   3.456  0.95  1.23           O ") != std::string::npos);
	}

	SECTION("Should handle negative coordinates in PDB format") {
		Atom a_neg(Vector(-1.234, -2.345, -3.456), 42, 16.0, -0.5, 0.1, 3.5, 8);
		std::string pdb = a_neg.toPDBFormat();
		REQUIRE(pdb.find("  -1.234  -2.345  -3.456  1.00  0.00           O ") != std::string::npos);
	}

	SECTION("Should handle Z=1 in PDB format") {
		Atom a_h(Vector(1.234, 2.345, 3.456), 42, 1.0, 0.417, 0.05, 2.5, 1);
		std::string pdb = a_h.toPDBFormat();
		REQUIRE(pdb.find("           H ") != std::string::npos);
	}

	SECTION("Should handle extreme values in PDB format") {
		Atom a_ext(Vector(std::numeric_limits<float>::max(), 0.0, 0.0), 42, 16.0, -0.5, 0.1, 3.5, 8);
		std::string pdb = a_ext.toPDBFormat();
		REQUIRE(pdb.find(" O ") != std::string::npos); // Ensure element symbol persists
	}
}

TEST_CASE("Molecule - Constructor and getters/setters", "[Molecule]") {
	SECTION("Should initialize with default values") {
		Molecule m;
		REQUIRE(m.getNAtoms() == 0);
		REQUIRE(m.getAtoms() == nullptr);
		REQUIRE(m.isWater() == false);
		REQUIRE(m.getClassification() == NOT_CLASSIFIED);
		REQUIRE(m.isClassified() == false);
	}

	SECTION("Should initialize with atoms") {
		Atom* atoms = new Atom[2];
		atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.417*2, 0.1, 3.5, 8);
		atoms[1] = Atom(Vector(1.0, 0.0, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
		Molecule m(42, atoms, 2);
		REQUIRE(m.getNAtoms() == 2);
		REQUIRE(m.getID() == 42);
		REQUIRE(m.getMass() == Approx(17.0));
		REQUIRE(m.getCharge() == Approx(-0.417));
		REQUIRE(m.getPosition().x == Approx(0.5)); // Center of mass
		REQUIRE(m.isWater() == false);
	}

	SECTION("Should copy correctly") {
		Atom* atoms = new Atom[1];
		atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
		Molecule m1(42, atoms, 1);
		Molecule m2(m1);
		REQUIRE(m2.getNAtoms() == 1);
		REQUIRE(m2.getID() == 42);
		REQUIRE(m2.getMass() == Approx(16.0));
		REQUIRE(m2.getAtoms()[0].getZ() == 8);
	}

	SECTION("Should assign correctly") {
		Atom* atoms = new Atom[1];
		atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
		Molecule m1(42, atoms, 1);
		Molecule m2;
		m2 = m1;
		REQUIRE(m2.getNAtoms() == 1);
		REQUIRE(m2.getID() == 42);
		REQUIRE(m2.getMass() == Approx(16.0));
		REQUIRE(m2.getAtoms()[0].getZ() == 8);
	}

	SECTION("Should handle classification") {
		Molecule m;
		m.setClassification(100);
		REQUIRE(m.getClassification() == 100);
		REQUIRE(m.isClassified() == true);
		m.removeClassification();
		REQUIRE(m.getClassification() == NOT_CLASSIFIED);
		REQUIRE(m.isClassified() == false);
	}

	SECTION("Should handle single atom") {
		Atom* atoms = new Atom[1];
		atoms[0] = Atom(Vector(1.0, 2.0, 3.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
		Molecule m(42, atoms, 1);
		REQUIRE(m.getNAtoms() == 1);
		REQUIRE(m.getMass() == Approx(16.0));
		REQUIRE(m.getCharge() == Approx(-0.834));
		REQUIRE(m.getPosition() == Vector(1.0, 2.0, 3.0));
	}

	SECTION("Should handle self-assignment") {
		Atom* atoms = new Atom[1];
		atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
		Molecule m(42, atoms, 1);
		m = m;
		REQUIRE(m.getNAtoms() == 1);
		REQUIRE(m.getID() == 42);
		REQUIRE(m.getMass() == Approx(16.0));
		REQUIRE(m.getAtoms()[0].getZ() == 8);
	}
}

TEST_CASE("Molecule - Center of Mass", "[Molecule]") {
	SECTION("Should compute center of mass correctly") {
		Atom atoms[] = {Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, 0.0, 0.1, 3.5, 8),
						Atom(Vector(2.0, 4.0, 6.0), 2, 1.0, 0.0, 0.05, 2.5, 1)};
		Vector com = Molecule::centerOfMass(atoms, 2);
		REQUIRE(com.x == Approx(1.0));
		REQUIRE(com.y == Approx(2.0));
		REQUIRE(com.z == Approx(3.0));
	}

	SECTION("Should compute center of mass with a single atom") {
		Atom atoms[] = {Atom(Vector(1.0, 2.0, 3.0), 1, 16.0, 0.0, 0.1, 3.5, 8)};
		Vector com = Molecule::centerOfMass(atoms, 1);
		REQUIRE(com.x == Approx(1.0));
		REQUIRE(com.y == Approx(2.0));
		REQUIRE(com.z == Approx(3.0));
	}

	SECTION("Should handle center of mass with zero mass atom") {
		Atom atoms[] = {Atom(Vector(1.0, 2.0, 3.0), 1, 0.0, 0.0, 0.1, 3.5, 8)};
		Vector com = Molecule::centerOfMass(atoms, 1);
		REQUIRE(com.x == Approx(1.0));
		REQUIRE(com.y == Approx(2.0));
		REQUIRE(com.z == Approx(3.0));
	}
}

TEST_CASE("Molecule - Nearby atoms", "[Molecule]") {
	Atom* atoms = new Atom[3];
	atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
	atoms[1] = Atom(Vector(1.0, 0.0, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
	atoms[2] = Atom(Vector(10.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
	Molecule m(42, atoms, 3);
	Vector box(100.0, 100.0, 100.0);

	SECTION("Should find nearby atoms") {
		auto nearby = m.findNearbyAtoms(1, 2.0, box);
		REQUIRE(nearby.size() == 1);
		REQUIRE(nearby[0] == 2);
	}

	SECTION("Should return empty vector if no nearby atoms") {
		auto nearby = m.findNearbyAtoms(1, 0.5, box);
		REQUIRE(nearby.empty());
	}

	SECTION("Should find nearby atoms with PBC") {
		Vector box_small(5.0, 5.0, 5.0);
		auto nearby = m.findNearbyAtoms(1, 12.0, box_small);
		REQUIRE(nearby.size() == 2); // Atom 2 (1.0 away), Atom 3 (0.0 away after PBC)
		REQUIRE(nearby[0] == 2);
		REQUIRE(nearby[1] == 3);
	}

	SECTION("Should handle zero distance") {
		auto nearby = m.findNearbyAtoms(1, 0.0, box);
		REQUIRE(nearby.empty());
	}
}
TEST_CASE("Water - Constructor and getters", "[Water]") {
	Atom* atoms = new Atom[3];
	atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
	atoms[1] = Atom(Vector(1.0, 0.0, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
	atoms[2] = Atom(Vector(10.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
	Water w(42, atoms, 3);

	SECTION("Should initialize correctly") {
		REQUIRE(w.getNAtoms() == 3);
		REQUIRE(w.getID() == 42);
		REQUIRE(w.isWater() == true);
		REQUIRE(w.getOxygen().getZ() == 8);
		REQUIRE(w.getHydrogen_1().getZ() == 1);
		REQUIRE(w.getHydrogen_2().getZ() == 1);
	}

	SECTION("Should copy correctly") {
		Water w2(w);
		REQUIRE(w2.getNAtoms() == 3);
		REQUIRE(w2.getID() == 42);
		REQUIRE(w2.isWater() == true);
		REQUIRE(w2.getOxygen().getZ() == 8);
	}

	SECTION("Should assign correctly") {
		Water w2;
		w2 = w;
		REQUIRE(w2.getNAtoms() == 3);
		REQUIRE(w2.getID() == 42);
		REQUIRE(w2.isWater() == true);
		REQUIRE(w2.getOxygen().getZ() == 8);
	}
}

TEST_CASE("Water - Hydrogen bond", "[Water]") {
	Atom* atoms1 = new Atom[3];
	atoms1[0] = Atom(Vector(38.69, 38.95, 20.04), 1, 16.0, -0.834, 0.1, 3.5, 8);
	atoms1[1] = Atom(Vector(39.54, 39.38, 20.05), 2, 1.0, 0.417, 0.05, 2.5, 1);
	atoms1[2] = Atom(Vector(38.50, 38.82, 19.11), 3, 1.0, 0.417, 0.05, 2.5, 1);

	Atom* atoms2 = new Atom[3];
	atoms2[0] = Atom(Vector(41.01, 40.49, 19.97), 4, 16.0, -0.834, 0.1, 3.5, 8);
	atoms2[1] = Atom(Vector(41.29, 40.68, 19.08), 5, 1.0, 0.417, 0.05, 2.5, 1);
	atoms2[2] = Atom(Vector(41.83, 40.43, 20.47), 6, 1.0, 0.417, 0.05, 2.5, 1);

	Water w1(1, atoms1, 3);
	Water w2(2, atoms2, 3);
	Vector box(100.0, 100.0, 100.0);

	SECTION("Should detect hydrogen bond") {
		REQUIRE(w1.isHB(w2, box, 3.5, 30.0) == true); // Mock getAngle returns 0.1 rad (~5.73 deg)
	}

	SECTION("Should not detect hydrogen bond if distance is too large") {
		Atom* atoms3 = new Atom[3];
		atoms3[0] = Atom(Vector(48.01, 40.49, 19.97), 4, 16.0, -0.834, 0.1, 3.5, 8);
		atoms3[1] = Atom(Vector(48.29, 40.68, 19.08), 5, 1.0, 0.417, 0.05, 2.5, 1);
		atoms3[2] = Atom(Vector(48.83, 40.43, 20.47), 6, 1.0, 0.417, 0.05, 2.5, 1);
		Water w3(3, atoms3, 3);
		REQUIRE(w1.isHB(w3, box, 3.5, 30.0) == false); // Distance > 3.5
	}

	SECTION("Should not detect hydrogen bond if angle is too large") {
		Atom* atoms4 = new Atom[3];
		atoms4[0] = Atom(Vector(41.01, 40.49, 19.97), 4, 16.0, -0.834, 0.1, 3.5, 8);
		atoms4[1] = Atom(Vector(41.29, 40.68, 19.08), 5, 1.0, 0.417, 0.05, 2.5, 1);
		atoms4[2] = Atom(Vector(41.83, 40.43, 20.47), 6, 1.0, 0.417, 0.05, 2.5, 1);
		Water w4(4, atoms4, 3);
		REQUIRE(w1.isHB(w4, box, 3.5, 0.1) == false); // Angle > 0.1 rad
	}

	SECTION("Should handle zero distance in hydrogen bond") {
		Atom* atoms5 = new Atom[3];
		atoms5[0] = Atom(Vector(38.69, 38.95, 20.04), 4, 16.0, -0.834, 0.1, 3.5, 8);
		atoms5[1] = Atom(Vector(39.54, 39.38, 20.05), 5, 1.0, 0.417, 0.05, 2.5, 1);
		atoms5[2] = Atom(Vector(38.50, 38.82, 19.11), 6, 1.0, 0.417, 0.05, 2.5, 1);
		Water w5(5, atoms5, 3);
		REQUIRE(w1.isHB(w5, box, 3.5, 30.0) == false); // Same position
	}
}

TEST_CASE("Water - Lennard-Jones and total potential", "[Water]") {
	Atom* atoms = new Atom[2];
	atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
	atoms[1] = Atom(Vector(0.757, 0.586, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
	Water w(1, atoms, 2);
	Vector box(10.0, 10.0, 10.0);

	SECTION("Should compute Lennard-Jones potential") {
		Atom a(Vector(1.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
		Real s = (3.5 + 2.5) / 2.0; // Lorentz-Berthelot
		Real e = std::sqrt(0.1 * 0.05);
		Real s_over_R = s / 1.0; // Distance = 1.0
		Real expected = 4 * e * (std::pow(s_over_R, 12) - std::pow(s_over_R, 6));
		REQUIRE(w.getLJPotential(a, s, e, box) == Approx(expected));
	}

	SECTION("Should compute total potential with another water molecule") {
		Atom* atoms2 = new Atom[2];
		atoms2[0] = Atom(Vector(1.0, 0.0, 0.0), 3, 16.0, -0.834, 0.1, 3.5, 8);
		atoms2[1] = Atom(Vector(1.757, 0.586, 0.0), 4, 1.0, 0.417, 0.05, 2.5, 1);
		Water w2(2, atoms2, 2);
		Real s = 3.5; // Oxygen-Oxygen
		Real e = 0.1;
		Real s_over_R = s / 1.0;
		Real lj = 4 * e * (std::pow(s_over_R, 12) - std::pow(s_over_R, 6));
		Real coulomb = 0.0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				coulomb += atoms[i].getCoulombPotential(atoms2[j], box);
		REQUIRE(w.potentialWith(w2, box) == Approx(lj + coulomb));
	}

	SECTION("Should compute total potential with an atom") {
		Atom a(Vector(1.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
		Real s = (3.5 + 2.5) / 2.0;
		Real e = std::sqrt(0.1 * 0.05);
		Real s_over_R = s / 1.0;
		Real lj = 4 * e * (std::pow(s_over_R, 12) - std::pow(s_over_R, 6));
		Real coulomb = 0.0;
		for (int i = 0; i < 2; i++)
			coulomb += atoms[i].getCoulombPotential(a, box);
		REQUIRE(w.potentialWith(a, box) == Approx(lj + coulomb));
	}

	SECTION("Should handle Lennard-Jones potential with zero sigma") {
		Atom a(Vector(1.0, 0.0, 0.0), 3, 1.0, 0.417, 0.0, 0.0, 1);
		Real s = (3.5 + 0.0) / 2.0;
		Real e = std::sqrt(0.1 * 0.0);
		Real lj = w.getLJPotential(a, s, e, box);
		REQUIRE(lj == Approx(0.0));
	}

	SECTION("Should handle Lennard-Jones potential with zero distance") {
		Atom a(Vector(0.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
		Real s = (3.5 + 2.5) / 2.0;
		Real e = std::sqrt(0.1 * 0.05);
		Real lj = w.getLJPotential(a, s, e, box);
		REQUIRE((std::isinf(lj) || std::isnan(lj)));
	}

	SECTION("Should handle total potential with identical molecule") {
		Water w2(w);
		Real lj = 4 * 0.1 * (std::pow(3.5 / 0.0, 12) - std::pow(3.5 / 0.0, 6));
		Real coulomb = 0.0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				coulomb += atoms[i].getCoulombPotential(w2.getAtoms()[j], box);
		REQUIRE((std::isinf(w.potentialWith(w2, box)) || std::isnan(w.potentialWith(w2, box))));
	}

	SECTION("Should handle total potential with PBC") {
		Atom* atoms2 = new Atom[2];
		atoms2[0] = Atom(Vector(11.0, 0.0, 0.0), 3, 16.0, -0.834, 0.1, 3.5, 8);
		atoms2[1] = Atom(Vector(11.757, 0.586, 0.0), 4, 1.0, 0.417, 0.05, 2.5, 1);
		Water w2(2, atoms2, 2);
		Real s = 3.5;
		Real e = 0.1;
		Real s_over_R = s / 1.0; // PBC distance = 1.0
		Real lj = 4 * e * (std::pow(s_over_R, 12) - std::pow(s_over_R, 6));
		Real coulomb = 0.0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				coulomb += atoms[i].getCoulombPotential(atoms2[j], box);
		REQUIRE(w.potentialWith(w2, box) == Approx(lj + coulomb));
	}

}
