#include "catch.hpp"
#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <sstream>
#include <cmath>

TEST_CASE("Particle - Constructor y getters/setters", "[Particle]") {
    SECTION("Debería inicializar con valores por defecto") {
        Particle p;
        REQUIRE(p.getPosition().x == Approx(0.0));
        REQUIRE(p.getPosition().y == Approx(0.0));
        REQUIRE(p.getPosition().z == Approx(0.0));
        REQUIRE(p.getID() == 0);
        REQUIRE(p.getMass() == Approx(0.0));
        REQUIRE(p.getCharge() == Approx(0.0));
    }

    SECTION("Debería inicializar con posición, ID, masa y carga") {
        Vector pos(1.0, 2.0, 3.0);
        Particle p(pos, 42, 18.0, 1.0);
        REQUIRE(p.getPosition().x == Approx(1.0));
        REQUIRE(p.getPosition().y == Approx(2.0));
        REQUIRE(p.getPosition().z == Approx(3.0));
        REQUIRE(p.getID() == 42);
        REQUIRE(p.getMass() == Approx(18.0));
        REQUIRE(p.getCharge() == Approx(1.0));
    }

    SECTION("Debería copiar correctamente") {
        Vector pos(1.0, 2.0, 3.0);
        Particle p1(pos, 42, 18.0, 1.0);
        Particle p2(p1);
        REQUIRE(p2.getPosition().x == Approx(1.0));
        REQUIRE(p2.getID() == 42);
        REQUIRE(p2.getMass() == Approx(18.0));
        REQUIRE(p2.getCharge() == Approx(1.0));
    }

    SECTION("Debería asignar correctamente") {
        Vector pos(1.0, 2.0, 3.0);
        Particle p1(pos, 42, 18.0, 1.0);
        Particle p2;
        p2 = p1;
        REQUIRE(p2.getPosition().x == Approx(1.0));
        REQUIRE(p2.getID() == 42);
        REQUIRE(p2.getMass() == Approx(18.0));
        REQUIRE(p2.getCharge() == Approx(1.0));
    }

    SECTION("Debería calcular distancia con PBC") {
        Vector pos1(0.0, 0.0, 0.0), pos2(1.0, 1.0, 1.0);
        Vector box(10.0, 10.0, 10.0);
        Particle p1(pos1, 1), p2(pos2, 2);
        Real dist = p1.distanceTo(p2, box);
        REQUIRE(dist == Approx(std::sqrt(3.0))); // Mock distancePBC uses Euclidean distance
    }
}

TEST_CASE("Atom - Constructor y getters/setters", "[Atom]") {
    SECTION("Debería inicializar con valores por defecto") {
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

    SECTION("Debería inicializar con parámetros completos") {
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

    SECTION("Debería establecer is_HAtom correctamente con Z=1") {
        Atom a;
        a.setZ(1);
        REQUIRE(a.getis_Hatom() == true);
        a.setZ(6);
        REQUIRE(a.getis_Hatom() == false);
    }

    SECTION("Debería copiar correctamente") {
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

    SECTION("Debería asignar correctamente") {
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
}

TEST_CASE("Atom - Coulomb Potential", "[Atom]") {
    Vector box(10.0, 10.0, 10.0);
    Atom a1(Vector(0.0, 0.0, 0.0), 1, 16.0, 1.0, 0.1, 3.5, 8);
    Atom a2(Vector(1.0, 0.0, 0.0), 2, 1.0, -1.0, 0.2, 2.5, 1);

    SECTION("Debería calcular el potencial de Coulomb correctamente") {
        Real potential = a1.getCoulombPotential(a2, box);
        Real expected = (1389.35458 * 1.0 * -1.0) / 1.0; // K_COULOMB * q1 * q2 / distance
        REQUIRE(potential == Approx(expected));
    }

    SECTION("Debería manejar carga cero") {
        a1.setCharge(0.0);
        Real potential = a1.getCoulombPotential(a2, box);
        REQUIRE(potential == Approx(0.0));
    }
}

TEST_CASE("Atom - Formato PDB", "[Atom]") {
    Atom a(Vector(1.234, 2.345, 3.456), 42, 16.0, -0.5, 0.1, 3.5, 8);

    SECTION("Debería generar formato PDB con valores por defecto") {
        std::string pdb = a.toPDBFormat();
        REQUIRE(pdb.find("ATOM      0 UNK  UNK     0    ") != std::string::npos);
        REQUIRE(pdb.find("   1.234   2.345   3.456  1.00  0.00           O ") != std::string::npos);
    }

    SECTION("Debería generar formato PDB con parámetros personalizados") {
        std::string pdb = a.toPDBFormat(1, 100, "O1", "WAT", "A", 0.95, 1.23);
        REQUIRE(pdb.find("ATOM      1 O1   WAT A 100    ") != std::string::npos);
        REQUIRE(pdb.find("   1.234   2.345   3.456  0.95  1.23           O ") != std::string::npos);
    }
}

TEST_CASE("Molecule - Constructor y getters/setters", "[Molecule]") {
    SECTION("Debería inicializar con valores por defecto") {
        Molecule m;
        REQUIRE(m.getNAtoms() == 0);
        REQUIRE(m.getAtoms() == nullptr);
        REQUIRE(m.isWater() == false);
        REQUIRE(m.getClassification() == NOT_CLASSIFIED);
        REQUIRE(m.isClassified() == false);
    }

    SECTION("Debería inicializar con átomos") {
        Atom* atoms= new Atom[2];
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

    SECTION("Debería copiar correctamente") {
        Atom* atoms= new Atom[1];
        atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
        Molecule m1(42, atoms, 1);
        Molecule m2(m1);
        REQUIRE(m2.getNAtoms() == 1);
        REQUIRE(m2.getID() == 42);
        REQUIRE(m2.getMass() == Approx(16.0));
        REQUIRE(m2.getAtoms()[0].getZ() == 8);
    }

    SECTION("Debería asignar correctamente") {
        Atom* atoms= new Atom[1];
        atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
        Molecule m1(42, atoms, 1);
        Molecule m2;
        m2 = m1;
        REQUIRE(m2.getNAtoms() == 1);
        REQUIRE(m2.getID() == 42);
        REQUIRE(m2.getMass() == Approx(16.0));
        REQUIRE(m2.getAtoms()[0].getZ() == 8);
    }

    SECTION("Debería manejar clasificación") {
        Molecule m;
        m.setClassification(100);
        REQUIRE(m.getClassification() == 100);
        REQUIRE(m.isClassified() == true);
        m.removeClassification();
        REQUIRE(m.getClassification() == NOT_CLASSIFIED);
        REQUIRE(m.isClassified() == false);
    }
}

TEST_CASE("Molecule - Centro de masa", "[Molecule]") {
    SECTION("Debería calcular el centro de masa correctamente") {
        Atom atoms[] = {Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, 0.0, 0.1, 3.5, 8),
                        Atom(Vector(2.0, 4.0, 6.0), 2, 1.0, 0.0, 0.05, 2.5, 1)};
        Vector com = Molecule::centerOfMass(atoms, 2);
        REQUIRE(com.x == Approx(1.0));
        REQUIRE(com.y == Approx(2.0));
        REQUIRE(com.z == Approx(3.0));
    }
}

TEST_CASE("Molecule - Átomos cercanos", "[Molecule]") {
    Atom* atoms= new Atom[3];
    atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
    atoms[1] = Atom(Vector(1.0, 0.0, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
    atoms[2] = Atom(Vector(10.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
    Molecule m(42, atoms, 3);
    Vector box(100.0, 100.0, 100.0);

    SECTION("Debería encontrar átomos cercanos") {
        auto nearby = m.findNearbyAtoms(1, 2.0, box);
        REQUIRE(nearby.size() == 1);
        REQUIRE(nearby[0] == 2);
    }

    SECTION("Debería devolver vector vacío si no hay átomos cercanos") {
        auto nearby = m.findNearbyAtoms(1, 0.5, box);
        REQUIRE(nearby.empty());
    }
}

TEST_CASE("Water - Constructor y getters", "[Water]") {
    Atom* atoms= new Atom[3];
    atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
    atoms[1] = Atom(Vector(1.0, 0.0, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
    atoms[2] = Atom(Vector(10.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
    Water w(42, atoms, 3);

    SECTION("Debería inicializar correctamente") {
        REQUIRE(w.getNAtoms() == 3);
        REQUIRE(w.getID() == 42);
        REQUIRE(w.isWater() == true);
        REQUIRE(w.getOxygen().getZ() == 8);
        REQUIRE(w.getHydrogen_1().getZ() == 1);
        REQUIRE(w.getHydrogen_2().getZ() == 1);
    }

    SECTION("Debería copiar correctamente") {
        Water w2(w);
        REQUIRE(w2.getNAtoms() == 3);
        REQUIRE(w2.getID() == 42);
        REQUIRE(w2.isWater() == true);
        REQUIRE(w2.getOxygen().getZ() == 8);
    }
}

TEST_CASE("Water - Enlace de hidrógeno", "[Water]") {
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

    SECTION("Debería detectar enlace de hidrógeno") {
        REQUIRE(w1.isHB(w2, box, 3.5, 30.0) == true); // Mock getAngle returns 0.1 rad (~5.73 deg)
    }

    SECTION("Debería no detectar enlace de hidrógeno si la distancia es demasiado grande") {
        Atom* atoms3 = new Atom[3];
        atoms3[0] = Atom(Vector(48.01, 40.49, 19.97), 4, 16.0, -0.834, 0.1, 3.5, 8);
        atoms3[1] = Atom(Vector(48.29, 40.68, 19.08), 5, 1.0, 0.417, 0.05, 2.5, 1);
        atoms3[2] = Atom(Vector(48.83, 40.43, 20.47), 6, 1.0, 0.417, 0.05, 2.5, 1);
        Water w3(3, atoms3, 3);
        REQUIRE(w1.isHB(w3, box, 3.5, 30.0) == false); // Distance > 3.5
    }
}

TEST_CASE("Water - Potencial Lennard-Jones y total", "[Water]") {
    Atom* atoms = new Atom[2];
    atoms[0] = Atom(Vector(0.0, 0.0, 0.0), 1, 16.0, -0.834, 0.1, 3.5, 8);
    atoms[1] = Atom(Vector(0.757, 0.586, 0.0), 2, 1.0, 0.417, 0.05, 2.5, 1);
    Water w(1, atoms, 2);
    Vector box(10.0, 10.0, 10.0);

    SECTION("Debería calcular potencial Lennard-Jones") {
        Atom a(Vector(1.0, 0.0, 0.0), 3, 1.0, 0.417, 0.05, 2.5, 1);
        Real s = (3.5 + 2.5) / 2.0; // Lorentz-Berthelot
        Real e = std::sqrt(0.1 * 0.05);
        Real s_over_R = s / 1.0; // Distance = 1.0
        Real expected = 4 * e * (std::pow(s_over_R, 12) - std::pow(s_over_R, 6));
        REQUIRE(w.getLJPotential(a, s, e, box) == Approx(expected));
    }

    SECTION("Debería calcular potencial total con otra molécula de agua") {
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

    SECTION("Debería calcular potencial total con átomo") {
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
}
