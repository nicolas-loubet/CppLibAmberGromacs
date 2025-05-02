#include <Vector.hpp>
#include <Particle.hpp>
#include <Atom.hpp>
#include <Molecule.hpp>
#include <iostream>
using namespace std;

int main(){
    Particle a;
    Particle b;
    Vector old_vectora(2,1,1);
    a.setPosition(old_vectora);
    Vector old_vectorb(5,5,2);
    b.setPosition(old_vectorb);
    Vector new_vector=a.getPosition();
    a.setID(10);
    a.setMass(12.01);
    //a.setZ(3);

    cout << "posicion: " << new_vector.x << endl;
    cout << "ID: " << a.getID() << endl;
    cout << "Mass: " << a.getMass() << endl;
    cout << "distance: " << a.distanceTo(b,{10,10,10}) << endl;
    
    Atom c;
    cout << "Z: " << c.getZ() << endl;
    cout << "s_AO: " << c.getSigma_AO() << endl;
    cout << "e_AO: " << c.getEpsilon_AO() << endl;
    cout << "charge: " << c.getCharge() << endl;

    Atom d({2.0,2.0,2.0}, 1);
    cout << "Z: " << d.getZ() << endl;
    cout << "s_AO: " << d.getSigma_AO() << endl;
    cout << "e_AO: " << d.getEpsilon_AO() << endl;
    cout << "charge: " << d.getCharge() << endl;

    Atom e({3.0,3.0,3.0}, 2, 1.0, 1.8, 2.0, 3.0, 2);
    cout << "Z: " << e.getZ() << endl;
    cout << "s_AO: " << e.getSigma_AO() << endl;
    cout << "e_AO: " << e.getEpsilon_AO() << endl;
    cout << "charge: " << e.getCharge() << endl;

    Atom* Atomos = new Atom[3]{c,d,e};
    Molecule f(1, Atomos, 3);
    Atom* conjunto=f.getAtoms();
    cout << "Z del atomo 1 de la molecula: " << conjunto[1].getZ() << endl;
    cout << "nAtoms: " << f.getNAtoms() << endl;
    return(0);
}