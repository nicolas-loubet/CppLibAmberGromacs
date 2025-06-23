#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <iostream>
using namespace std;

int main() {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS)->readTopology("../files/lammps_system.data");
    string f= "../files/lammps_traj_minimized_9A.xyz";

    //cout << ti << endl;

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::LAMMPS);
    
    for(int i_conf= 0; i_conf < 5; i_conf++) {
        Configuration conf= Configuration(cr, f, ti);
        cr->incFrame();
        for(int i_molec= ti.num_solutes+1; i_molec <= conf.getNMolec(); i_molec++) {
            vector<float> vis= conf.getInteractionsPerSite(i_molec);
            cout << i_conf << "\t" << i_molec << "\t" << vis[0] << "\t" << vis[1] << "\t" << vis[2] << "\t" << vis[3] << endl;
        }
    }

    return 0;
}