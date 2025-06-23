#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <iostream>
using namespace std;

int main() {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::LAMMPS)->readTopology("../files/lammps_system.data");
    string f= "../files/lammps_traj_minimized_9A.xyz";

    

    return 0;
}