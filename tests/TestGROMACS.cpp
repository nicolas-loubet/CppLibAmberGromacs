#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <iostream>
using namespace std;

int main() {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS)->readTopology("../files/gromacs.top");
    string f= "../files/gromacs.gro";

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
    ToolKit::takeTime([&ti,&f,&cr]() {
        for(int i= 0; i < 5000; i++) {
            Configuration conf= Configuration(cr, f, ti);
            //vector<float> vis= conf.getInteractionsPerSite(ti.num_solutes+5);
            //cout << vis[0] << "\t" << vis[1] << "\t" << vis[2] << "\t" << vis[3] << endl;
        }
    });

    return(0);
}