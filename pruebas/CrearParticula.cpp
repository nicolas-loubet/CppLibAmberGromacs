#include "../cppambergromacs/ReaderFactory.hpp"
#include <iostream>
using namespace std;

int main() {
    ToolKit::takeTime([]() {
        for(int i= 0; i < 50000; i++) //50000
            TopolInfo topolInfo= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS)->readTopology("../archivos/gromacs_topol.top");
    });

    return(0);
}