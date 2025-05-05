#include "../cppambergromacs/ReaderFactory.hpp"
#include <iostream>
using namespace std;

int main() {
    ToolKit::takeTime([]() {
        for(int i= 0; i < 1; i++) //50000
            Configuration::TopolInfo topolInfo= ReaderFactory::createTopologyInfo(ReaderFactory::GROMACS, "../archivos/gromacs_topol.top");
    });

    return(0);
}