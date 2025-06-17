#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <iostream>
using namespace std;

int main() {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("../files/amber.prmtop");
    vector<pair<int,string>> files= CoordinateReader::getFileIterator("../files/", "amber_frame_*.pdb");
    string f= "../files/"+files[0].second;

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    ToolKit::takeTime([&ti,&f,&cr]() {
        for(int i= 0; i < 100; i++) {
            Configuration conf= Configuration(cr, f, ti);
        }
    });

    return(0);
}