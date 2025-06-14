#include "../cppambergromacs/ReaderFactory.hpp"
#include "../cppambergromacs/Configurations/Configuration.hpp"
#include <iostream>
using namespace std;

int main() {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("../archivos/jj.prmtop");
    vector<pair<int,string>> files= CoordinateReader::getFileIterator("../archivos/", "frame_*.pdb");
    string f= "../archivos/"+files[0].second;
    ToolKit::takeTime([&ti,&f]() {
        for(int i= 0; i < 1000; i++) {
            Configuration conf= Configuration(ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER), f, ti);
        }
    });

    return(0);
}