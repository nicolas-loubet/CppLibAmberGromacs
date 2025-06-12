#include "../cppambergromacs/ReaderFactory.hpp"
#include "../cppambergromacs/Configurations/Configuration.hpp"
#include <iostream>
using namespace std;

int main() {
    ToolKit::takeTime([]() {
        for(int i= 0; i < 1; i++) {
            TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("../archivos/jj.prmtop");
            vector<pair<int,string>> files= CoordinateReader::getFileIterator("../archivos/", "frame_*.pdb");
            for(auto& f: files)
                cout << f.first << " " << f.second << endl;
            string f= files[0].second;
            Configuration conf= Configuration(ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER), f, ti);
            cout << conf.getNMolec() << endl;
        }
    });

    return(0);
}