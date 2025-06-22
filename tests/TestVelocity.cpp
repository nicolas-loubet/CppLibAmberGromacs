#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <iostream>
using namespace std;

int main() {
    /*TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("../files/amber.prmtop");
    vector<pair<int,string>> files= CoordinateReader::getFileIterator("../files/", "amber_frame_*.pdb");
    string f= "../files/"+files[0].second;

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    ToolKit::takeTime([&ti,&f,&cr]() {
        for(int i= 0; i < 1; i++) {
            Configuration conf= Configuration(cr, f, ti);
            vector<float> vis= conf.getInteractionsPerSite(13722);
            cout << vis[0] << "\t" << vis[1] << "\t" << vis[2] << "\t" << vis[3] << endl;
        }
    });*/

    cout << "Leo ti" << endl;
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("../files/control_bulk_amber.prmtop");
    cout << "Leido" << endl;
    string f= "../files/control_bulk_amber.pdb";

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    ToolKit::takeTime([&ti,&f,&cr]() {
        for(int i= 0; i < 1; i++) {
            cout << "Leo conf" << endl;
            Configuration conf= Configuration(cr, f, ti);
            cout << "Leida" << endl;
            vector<float> vis= conf.getInteractionsPerSite(1);
            cout << vis[0] << "\t" << vis[1] << "\t" << vis[2] << "\t" << vis[3] << endl;
        }
    });

    return(0);
}