#include "../cppambergromacs/CppLibAmberGromacs.hpp"
#include <iostream>
using namespace std;

int main() {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::GROMACS)->readTopology("../files/gromacs.top");
    string f= "../files/gromacs.gro";

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::GROMACS);
    ToolKit::takeTime([&ti,&f,&cr]() {
        for(int i= 0; i < 1; i++) {
            Configuration conf= Configuration(cr, f, ti);
            conf.classifyMolecules_includePentacoordinated();

            int cuenta_D3= 0, cuenta_D5= 0, cuenta_TA= 0, cuenta_TB= 0, cuenta_total= 0;
            for(int m= ti.num_solutes+1; m < conf.getNMolec(); m++) {
                if(conf.getMolec(m).getPosition().z > 5 &&
                   conf.getMolec(m).getPosition().z < conf.getBounds().z-5) continue;
                
                if(conf.getMolec(m).getClassification() == Configuration::CLASSIFICATION_D3_MOLECULE) cuenta_D3++;
                if(conf.getMolec(m).getClassification() == Configuration::CLASSIFICATION_D5_MOLECULE) cuenta_D5++;
                if(conf.getMolec(m).getClassification() == Configuration::CLASSIFICATION_TA_MOLECULE) cuenta_TA++;
                if(conf.getMolec(m).getClassification() == Configuration::CLASSIFICATION_TB_MOLECULE) cuenta_TB++;
                cuenta_total++;
            }

            cout << "f_D3=" << float(cuenta_D3)/float(cuenta_total)
                 << "\tf_D5=" << float(cuenta_D5)/float(cuenta_total)
                 << "\tf_TA=" << float(cuenta_TA)/float(cuenta_total)
                 << "\tf_TB=" << float(cuenta_TB)/float(cuenta_total) << endl;
        }
    });

    return(0);
}