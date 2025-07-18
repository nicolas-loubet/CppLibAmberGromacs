
#include "../cppambergromacs/CppLibAmberGromacs.hpp"

int main() {
    cout<<"0"<<endl;
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("step5_input.parm7");
    cout<<"1"<<endl;
    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    cout<<"2"<<endl;
    ConfigurationLipid c= ConfigurationLipid(cr, "frame1.pdb", ti);
    cout<<"3"<<endl;
    
    vector<float> resultado(20);
    for(int i=1; i<=ti.num_solutes;i++)
        {
        vector<float> resultado_molecula = c.orderParameter(i);
        for(int j=0; j<resultado_molecula.size();j++)
            {
            resultado[j]+=resultado_molecula[j]/ti.num_solutes;
            }
        }
    for(int i=0; i<resultado.size();i++)
        {
        cout<<"Atom C: "<<i<<" Order parameter: "<<resultado[i]<<endl;
        }
    cout<<"4"<<endl;
    return 0;
}