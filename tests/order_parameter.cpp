
#include "../cppambergromacs/CppLibAmberGromacs.hpp"

int main() {

    TopolInfo ti1= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("step5_input.parm7");
    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    ConfigurationLipid c1= ConfigurationLipid(cr, "frame1.pdb", ti1);

    TopolInfo ti2= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("DMPC/jj.prmtop");
    ConfigurationLipid c2= ConfigurationLipid(cr, "DMPC/frame1.pdb", ti2);
    
    vector<float> resultado1(23);
    vector<float> resultado2(23);

    for(int i=1; i<=ti1.num_solutes;i++)
        {
        vector<float> resultado_molecula = c1.orderParameter(i);
        for(int j=0; j<resultado_molecula.size();j++)
            {
            resultado1[j]+=resultado_molecula[j]/ti1.num_solutes;
            }
        }
    for(int i=1; i<=ti2.num_solutes;i++)
        {
        vector<float> resultado_molecula = c2.orderParameter(i);
        for(int j=0; j<resultado_molecula.size();j++)
            {
            resultado2[j]+=resultado_molecula[j]/ti2.num_solutes;
            }
        }
    
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("DMPC/jj.prmtop");
    string file= "DMPC/frame1.pdb";
    ToolKit::takeTime([&ti,&file,&cr]() {
        vector<float> resultado(24);
        for(int i= 0; i < 100; i++) {
            ConfigurationLipid c= ConfigurationLipid(cr, file, ti);
            for(int i=1; i<=ti.num_solutes;i++)
                {
                vector<float> resultado_molecula = c.orderParameter(i);
                for(int j=0; j<resultado_molecula.size();j++)
                    {
                    resultado[j]=resultado_molecula[j]/ti.num_solutes;
                    }
                }
            
        }
    });
    


    ofstream f("parametro_de_orden.csv");
    f<<"Carbono,"<< "Orden_gel,"<<"Orden_cristal" <<endl; 
    for(int i=0; i<resultado1.size();i++)
        {
        //cout<<"Atom C: "<<i<<" Order parameter: "<<resultado1[i]<<endl;
        if(resultado1[i]!=0 || resultado2[i]!=0)
            {
            f<<i+1<<","<< resultado1[i] << "," << resultado2[i] << endl; 
            }
        }

    f.flush(); f.close();
    return 0;
}