
#include "../cppambergromacs/CppLibAmberGromacs.hpp"

int main() {

    TopolInfo ti1= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("step5_input.parm7");
    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    ConfigurationLipid c1= ConfigurationLipid(cr, "frame1.pdb", ti1);

    TopolInfo ti2= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("DMPC/jj.prmtop");
    ConfigurationLipid c2= ConfigurationLipid(cr, "DMPC/frame1.pdb", ti2);
    
    vector<vector<float>> resultado1(3,vector<float>(25,0.0f));
    vector<vector<float>> resultado2(3,vector<float>(25,0.0f));

    for(int i=1; i<=ti1.num_solutes;i++)
        {
        vector<vector<float>> resultado_molecula = c1.orderParameter(i);
        for(int j=0; j<resultado_molecula[1].size();j++)
            {
            resultado1[1][j]+=resultado_molecula[1][j]/ti1.num_solutes;
            resultado1[2][j]+=resultado_molecula[2][j]/ti1.num_solutes;
            //cout<<resultado_molecula[2][j]<<",";
            }
            //cout<<endl;
        }
    for(int i=1; i<=ti2.num_solutes;i++)
        {
        vector<vector<float>> resultado_molecula = c2.orderParameter(i);
        for(int j=0; j<resultado_molecula[1].size();j++)
            {
            resultado2[1][j]+=resultado_molecula[1][j]/ti2.num_solutes;
            resultado2[2][j]+=resultado_molecula[2][j]/ti2.num_solutes;
            }
        }
    /*
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology("DMPC/jj.prmtop");
    string file= "DMPC/frame1.pdb";
    ToolKit::takeTime([&ti,&file,&cr]() {
        vector<vector<float>> resultado(24);
        for(int i= 0; i < 100; i++) {
            ConfigurationLipid c= ConfigurationLipid(cr, file, ti);
            for(int i=1; i<=ti.num_solutes;i++)
                {
                vector<vector<float>> resultado_molecula = c.orderParameter(i);
                for(int j=0; j<resultado_molecula.size();j++)
                    {
                    resultado[1][j]=resultado_molecula[1][j]/ti.num_solutes;
                    resultado[2][j]=resultado_molecula[2][j]/ti.num_solutes;
                    }
                }
            
        }
    });*/
    

    ofstream f("parametro_de_orden.csv");
    f<<"Carbono,"<< "Orden_gel_SN1,"<<"Orden_gel_SN2,"<<"Orden_cristal_SN1," <<"Orden_cristal_SN2" <<endl; 
    for(int i=0; i<resultado1[1].size();i++)
        {
        //cout<<"Atom C: "<<i<<" Order parameter: "<<resultado1[i]<<endl;
        if(resultado1[1][i]!=0 || resultado2[1][i]!=0)
            {
            f<<i+1<<","<< resultado1[2][i] << "," << resultado1[1][i] << "," << resultado2[2][i] << "," << resultado2[1][i] << endl; 
            }
        }

    f.flush(); f.close();
    return 0;
}