
#include <CppLibAmberGromacs.hpp>

void calcularOrden(string carpeta, vector<vector<float>>& resultado) 
    {
    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology(carpeta+"step5_input.parm7");
    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);
    
    for(int i=1;i<=1501;i++)
        {
        ConfigurationLipid c= ConfigurationLipid(cr, carpeta+"separados/frame"+to_string(i)+".pdb", ti);
        for(int i=1; i<=ti.num_solutes;i++)
            {
            vector<vector<float>> resultado_molecula = c.orderParameter(i);
            for(int j=0; j<resultado_molecula[1].size();j++)
                {
                resultado[1][j]+=resultado_molecula[1][j]/ti.num_solutes;
                resultado[2][j]+=resultado_molecula[2][j]/ti.num_solutes;
                }
            }
        }
    
    for(int i=1; i<=resultado[1].size();i++)
        {
        resultado[1][i]=resultado[1][i]/(1501);
        resultado[2][i]=resultado[2][i]/(1501);
        }
    }
void escribirOrden(vector<vector<float>>& resultado)
    {
    ofstream f("parametro_de_orden_DMPC.csv");
    f<<"Carbono,"<< "Orden_Parameter_SN2," << "Orden_Parameter_SN1" <<endl; 
    for(int i=0; i<resultado[1].size();i++)
        {
        if(resultado[1][i]!=0)
            {
            f<<i+1<<","<< resultado[1][i] << ","<< resultado[2][i] <<endl; 
            }
        }
    f.flush(); f.close();
    }

int main() {
    
    ToolKit::takeTime([]() {
        vector<vector<float>> resultado(3,vector<float>(23,0.0f));
        //vector<string> carpeta= {"../DMPC/LIPID21/"};
        string carpeta="../../DMPC/LIPID21/";
        calcularOrden(carpeta, resultado);
        //ToolKit::parallel(calcularOrden,carpetas,resultado);
        escribirOrden(resultado);
    });
    
 
    return 0;
}
