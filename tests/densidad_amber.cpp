
#include <CppLibAmberGromacs.hpp>

string getPath(const string path) { return "../"+path+"/"; }
string getTopology(const string path) { return getPath(path)+"jj.prmtop"; }
string getRealConfigurationPath(const string path) { return getPath(path)+"separados/"; }
string getMinConfigurationPath(const string path) { return getPath(path)+"min/"; }

struct Salida {
    float** densidad;
    int contador_moleculas;
};

struct LimitesZ {
    float minz, maxz;
    float boundz;
};

bool estaDentro(float z, LimitesZ lim) {
    if(lim.minz < lim.maxz)
        return lim.minz <= z && z <= lim.maxz;
    else {
        if(z < 0)
            z+= lim.boundz;
        else if(z > lim.boundz)
            z-= lim.boundz;

        return (lim.minz <= z && z <= lim.boundz) || (0 <= z && z <= lim.maxz);
    }
}

LimitesZ calcularLimitesZ(Configuration& c, const float NUM_SOLUTES) {
    LimitesZ lim;
    lim.boundz= c.getBounds().z;
    lim.minz= lim.boundz; lim.maxz= 0;
    for(int m= 1; m <= NUM_SOLUTES; m++) {
        for(int a= 1; a <= c.getMolec(m).getNAtoms(); a++) {
            Atom atom= c.getMolec(m).getAtoms()[a-1];
            if(atom.getPosition().z < lim.minz) lim.minz= atom.getPosition().z;
            if(atom.getPosition().z > lim.maxz) lim.maxz= atom.getPosition().z;
        }
    }
    lim.minz-= 25;
    lim.maxz+= 25;
    return lim;
}

float centro_de_masa_z_molec(Configuration& conf,const int numero_lipidos){
    float z_medio=0.0;
    for(int _i = 1; _i<=numero_lipidos; _i++)
        {
        Vector vector_lipido= conf.getMolec(_i).getPosition();
        float z_lipid=vector_lipido.z;
        z_medio+=z_lipid;
        }
    z_medio /= numero_lipidos;
    return z_medio;
}

float centro_de_masa_z_atom(Configuration& conf,const int numero_lipidos){
    float z_medio=0.0;
    int number_of_total_atom=0;
    for(int _i = 1; _i<=numero_lipidos; _i++)
        {
        int number_of_atom=conf.getMolec(_i).getNAtoms();
        for(int _j = 1; _j<=number_of_atom; _j++)
            {
            Vector vector_atomo= conf.getMolec(_i).getAtom(_j).getPosition();
            float z_atom=vector_atomo.z;
            z_medio+=z_atom;
            number_of_total_atom+=1;
            }
        }
    z_medio /= number_of_total_atom;
    return z_medio;
}


void calculoPrincipal(string carpeta, Salida& salida) {
    string carpeta_configuraciones= getRealConfigurationPath(carpeta);
    cout << "Estudiando carpeta: "+carpeta+"\n" << flush;

    TopolInfo ti= ReaderFactory::createTopologyReader(ReaderFactory::ProgramFormat::AMBER)->readTopology(getTopology(carpeta));
    vector<pair<int,string>> files= CoordinateReader::getFileIterator(carpeta_configuraciones, "frame*.pdb");
    
    const int N_BINS= 2000;
    const float MIN_BINS= 0;
    const float MAX_BINS= 100;

    salida.densidad= new float*[5];
    for(int group= 0; group < 5; group++) {
        salida.densidad[group]= new float[N_BINS];
        for(int i= 0; i < N_BINS; i++)
            salida.densidad[group][i]= 0;
    }
    salida.contador_moleculas= 0;

    CoordinateReader* cr= ReaderFactory::createCoordinateReader(ReaderFactory::ProgramFormat::AMBER);

    for(pair<int,string> f: files) {
        Configuration c= Configuration(cr,carpeta_configuraciones+f.second,ti);
        LimitesZ lim= calcularLimitesZ(c,ti.num_solutes);
        //cout << "    frame "+to_string(f.first)+": "+f.second+"\n" << flush;

        float centro_z = centro_de_masa_z_molec(c, ti.num_solutes);
        
        vector<float> bins_agua(N_BINS);
        vector<float> bins_CO_sn1(N_BINS);
        vector<float> bins_CO_sn2(N_BINS);
        vector<float> bins_P(N_BINS);
        vector<float> bins_N(N_BINS);

        float diferencia=(MAX_BINS-MIN_BINS)/N_BINS;
        float volumen=diferencia*c.getBounds().x*c.getBounds().y;
        float factor_volumen=1/(2*6.02e-1*volumen);

        for(int a =ti.num_solutes; a<=c.getNMolec();a++)
            {
            Vector vector_agua = c.getMolec(a).getPosition();
            float z=vector_agua.z-centro_z;
            float z_PBC= abs(z - c.getBounds().z*round(z/c.getBounds().z));  
            
            int pos=ToolKit::getBinPosition(z_PBC ,MIN_BINS,MAX_BINS,N_BINS,true);
            bins_agua[pos]+=1;  
            }

        for(int a =1; a<=c.getNMolec();a++)
            {
            int n_atoms=c.getMolec(a).getNAtoms();
            for(int b =1;b<=n_atoms;b++)
                {
                string valor = get<1>(ti.atom_type_name_charge_mass[a-1].at(b-1));
                if(valor=="C32")
                    {
                    Vector vector_c32 = c.getMolec(a).getAtom(b).getPosition();
                    float z=vector_c32.z-centro_z;
                    float z_PBC= abs(z - c.getBounds().z*round(z/c.getBounds().z));  
                    int pos=ToolKit::getBinPosition(z_PBC ,MIN_BINS,MAX_BINS,N_BINS,true);
                    bins_CO_sn2[pos]+=1; 
                    }    
                if(valor=="C15")
                    {
                    Vector vector_c15 = c.getMolec(a).getAtom(b).getPosition();
                    float z=vector_c15.z-centro_z;
                    float z_PBC= abs(z - c.getBounds().z*round(z/c.getBounds().z));  
                    int pos=ToolKit::getBinPosition(z_PBC ,MIN_BINS,MAX_BINS,N_BINS,true);
                    bins_CO_sn1[pos]+=1; 
                    }    
                if(valor=="P8")
                    {
                    Vector vector_P8 = c.getMolec(a).getAtom(b).getPosition();
                    float z=vector_P8.z-centro_z;
                    float z_PBC= abs(z - c.getBounds().z*round(z/c.getBounds().z));  
                    int pos=ToolKit::getBinPosition(z_PBC ,MIN_BINS,MAX_BINS,N_BINS,true);
                    bins_P[pos]+=1; 
                    } 
                if(valor=="N4")
                    {
                    Vector vector_N4 = c.getMolec(a).getAtom(b).getPosition();
                    float z=vector_N4.z-centro_z;
                    float z_PBC= abs(z - c.getBounds().z*round(z/c.getBounds().z));  
                    int pos=ToolKit::getBinPosition(z_PBC ,MIN_BINS,MAX_BINS,N_BINS,true);
                    bins_N[pos]+=1; 
                    } 
                
                }
            }
        for(int a =0; a<N_BINS;a++)
            {
            salida.densidad[0][a]+=bins_agua[a]*18.02*factor_volumen/files.size();
            salida.densidad[1][a]+=bins_CO_sn1[a]*(12.011+15.9994*2)*factor_volumen/files.size();
            salida.densidad[2][a]+=bins_CO_sn2[a]*(12.011+15.9994*2)*factor_volumen/files.size();
            salida.densidad[3][a]+=bins_P[a]*(30.974+15.9994*4)*factor_volumen/files.size();
            salida.densidad[4][a]+=bins_N[a]*(14.007+3*(12.0110+3*1.0080))*factor_volumen/files.size();
            }
            /*
        for(int a =0; a<N_BINS;a++)
            {
            salida.densidad[0][a]/=files.size();
            salida.densidad[1][a]/=files.size();
            salida.densidad[2][a]/=files.size();
            salida.densidad[3][a]/=files.size();
            salida.densidad[4][a]/=files.size();
            }   
            */
    }
}

void escribirResultado(vector<Salida> salidas) {

    cout << "Tamaño de salidas: " << salidas.size() << endl;
    for (int rep = 0; rep < salidas.size(); rep++) {
        cout << "salidas[" << rep << "].vis_bins es " << (salidas[rep].densidad == nullptr ? "nullptr" : "OK") << endl;
    }
    const int N_BINS= 2000;
    const float MIN_BINS= 0;
    const float MAX_BINS= 100;
    const float DX= (MAX_BINS-MIN_BINS)/N_BINS;

    float** densidad_promedio= new float*[5];
    for(int group= 0; group < 5; group++) {
        densidad_promedio[group]= new float[N_BINS];
        for(int i= 0; i < N_BINS; i++) {
            densidad_promedio[group][i]= 0.0;
            for(int rep= 0; rep < salidas.size(); rep++){
                if (salidas[rep].densidad == nullptr) {
                cerr << "Error: vis_bins de salidas[" << rep << "] no está inicializado." << endl;
                exit(EXIT_FAILURE);
                }
                densidad_promedio[group][i]+= salidas[rep].densidad[group][i];
            }
        }
    }
    for(int group= 0; group < 5; group++)
        {
        for(int i= 0; i < N_BINS; i++)
            {
            densidad_promedio[group][i]/= salidas.size();
            }
        }

    ofstream f("densidad.csv");
    const string CSV_SEP= ",";
    f<<"bin" <<","<< "water"<<","<< "carboni1"<<","<< "carboni2"<<"," << "phosp"<<","<< "nitrogen"<<endl; 
    
    for(int i= 0; i < N_BINS; i++) {
        f << MIN_BINS + i*DX + DX/2;
        for(int group= 0; group < 5; group++)
            f << CSV_SEP << densidad_promedio[group][i];
        f << endl;
    }

    f.flush(); f.close();
}

int main() {
    ToolKit::takeTime([]() {
        //vector<string> carpetas= {"DMPC/dmpc323-3"};
        vector<string> carpetas= {"DMPC/dmpc323-1","DMPC/dmpc323-2","DMPC/dmpc323-3"};
        vector<Salida> salidas(carpetas.size());
        ToolKit::parallel(calculoPrincipal, carpetas, salidas);
        
        escribirResultado(salidas);
    });
    return 0;
}