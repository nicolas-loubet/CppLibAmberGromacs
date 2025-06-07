#ifndef AMBER_READERS_HPP
#define AMBER_READERS_HPP

/**
 * Version: June 2025
 * Author: Ezequiel Cuenca
 */

#include "ReaderInterfaces.hpp"
#include <string>
#include <fstream>
#include <sstream>

class AmberTopologyReader : public TopologyReader {
    private:
        /**
         * Given a file, find the position of the flag in the file.
         * The returned map contains the flag name as the key and the position as the value.
         * The position is the byte position of the first character after the flag.
         * @param file The file to search
         * @return A map with the position of the flag.
         */
        static map<string,int> flag_position(ifstream &file)
        {
            string line= "";
            map<string,int> flag;
            while(getline(file, line)) {
                if(line.find("%FLAG") != string::npos) {
                    line= ToolKit::strip(line.substr(5, 80));
                    flag[line]= file.tellg();
                    flag[line]-= 81;
                }
            }
            return(flag);
        }

    public:
        AmberTopologyReader()= default;
        TopolInfo readTopology(const std::string& filename) const override {
            // Implementación para leer prmtop de AMBER-EZEQUIEL
            return TopolInfo();
        }
};

class AmberCoordinateReader : public CoordinateReader {
    public:
        AmberCoordinateReader()= default;
        bool readCoordinates(const string& filename, const TopolInfo& topol_info, Molecule** molecs) const override {
            // Implementación para leer mdcrd de AMBER
            // Los datos se guardan en el Molecule**
            return true;
        }
};

#endif
