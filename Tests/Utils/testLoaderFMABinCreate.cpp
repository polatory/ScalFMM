// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FParameters.hpp"

// This file can generate basic particles files in the FMA format

int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable can create a FMA particles files\n";
    std::cout << ">> You can pass a filename in parameter else the program will use\n";
    std::cout << ">> a default filename.\n";
    std::cout << ">> The format of the file is : \n";
    std::cout << ">> [number of particles] \n";
    std::cout << ">> [boxe width] [boxe x center] [boxe y center] [boxe z center]\n";
    std::cout << ">> [x] [y] [z] [physical value]...\n";
    //////////////////////////////////////////////////////////////

    // Nb of particles
    const FSize NbParticles = FParameters::getValue(argc,argv,"-nb", FSize(20000));
    const FReal physicalValue = FParameters::getValue(argc,argv,"-pv", FReal(0.1));
    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
                                    "../../Data/test20k.bin.fma.single":
                                    "../../Data/test20k.bin.fma.double";
    const char* const Output = FParameters::getStr(argc,argv,"-f", defaultFilename);
    std::cout << "Creating : " << Output << "\n";

    // Create file
    FILE* const myfile = fopen(Output, "wb");
    if(!myfile){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }

    std::cout << "Generating " << NbParticles << " in " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    const int sizeOfFreal = int(sizeof(FReal));
    const FReal FRandMax = FReal(RAND_MAX);

    fwrite(&sizeOfFreal, sizeof(int),   1, myfile);
    fwrite(&NbParticles, sizeof(FSize), 1, myfile);

    fwrite(&BoxWidth,   sizeof(FReal), 1, myfile);
    fwrite(&XCenter,    sizeof(FReal), 1, myfile);
    fwrite(&YCenter,    sizeof(FReal), 1, myfile);
    fwrite(&ZCenter,    sizeof(FReal), 1, myfile);

    const FReal Limite = FReal(0.00001);
    const FReal LimitedBoxWidth = BoxWidth - Limite * 2;

    FReal data[4];
    data[3] = physicalValue;
    // Generate particles
    for( FSize idx = 0 ; idx < NbParticles ; ++idx ){
        data[0] = ((FReal(rand())/FRandMax) * LimitedBoxWidth * 2) + XCenter - BoxWidth + Limite;
        data[1] = ((FReal(rand())/FRandMax) * LimitedBoxWidth * 2) + YCenter - BoxWidth + Limite;
        data[2] = ((FReal(rand())/FRandMax) * LimitedBoxWidth * 2) + ZCenter - BoxWidth + Limite;

        fwrite(&data, sizeof(FReal), 4, myfile);
    }

    fclose(myfile);

    std::cout << "Done\n";

    return 0;
}



