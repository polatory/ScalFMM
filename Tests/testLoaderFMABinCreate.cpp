// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <cstdio>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FParameters.hpp"

// This file can generate basic particles files in the FMA format
// g++ testLoaderFMACreate.cpp -o testLoaderFMACreate.exe

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

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char defaultFilename[] = "testLoaderFMA.fma";

    const char* Output;

    if(argc == 1){
        std::cout << "You have to give a filename in argument.\n";
        std::cout << "The program will create one with a default name : " << defaultFilename << "\n";
        Output = defaultFilename;
    }
    else{
        Output = argv[1];
        std::cout << "Creating : " << Output << "\n";
    }

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

    FReal data[4];
    data[3] = FReal(0.1);
    // Generate particles
    for( FSize idx = 0 ; idx < NbParticles ; ++idx ){
        data[0] = ((FReal(rand())/FRandMax) * BoxWidth * 2) + XCenter - BoxWidth;
        data[1] = ((FReal(rand())/FRandMax) * BoxWidth * 2) + YCenter - BoxWidth;
        data[2] = ((FReal(rand())/FRandMax) * BoxWidth * 2) + ZCenter - BoxWidth;

        /*data[0] = ((FReal(idx)/NbParticles) * BoxWidth * 2) + XCenter - BoxWidth;
        data[1] = ((FReal(idx)/NbParticles) * BoxWidth * 2) + YCenter - BoxWidth;
        data[2] = ((FReal(idx)/NbParticles) * BoxWidth * 2) + ZCenter - BoxWidth;*/

        fwrite(&data, sizeof(FReal), 4, myfile);
    }

    fclose(myfile);

    std::cout << "Done\n";

    return 0;
}


// [--LICENSE--]
