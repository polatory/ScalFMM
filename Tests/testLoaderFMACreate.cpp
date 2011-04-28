// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FGlobal.hpp"

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
    const long NbParticles = 20;

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char * const defaultFilename = "testLoaderFMA.fma";

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
    std::ofstream myfile;
    myfile.open (Output);

    if(!myfile.is_open()){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }

    std::cout << "Generating " << NbParticles << " in " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticles << "\n";
    myfile << BoxWidth << "\t" << XCenter << "\t" << YCenter << "\t" << ZCenter;

    // Generate particles
    for( long idx = 0 ; idx < NbParticles ; ++idx ){
        const FReal px = ((FReal(rand())/RAND_MAX) * BoxWidth * 2) + XCenter - BoxWidth;
        const FReal py = ((FReal(rand())/RAND_MAX) * BoxWidth * 2) + YCenter - BoxWidth;
        const FReal pz = ((FReal(rand())/RAND_MAX) * BoxWidth * 2) + ZCenter - BoxWidth;

        myfile << "\n" << px << "\t" << py << "\t" <<  pz << "\t" << (0.01);
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--LICENSE--]
