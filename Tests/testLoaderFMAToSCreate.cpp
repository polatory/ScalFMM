// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Sources/Utils/FGlobal.hpp"

// This file can generate basic particules files in the FMA format
// g++ testLoaderFMAToSCreate.cpp -o testLoaderFMAToSCreate.exe

int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable can create a FMA-like particules files (this is not a fma file!)";
    std::cout << ">> You can pass a filename in parameter else the program will use\n";
    std::cout << ">> a default filename.\n";
    std::cout << ">> The format of the file is : \n";
    std::cout << ">> [number of particules] \n";
    std::cout << ">> [boxe width] [boxe x center] [boxe y center] [boxe z center]\n";
    std::cout << ">> [x] [y] [z] [physical value] [1 if target 0 if source]...\n";
    //////////////////////////////////////////////////////////////
    // Nb of particules
    const long NbParticules = 200000;

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char * const defaultFilename = "testLoaderFMA.tor.fma";

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

    std::cout << "Generating " << NbParticules << " in " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticules << "\n";
    myfile << BoxWidth << "\t" << XCenter << "\t" << YCenter << "\t" << ZCenter;

    // Generate particules
    for( long idx = 0 ; idx < NbParticules ; ++idx ){
        const FReal px = ((FReal(rand())/RAND_MAX) * BoxWidth * 2) + XCenter - BoxWidth;
        const FReal py = ((FReal(rand())/RAND_MAX) * BoxWidth * 2) + YCenter - BoxWidth;
        const FReal pz = ((FReal(rand())/RAND_MAX) * BoxWidth * 2) + ZCenter - BoxWidth;

        const int isTarget = rand() > RAND_MAX/2 ? 1 : 0;

        myfile << "\n" << px << "\t" << py << "\t" <<  pz << "\t" << (0.01) << "\t" << isTarget;
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--LICENSE--]
