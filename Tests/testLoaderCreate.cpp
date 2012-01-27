// [--License--]

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FParameters.hpp"

// This file can generate basic particles files to load with basic loader
// g++ testLoaderCreate.cpp -o testLoaderCreate.exe

int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable can create a basic particles files in custom format\n";
    std::cout << ">> You can pass a filename in parameter else the program will use\n";
    std::cout << ">> a default filename.\n";
    std::cout << ">> The format of the file is : \n";
    std::cout << ">> [number of particles] [boxe width] [boxe x center] [boxe y center] [boxe z center]\n";
    std::cout << ">> [x] [y] [z]...\n";
    //////////////////////////////////////////////////////////////

    // Nb of particles
    const FSize NbParticles = FParameters::getValue(argc,argv,"-nb", FSize(20000));
    const FReal FRandMax = FReal(RAND_MAX);

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;
    // Box width
    const FReal BoxWidth = 1.0;
    // Output file please let .temp extension
    const char defaultFilename[] = "../Data/test20k.basic";
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

    std::cout << "Creating " << NbParticles << " particles at " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticles << " " << BoxWidth << " " << XCenter << " " << YCenter << " " << ZCenter;

    // Generate particles
    for( long idx = 0 ; idx < NbParticles ; ++idx ){
        const FReal px = ((FReal(rand())/FRandMax) * BoxWidth) + XCenter - (BoxWidth/FReal(2.0));
        const FReal py = ((FReal(rand())/FRandMax) * BoxWidth) + YCenter - (BoxWidth/FReal(2.0));
        const FReal pz = ((FReal(rand())/FRandMax) * BoxWidth) + ZCenter - (BoxWidth/FReal(2.0));

        myfile << " \n" << px << " " << py << " " <<  pz;
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--END--]
