// [--License--]

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include <cmath>

#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FMath.hpp"

#include "../Src/Utils/FParameters.hpp"

// This file can generate basic particles files to load with basic loader

int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable can create a FMA particles files in a spherical scattering\n";
    std::cout << ">> You can pass a filename in parameter else the program will use\n";
    std::cout << ">> a default filename.\n";
    std::cout << ">> The format of the file is : \n";
    std::cout << ">> [number of particles] \n";
    std::cout << ">> [boxe width] [boxe x center] [boxe y center] [boxe z center]\n";
    std::cout << ">> [x] [y] [z] [physical value]...\n";
    //////////////////////////////////////////////////////////////

    // Nb of particles
    const long NbParticles = FParameters::getValue(argc,argv,"-nb", long(20000));

    // Center of the box
    const FReal XCenter = 0.;
    const FReal YCenter = 0.;
    const FReal ZCenter = 0.;
    // Box width
    const FReal BoxWidth = 1.0;
    // Output file please let .temp extension
    const char defaultFilename[] = "../Data/testSphere20k.fma";

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


    const FReal FRandMax = FReal(RAND_MAX);
    const FReal rayon = FReal(0.4);
    const FReal thresh = FReal(0.15);
    const FReal threshDiv2 = thresh/2;

    // Generate particles
    for( long idx = 0 ; idx < NbParticles ; ++idx ){
        const FReal theta = (FReal(rand())/FRandMax) * FMath::FPi;
        const FReal omega = (FReal(rand())/FRandMax) * FMath::FPi * FReal(2);

        const FReal px = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter + thresh * (FReal(rand())/FRandMax) - threshDiv2;
        const FReal py = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter + thresh * (FReal(rand())/FRandMax) - threshDiv2;
        const FReal pz = rayon * FMath::Cos(theta) + ZCenter + thresh * (FReal(rand())/FRandMax) - threshDiv2;

        myfile << " \n" << px << " " << py << " " <<  pz << " 0.01";
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--END--]
