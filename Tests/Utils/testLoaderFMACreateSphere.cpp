// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================

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
    const char* const Output = FParameters::getStr(argc,argv,"-f", "../Data/testSphere20k.fma");
    std::cout << "Creating : " << Output << "\n";

    // Create file
    std::ofstream myfile;
    myfile.open (Output);

    if(!myfile.is_open()){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }

    std::cout << "Creating " << NbParticles << " particles at " << Output << "\n";
    std::cout << "Working...\n";

    const FReal FRandMax = FReal(RAND_MAX);
    // Center of the box
    const FReal XCenter = 0.;
    const FReal YCenter = 0.;
    const FReal ZCenter = 0.;
    // Box width
    const FReal BoxWidth = 1.0;

    // System properties
    myfile << NbParticles << " " << BoxWidth << " " << XCenter << " " << YCenter << " " << ZCenter;

    if( FParameters::findParameter(argc,argv,"-double") == FParameters::NotFound){
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
    }
    else{
        const FReal rayon = FReal(0.2);
        const FReal thresh = FReal(0.03);
        const FReal threshDiv2 = thresh/2;

        const FReal offset = FReal(0.25);

        // Generate particles
        for( long idx = 0 ; idx < NbParticles/2 ; ++idx ){
            const FReal theta = (FReal(rand())/FRandMax) * FMath::FPi;
            const FReal omega = (FReal(rand())/FRandMax) * FMath::FPi * FReal(2);

            const FReal px = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter - offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            const FReal py = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter - offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            const FReal pz = rayon * FMath::Cos(theta) + ZCenter - offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;

            myfile << " \n" << px << " " << py << " " <<  pz << " 0.01";
        }

        for( long idx = 0 ; idx < NbParticles/2 ; ++idx ){
            const FReal theta = (FReal(rand())/FRandMax) * FMath::FPi;
            const FReal omega = (FReal(rand())/FRandMax) * FMath::FPi * FReal(2);

            const FReal px = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter + offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            const FReal py = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter + offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            const FReal pz = rayon * FMath::Cos(theta) + ZCenter + offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;

            myfile << " \n" << px << " " << py << " " <<  pz << " 0.01";
        }
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}



