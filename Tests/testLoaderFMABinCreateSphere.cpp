
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

#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FMath.hpp"
#include "../Src/Utils/FParameters.hpp"

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

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char* const Output = FParameters::getStr(argc,argv,"-f", "../Data/stest20k.bin.fma");
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


    FReal data[4];
    data[3] = FReal(0.1);

    if( !FParameters::findParameter(argc,argv,"-double") ){
        const FReal rayon = FReal(0.4);
        const FReal thresh = FReal(0.15);
        const FReal threshDiv2 = thresh/2;

        // Generate particles
        for( long idx = 0 ; idx < NbParticles ; ++idx ){
            const FReal theta = (FReal(rand())/FRandMax) * FMath::FPi;
            const FReal omega = (FReal(rand())/FRandMax) * FMath::FPi * FReal(2);

            data[0] = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            data[1] = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            data[2] = rayon * FMath::Cos(theta) + ZCenter + thresh * (FReal(rand())/FRandMax) - threshDiv2;

            fwrite(&data, sizeof(FReal), 4, myfile);
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

            data[0] = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter - offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            data[1] = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter - offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            data[2] = rayon * FMath::Cos(theta) + ZCenter - offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;

            fwrite(&data, sizeof(FReal), 4, myfile);
        }

        for( long idx = 0 ; idx < NbParticles/2 ; ++idx ){
            const FReal theta = (FReal(rand())/FRandMax) * FMath::FPi;
            const FReal omega = (FReal(rand())/FRandMax) * FMath::FPi * FReal(2);

            data[0] = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter + offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            data[1] = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter + offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;
            data[2] = rayon * FMath::Cos(theta) + ZCenter + offset + thresh * (FReal(rand())/FRandMax) - threshDiv2;

            fwrite(&data, sizeof(FReal), 4, myfile);
        }
    }

    fclose(myfile);

    std::cout << "Done\n";

    return 0;
}



