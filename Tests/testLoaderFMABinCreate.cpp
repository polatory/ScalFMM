// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../Src/Utils/FGlobal.hpp"
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
    const char* const Output = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma");
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
    // Generate particles
    for( FSize idx = 0 ; idx < NbParticles ; ++idx ){
        data[0] = ((FReal(rand())/FRandMax) * BoxWidth * 2) + XCenter - BoxWidth;
        data[1] = ((FReal(rand())/FRandMax) * BoxWidth * 2) + YCenter - BoxWidth;
        data[2] = ((FReal(rand())/FRandMax) * BoxWidth * 2) + ZCenter - BoxWidth;

        fwrite(&data, sizeof(FReal), 4, myfile);
    }

    fclose(myfile);

    std::cout << "Done\n";

    return 0;
}



