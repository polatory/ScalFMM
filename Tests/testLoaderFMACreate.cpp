// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// This file can generate basic particules files in the FMA format
// g++ testLoaderFMACreate.cpp -o testLoaderFMACreate.exe

int main(int , char ** ){
    // Nb of particules
    const long NbParticules = 200000;

    // Center of the box
    const FReal XCenter = 0.5;
    const FReal YCenter = 0.5;
    const FReal ZCenter = 0.5;

    // Box width
    const FReal BoxWidth = 1.0/2;
    // Output file please let .temp extension
    const char * const Output = "testLoaderFMA.fma";

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

        myfile << "\n" << px << "\t" << py << "\t" <<  pz << "\t" << (0.01);
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--LICENSE--]
