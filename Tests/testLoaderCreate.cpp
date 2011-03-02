// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// This file can generate basic particules files to load with basic loader
// g++ testLoaderCreate.cpp -o testLoaderCreate.exe

int main(int , char ** ){
    // Nb of particules
    const long NbParticules = 200000;

    // Center of the box
    const double XCenter = 0.5;
    const double YCenter = 0.5;
    const double ZCenter = 0.5;
    // Box width
    const double BoxWidth = 1.0;
    // Output file please let .temp extension
    const char * const Output = "testLoader.basic.temp";

    // Create file
    std::ofstream myfile;
    myfile.open (Output);

    if(!myfile.is_open()){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticules << " " << BoxWidth << " " << XCenter << " " << YCenter << " " << ZCenter;

    // Generate particules
    for( long idx = 0 ; idx < NbParticules ; ++idx ){
        const double px = ((double(rand())/RAND_MAX) * BoxWidth) + XCenter - (BoxWidth/2);
        const double py = ((double(rand())/RAND_MAX) * BoxWidth) + YCenter - (BoxWidth/2);
        const double pz = ((double(rand())/RAND_MAX) * BoxWidth) + ZCenter - (BoxWidth/2);

        myfile << " \n" << px << " " << py << " " <<  pz;
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--LICENSE--]
