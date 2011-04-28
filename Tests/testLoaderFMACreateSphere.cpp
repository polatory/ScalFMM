// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>

#include "../Src/Utils/FGlobal.hpp"

// This file can generate basic particles files to load with basic loader
// g++ testLoaderCreateSphere.cpp -O2 -o testLoaderCreateSphere.exe

int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable can create a FMA particles files in a spherical scattering";
    std::cout << ">> You can pass a filename in parameter else the program will use\n";
    std::cout << ">> a default filename.\n";
    std::cout << ">> The format of the file is : \n";
    std::cout << ">> [number of particles] \n";
    std::cout << ">> [boxe width] [boxe x center] [boxe y center] [boxe z center]\n";
    std::cout << ">> [x] [y] [z] [physical value]...\n";
    //////////////////////////////////////////////////////////////

    // Nb of particles
    const long NbParticles = 10;

    // Center of the box
    const float XCenter = 0.5;
    const float YCenter = 0.5;
    const float ZCenter = 0.5;
    // Box width
    const float BoxWidth = 1.0;
    // Output file please let .temp extension
    const char * const defaultFilename = "Sphere.fma";

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


    const float rayon = 0.4;
    const float thresh = 0.2;

    // Generate particles
    for( long idx = 0 ; idx < NbParticles ; ++idx ){
        const float phi = ((float(rand())/RAND_MAX) * thresh - (thresh/2)) + rayon;
        const float theta = M_PI * (float(rand())/RAND_MAX);
        const float omega = 2 * M_PI * (float(rand())/RAND_MAX);

        const float px = phi*cos(omega)*sin(theta) + XCenter;
        const float py = phi*sin(omega)*cos(theta) + YCenter;
        const float pz = phi*cos(theta) + ZCenter;

        myfile << " \n" << px << " " << py << " " <<  pz << " 0.01";
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}


// [--LICENSE--]
