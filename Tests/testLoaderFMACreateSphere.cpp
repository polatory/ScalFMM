// /!\ Please, you must read the license at the bottom of this page

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <cmath>

// This file can generate basic particules files to load with basic loader
// g++ testLoaderCreateSphere.cpp -O2 -o testLoaderCreateSphere.exe

int main(int , char ** ){
    // Nb of particules
    const long NbParticules = 2000000;

    // Center of the box
    const float XCenter = 0.5;
    const float YCenter = 0.5;
    const float ZCenter = 0.5;
    // Box width
    const float BoxWidth = 1.0;
    // Output file please let .temp extension
    const char * const Output = "Sphere.fma";

    // Create file
    std::ofstream myfile;
    myfile.open (Output);

    if(!myfile.is_open()){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }

    std::cout << "Creating " << NbParticules << " particules at " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticules << " " << BoxWidth << " " << XCenter << " " << YCenter << " " << ZCenter;


    const float rayon = 0.4;
    const float thresh = 0.2;

    // Generate particules
    for( long idx = 0 ; idx < NbParticules ; ++idx ){
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
