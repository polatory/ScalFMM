// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FParameters.hpp"

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
    const long NbParticles = FParameters::getValue(argc,argv,"-nb", long(20000));
    const FReal physicalValue = FParameters::getValue(argc,argv,"-pv", FReal(0.1));

    // Box width
    const FReal BoxWidth = FParameters::getValue(argc,argv,"-width", FReal(1.0/2.0));

    // Center of the box
    const FReal XCenter = BoxWidth;
    const FReal YCenter = BoxWidth;
    const FReal ZCenter = BoxWidth;

    // Output file please let .temp extension
    const char* const Output = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    std::cout << "Creating : " << Output << "\n";

    // Create file
    std::ofstream myfile;
    myfile.open (Output);

    if(!myfile.is_open()){
        std::cout << "Cannot create " << Output << "\n";
        return 1;
    }

    std::cout << "Generating " << NbParticles << " in " << Output << "\n";
    std::cout << "Working...\n";

    // System properties
    myfile << NbParticles << "\n";
    myfile << BoxWidth << "\t" << XCenter << "\t" << YCenter << "\t" << ZCenter;

    const FReal Limite = FReal(0.00001);
    const FReal LimitedBoxWidth = BoxWidth - Limite * 2;

    // Generate particles
    for( long idx = 0 ; idx < NbParticles ; ++idx ){
        const FReal px = ((FReal(drand48())) * LimitedBoxWidth * 2) + XCenter - BoxWidth + Limite;
        const FReal py = ((FReal(drand48())) * LimitedBoxWidth * 2) + YCenter - BoxWidth + Limite;
        const FReal pz = ((FReal(drand48())) * LimitedBoxWidth * 2) + ZCenter - BoxWidth + Limite;

        myfile << "\n" << px << "\t" << py << "\t" <<  pz << "\t" << physicalValue;
    }

    myfile.close();

    std::cout << "Done\n";

    return 0;
}



