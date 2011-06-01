// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FList.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"


#include "../Src/Components/FBasicKernels.hpp"

// Compile by : g++ testFMMAlgorithm.cpp ../Src/Utils/FAssertable.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFMMAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 6;//10;
    const int SizeSubLevels = 3;//3
    const long NbPart = 200000;//2000000

    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FOctree<FTestParticle, FTestCell, FSimpleLeaf> tree(NbLevels, SizeSubLevels,1.0,F3DPosition(0.5,0.5,0.5));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    counter.tic();

    for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        FTestParticle particleToFill;
        particleToFill.setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
        tree.insert(particleToFill);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    // FTestKernels FBasicKernels
    FTestKernels<FTestParticle, FTestCell> kernels;
    //FFmmAlgorithm FFmmAlgorithmThread
    FFmmAlgorithmThread<FTestKernels, FTestParticle, FTestCell, FSimpleLeaf> algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgo(&tree);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
