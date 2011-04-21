// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Components/FTypedLeaf.hpp"

#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Components/FTestParticle.hpp"
#include "../Sources/Components/FTestCell.hpp"
#include "../Sources/Components/FTestKernels.hpp"

#include "../Sources/Extenssions/FExtendParticleType.hpp"
#include "../Sources/Extenssions/FExtendCellType.hpp"

#include "../Sources/Core/FFmmAlgorithmTsm.hpp"
#include "../Sources/Core/FFmmAlgorithmThreadTsm.hpp"

#include "../Sources/Components/FBasicKernels.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

class FTestParticleTsm : public FTestParticle, public FExtendParticleType {
};

class FTestCellTsm: public FTestCell , public FExtendCellType{
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 10;//10;
    const int SizeSubLevels = 3;//3
    const long NbPart = 2000000;//2000000
    FTestParticleTsm* particles = new FTestParticleTsm[NbPart];
    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating " << NbPart << " particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        particles[idxPart].setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
        if(rand() > RAND_MAX/2) particles[idxPart].setAsTarget();
        else particles[idxPart].setAsSource();
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FOctree<FTestParticleTsm, FTestCellTsm, FTypedLeaf, NbLevels, SizeSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Inserting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        tree.insert(&particles[idxPart]);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    // FTestKernels FBasicKernels
    FTestKernels<FTestParticleTsm, FTestCellTsm, NbLevels> kernels;
    //FFmmAlgorithmTsm FFmmAlgorithmThreadTsm
    FFmmAlgorithmThreadTsm<FTestKernels, FTestParticleTsm, FTestCellTsm, FTypedLeaf, NbLevels, SizeSubLevels> algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgo(&tree);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Deleting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        particles[idxPart].~FTestParticleTsm();
    }
    delete [] particles;
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
