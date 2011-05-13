// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FList.hpp"

#include "../Src/Components/FTypedLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"

#include "../Src/Extensions/FExtendParticleType.hpp"
#include "../Src/Extensions/FExtendCellType.hpp"

#include "../Src/Core/FFmmAlgorithmTsm.hpp"
#include "../Src/Core/FFmmAlgorithmThreadTsm.hpp"

#include "../Src/Components/FBasicKernels.hpp"

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
    std::cout << "Done  " << "(@Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

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
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

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
