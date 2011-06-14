// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"
#include "../Src/Utils/FParameters.hpp"

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

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 9);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbPart = 2000000;//2000000
    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FOctree<FTestParticleTsm, FTestCellTsm, FVector, FTypedLeaf> tree(NbLevels, SizeSubLevels,1.0,F3DPosition(0.5,0.5,0.5));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating " << NbPart << " particles ..." << std::endl;
    counter.tic();
    {
        FTestParticleTsm particle;
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particle.setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
            if(rand() > RAND_MAX/2) particle.setAsTarget();
            else particle.setAsSource();

            tree.insert(particle);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    // FTestKernels FBasicKernels
    FTestKernels<FTestParticleTsm, FTestCellTsm, FVector> kernels;
    //FFmmAlgorithmTsm FFmmAlgorithmThreadTsm
    FFmmAlgorithmThreadTsm<FTestKernels, FTestParticleTsm, FTestCellTsm, FVector, FTypedLeaf> algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgo(&tree);


    return 0;
}


// [--LICENSE--]
