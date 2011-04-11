// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Components/FTypedLeaf.hpp"

#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Components/FTestParticule.hpp"
#include "../Sources/Components/FTestCell.hpp"
#include "../Sources/Components/FTestKernels.hpp"

#include "../Sources/Extenssions/FExtendParticuleType.hpp"
#include "../Sources/Extenssions/FExtendCellType.hpp"

#include "../Sources/Core/FFmmAlgorithmTsm.hpp"
#include "../Sources/Core/FFmmAlgorithm.hpp"
#include "../Sources/Core/FFmmAlgorithmArray.hpp"
#include "../Sources/Core/FFmmAlgorithmThreaded.hpp"
#include "../Sources/Core/FFmmAlgorithmTask.hpp"

#include "../Sources/Components/FBasicKernels.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particules is impacted each other particules
  */

class FTestParticuleTsm : public FTestParticule, public FExtendParticuleType {
};

class FTestCellTsm: public FTestCell , public FExtendCellType{
};


// Simply create particules and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 10;//10;
    const int SizeSubLevels = 3;//3
    const long NbPart = 2000000;//2000000
    FTestParticuleTsm* particules = new FTestParticuleTsm[NbPart];
    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating " << NbPart << " particules ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        particules[idxPart].setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
        if(rand() > RAND_MAX/2) particules[idxPart].setAsTarget();
        else particules[idxPart].setAsSource();
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FOctree<FTestParticuleTsm, FTestCellTsm, FTypedLeaf, NbLevels, SizeSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Inserting particules ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        tree.insert(&particules[idxPart]);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particules ..." << std::endl;
    counter.tic();

    // FTestKernels FBasicKernels
    FTestKernels<FTestParticuleTsm, FTestCellTsm, NbLevels> kernels;
    //FFMMAlgorithm FFMMAlgorithmThreaded FFMMAlgorithmArray FFMMAlgorithmTask
    FFmmAlgorithm<FTestKernels, FTestParticuleTsm, FTestCellTsm, FTypedLeaf, NbLevels, SizeSubLevels> algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgo(&tree);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Deleting particules ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        particules[idxPart].~FTestParticuleTsm();
    }
    delete [] particules;
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
