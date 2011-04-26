// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FList.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Extenssions/FExtendForces.hpp"
#include "../Src/Extenssions/FExtendPotential.hpp"

#include "../Src/Extenssions/FExtendParticleType.hpp"
#include "../Src/Extenssions/FExtendCellType.hpp"


#include "../Src/Components/FBasicCell.hpp"
#include "../Src/Fmb/FExtendFmbCell.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmTsm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithmThreadTsm.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Components/FTypedLeaf.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"


// With openmp : g++ testFmbTsmAlgorithm.cpp ../Src/Utils/FAssertable.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmbTsmAlgorithm.exe
// icpc -openmp -openmp-lib=compat testFmbTsmAlgorithm.cpp ../Src/Utils/FAssertable.cpp ../Src/Utils/FDebug.cpp -O2 -o testFmbTsmAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FFmaParticle}
  */
class FmbParticle : public FFmaParticle, public FExtendForces, public FExtendPotential {
public:
};

class FmbParticleTyped : public FmbParticle, public FExtendParticleType {
public:
};

/** Custom cell
  *
  */
class FmbCell : public FBasicCell, public FExtendFmbCell {
public:
};

class FmbCellTyped : public FmbCell, public FExtendCellType {
public:
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Fmb on a Tsm system.\n";
    std::cout << ">> It compares the results between Tms and no Tms (except P2P & L2P).\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 9;//10;
    const int SizeSubLevels = 3;//3
    FTic counter;
    const long NbPart = 200000;//2000000
    const double BoxWidth = 1.0;
    const F3DPosition CenterOfBox(0.5,0.5,0.5);

    // -----------------------------------------------------

    std::cout << "Creating " << NbPart << " particles ..." << std::endl;
    counter.tic();

    FmbParticle* particles = new FmbParticle[NbPart];
    FmbParticleTyped* particlesTyped = new FmbParticleTyped[NbPart*2];

    for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        const double x = FReal(rand())/RAND_MAX;
        const double y = FReal(rand())/RAND_MAX;
        const double z = FReal(rand())/RAND_MAX;
        // Particle for standart model
        particles[idxPart].setPosition(x,y,z);
        particles[idxPart].setPhysicalValue(1);

        // Create a clone for typed (Tsm) version
        particlesTyped[idxPart*2].setPosition(x,y,z);
        particlesTyped[idxPart*2+1].setPosition(x,y,z);

        particlesTyped[idxPart*2].setPhysicalValue(1);
        particlesTyped[idxPart*2+1].setPhysicalValue(1);

        particlesTyped[idxPart*2].setAsSource();
        particlesTyped[idxPart*2+1].setAsTarget();
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    FOctree<FmbParticle, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels> tree(BoxWidth,CenterOfBox);
    FOctree<FmbParticleTyped, FmbCellTyped, FTypedLeaf, NbLevels, SizeSubLevels> treeTyped(BoxWidth,CenterOfBox);


    std::cout << "Inserting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        tree.insert(&particles[idxPart]);
        treeTyped.insert(&particlesTyped[idxPart*2]);
        treeTyped.insert(&particlesTyped[idxPart*2+1]);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    FFmbKernels<FmbParticle, FmbCell, NbLevels> kernels(BoxWidth);
    FFmbKernels<FmbParticleTyped, FmbCellTyped, NbLevels> kernelsTyped(BoxWidth);
    //FFmmAlgorithm FFmmAlgorithmThread
    FFmmAlgorithmThread<FFmbKernels, FmbParticle, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels> algo(&tree,&kernels);
    FFmmAlgorithmThreadTsm<FFmbKernels, FmbParticleTyped, FmbCellTyped, FTypedLeaf, NbLevels, SizeSubLevels> algoTyped(&treeTyped,&kernelsTyped);
    algo.execute();
    algoTyped.execute();

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    // Here we compare the cells of the trees that must contains the same values

    std::cout << "Start checking ..." << std::endl;
    {
        FOctree<FmbParticle, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        FOctree<FmbParticleTyped, FmbCellTyped, FTypedLeaf, NbLevels, SizeSubLevels>::Iterator octreeIteratorTyped(&treeTyped);
        octreeIteratorTyped.gotoBottomLeft();

        for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
            std::cout << "\t test level " << idxLevel << "\n";

            do{
                bool poleDiff = false;
                bool localDiff = false;
                for(int idxValues = 0 ; idxValues < FExtendFmbCell::MultipoleSize && !(poleDiff && localDiff); ++idxValues){
                    const FComplexe pole = octreeIterator.getCurrentCell()->getMultipole()[idxValues];
                    const FComplexe poleTyped = octreeIteratorTyped.getCurrentCell()->getMultipole()[idxValues];
                    if(!FMath::LookEqual(pole.getImag(),poleTyped.getImag()) || !FMath::LookEqual(pole.getReal(),poleTyped.getReal())){
                        poleDiff = true;
                        printf("Pole diff imag( %.15e , %.15e ) real( %.15e , %.15e)\n",
                               pole.getImag(),poleTyped.getImag(),pole.getReal(),poleTyped.getReal());
                    }
                    const FComplexe local = octreeIterator.getCurrentCell()->getLocal()[idxValues];
                    const FComplexe localTyped = octreeIteratorTyped.getCurrentCell()->getLocal()[idxValues];
                    if(!FMath::LookEqual(local.getImag(),localTyped.getImag()) || !FMath::LookEqual(local.getReal(),localTyped.getReal())){
                        localDiff = true;
                        printf("Pole diff imag( %.15e , %.15e ) real( %.15e , %.15e)\n",
                               local.getImag(),localTyped.getImag(),local.getReal(),localTyped.getReal());
                    }
                }
                if(poleDiff){
                    std::cout << "Multipole error at level " << idxLevel << "\n";
                }
                if(localDiff){
                    std::cout << "Locale error at level " << idxLevel << "\n";
                }
            } while(octreeIterator.moveRight() && octreeIteratorTyped.moveRight());

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorTyped.moveUp();
            octreeIteratorTyped.gotoLeft();
        }
    }

    std::cout << "Done ..." << std::endl;

    delete [] particles;
    delete [] particlesTyped;

    return 0;
}


// [--LICENSE--]
