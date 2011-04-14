// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Components/FFmaParticle.hpp"
#include "../Sources/Extenssions/FExtendForces.hpp"
#include "../Sources/Extenssions/FExtendPotential.hpp"

#include "../Sources/Extenssions/FExtendParticleType.hpp"

#include "../Sources/Components/FBasicCell.hpp"
#include "../Sources/Fmb/FExtendFmbCell.hpp"

#include "../Sources/Core/FFmmAlgorithm.hpp"
#include "../Sources/Core/FFmmAlgorithmArray.hpp"
#include "../Sources/Core/FFmmAlgorithmThreaded.hpp"
#include "../Sources/Core/FFmmAlgorithmTask.hpp"

#include "../Sources/Components/FSimpleLeaf.hpp"

#include "../Sources/Fmb/FFmbKernelsPotentialForces.hpp"
#include "../Sources/Fmb/FFmbKernelsForces.hpp"
#include "../Sources/Fmb/FFmbKernelsPotential.hpp"

#include "../Sources/Files/FFMATsmLoader.hpp"

// With openmp : g++ testFmbTsmAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp ../Sources/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmbTsmAlgorithm.exe
// icpc -openmp -openmp-lib=compat testFmbTsmAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -O2 -o testFmbTsmAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FFmaParticle}
  */
class FmbParticle : public FFmaParticle, public FExtendParticleType, public FExtendForces, public FExtendPotential {
public:
};

/** Custom cell
  *
  */
class FmbCell : public FBasicCell, public FExtendFmbCell {
public:
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Fmb on a Tsm system.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 9;//10;
    const int SizeSubLevels = 3;//3
    FTic counter;
    const char* const defaultFilename = "testLoaderFMA.tor.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
    const char* filename;

    if(argc == 1){
        std::cout << "You have to give a .fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    FFMATsmLoader<FmbParticle> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    FOctree<FmbParticle, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    counter.tic();

    FmbParticle* particles = new FmbParticle[loader.getNumberOfParticles()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart]);
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Inserting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        tree.insert(&particles[idxPart]);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    //FFmbKernelsPotentialForces FFmbKernelsForces FFmbKernelsPotential
    FFmbKernelsPotentialForces<FmbParticle, FmbCell, NbLevels> kernels(loader.getBoxWidth());
    //FFmmAlgorithm FFmmAlgorithmThreaded FFmmAlgorithmArray FFmmAlgorithmTask
    FFmmAlgorithm<FFmbKernelsPotentialForces, FmbParticle, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels> algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //std::cout << "Foces Sum  x = " << kernels.getForcesSum().getX() << " y = " << kernels.getForcesSum().getY() << " z = " << kernels.getForcesSum().getZ() << std::endl;
    //std::cout << "Potential = " << kernels.getPotential() << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
        FOctree<FmbParticle, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            FList<FmbParticle*>::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.isValide() ){
                potential += iter.value()->getPotential() * iter.value()->getPhysicalValue();
                forces += iter.value()->getForces();

                //printf("x = %e y = %e z = %e \n",iter.value()->getPosition().getX(),iter.value()->getPosition().getY(),iter.value()->getPosition().getZ());
                //printf("\t fx = %e fy = %e fz = %e \n",iter.value()->getForces().getX(),iter.value()->getForces().getY(),iter.value()->getForces().getZ());

                //printf("\t\t Sum Forces ( %e , %e , %e)\n",
                //forces.getX(),forces.getY(),forces.getZ());

                iter.progress();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }


    // -----------------------------------------------------

    std::cout << "Deleting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        particles[idxPart].~FmbParticle();
    }
    delete [] particles;
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
