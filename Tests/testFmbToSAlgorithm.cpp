// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Components/FFmaParticule.hpp"
#include "../Sources/Extenssions/FExtendForces.hpp"
#include "../Sources/Extenssions/FExtendPotential.hpp"

#include "../Sources/Extenssions/FExtendParticuleType.hpp"

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

#include "../Sources/Files/FFMAToSLoader.hpp"

// With openmp : g++ testFmbToSAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp ../Sources/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmbToSAlgorithm.exe
// icpc -openmp -openmp-lib=compat testFmbToSAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -O2 -o testFmbToSAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particules is little or longer
  * related that each other
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FFmaParticule}
  */
class FmbParticule : public FFmaParticule, public FExtendParticuleType, public FExtendForces, public FExtendPotential {
public:
};

/** Custom cell
  *
  */
class FmbCell : public FBasicCell, public FExtendFmbCell {
public:
};


// Simply create particules and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Fmb on a ToS system.\n";
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

    FFMAToSLoader<FmbParticule> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    FOctree<FmbParticule, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating " << loader.getNumberOfParticules() << " particules ..." << std::endl;
    counter.tic();

    FmbParticule* particules = new FmbParticule[loader.getNumberOfParticules()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticules() ; ++idxPart){
        loader.fillParticule(&particules[idxPart]);
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Inserting particules ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < loader.getNumberOfParticules() ; ++idxPart){
        tree.insert(&particules[idxPart]);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particules ..." << std::endl;
    counter.tic();

    //FFmbKernelsPotentialForces FFmbKernelsForces FFmbKernelsPotential
    FFmbKernelsPotentialForces<FmbParticule, FmbCell, NbLevels> kernels(loader.getBoxWidth());
    //FFMMAlgorithm FFMMAlgorithmThreaded FFMMAlgorithmArray FFMMAlgorithmTask
    FFmmAlgorithm<FFmbKernelsPotentialForces, FmbParticule, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels> algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //std::cout << "Foces Sum  x = " << kernels.getForcesSum().getX() << " y = " << kernels.getForcesSum().getY() << " z = " << kernels.getForcesSum().getZ() << std::endl;
    //std::cout << "Potential = " << kernels.getPotential() << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
        FOctree<FmbParticule, FmbCell, FSimpleLeaf, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            FList<FmbParticule*>::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
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

    std::cout << "Deleting particules ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < loader.getNumberOfParticules() ; ++idxPart){
        particules[idxPart].~FmbParticule();
    }
    delete [] particules;
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
