// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Core/FFmaParticule.hpp"
#include "../Sources/Core/FExtendForces.hpp"
#include "../Sources/Core/FExtendPotential.hpp"

#include "../Sources/Core/FBasicCell.hpp"
#include "../Sources/Fmb/FExtendFmbCell.hpp"

#include "../Sources/Core/FFMMAlgorithm.hpp"
#include "../Sources/Core/FFMMAlgorithmArray.hpp"
#include "../Sources/Core/FFMMAlgorithmThreaded.hpp"
#include "../Sources/Core/FFMMAlgorithmTask.hpp"

#include "../Sources/Fmb/FFmbKernelsPotentialForces.hpp"
#include "../Sources/Fmb/FFmbKernelsForces.hpp"
#include "../Sources/Fmb/FFmbKernelsPotential.hpp"

#include "../Sources/Files/FFMALoader.hpp"

// With openmp : g++ testFmbAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp ../Sources/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmbAlgorithm.exe
// icpc -openmp -openmp-lib=compat testFmbAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -O2 -o testFmbAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particules is little or longer
  * related that each other
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FFmaParticule}
  */
class FmbParticule : public FFmaParticule, public FExtendForces, public FExtendPotential {
public:
};

/** Custom cell
  *
  */
class FmbCell : public FBasicCell, public FExtendFmbCell {
public:
};


// Simply create particules and try the kernels
int main(int , char ** ){
        const int NbLevels = 9;//10;
        const int SizeSubLevels = 3;//3
        FTic counter;
        const char* const filename = "../../Data/testLoaderFMA.fma"; //"testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma

        FFMALoader<FmbParticule> loader(filename);
        if(!loader.isValide()){
            std::cout << "Loader Error, " << filename << " is missing\n";
            return 1;
        }

        // -----------------------------------------------------

        FOctree<FmbParticule, FmbCell, NbLevels, SizeSubLevels> tree(loader.getBoxWidth(),loader.getCenterOfBox());

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
        FFMMAlgorithmTask<FFmbKernelsPotentialForces, FmbParticule, FmbCell, NbLevels, SizeSubLevels> algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //std::cout << "Foces Sum  x = " << kernels.getForcesSum().getX() << " y = " << kernels.getForcesSum().getY() << " z = " << kernels.getForcesSum().getZ() << std::endl;
        //std::cout << "Potential = " << kernels.getPotential() << std::endl;

        { // get sum forces&potential
            FReal potential = 0;
            F3DPosition forces;
            FOctree<FmbParticule, FmbCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                FList<FmbParticule*>::BasicIterator iter(*octreeIterator.getCurrentList());
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
