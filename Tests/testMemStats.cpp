// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FParameters.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"


// Compile by : g++ testFmmAlgorithm.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    {
        typedef FTestParticle               ParticleClass;
        typedef FTestCell                   CellClass;
        typedef FVector<ParticleClass>      ContainerClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FTestKernels<ParticleClass, CellClass, ContainerClass >         KernelClass;

        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
        ///////////////////////What we do/////////////////////////////
        std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
        //////////////////////////////////////////////////////////////

        const int NbLevels      = FParameters::getValue(argc,argv,"-h", 9);
        const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
        const long NbPart       = FParameters::getValue(argc,argv,"-nb", 2000000);
        const FReal FRandMax    = FReal(RAND_MAX);

        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        OctreeClass tree(NbLevels, SizeSubLevels, 1.0, F3DPosition(0.5,0.5,0.5));

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
        std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
        counter.tic();

        {
            FTestParticle particleToFill;
            for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
                particleToFill.setPosition(FReal(rand())/FRandMax,FReal(rand())/FRandMax,FReal(rand())/FRandMax);
                tree.insert(particleToFill);
            }
        }

        counter.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Working on particles ..." << std::endl;
        counter.tic();

        // FTestKernels FBasicKernels
        KernelClass kernels;
        FmmClass algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////
    }

    std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
    std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
    std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";

    return 0;
}


// [--LICENSE--]
