// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestParticle.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    {
        //typedef FTestParticle               ParticleClass;
        //typedef FTestCell                   CellClass;
        typedef FSphericalParticle             ParticleClass;
        typedef FSphericalCell                 CellClass;

        typedef FVector<ParticleClass>      ContainerClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
        //typedef FTestKernels<ParticleClass, CellClass, ContainerClass >         KernelClass;
        typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
        ///////////////////////What we do/////////////////////////////
        std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
        //////////////////////////////////////////////////////////////

        const int NbLevels      = FParameters::getValue(argc,argv,"-h", 5);
        const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
        const int NbPart       = FParameters::getValue(argc,argv,"-nb", 2000000);
        const int DevP         = FParameters::getValue(argc,argv,"-p", 5);
        const FReal FRandMax    = FReal(RAND_MAX);
        const FPoint centerOfBox = FPoint(0.5,0.5,0.5);
        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        CellClass::Init(DevP);
        const FReal boxWidth = 1.0;
        OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, centerOfBox);

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
        std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
        counter.tic();

        {
            ParticleClass particleToFill;
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
        //KernelClass kernels;
        KernelClass kernels(DevP, NbLevels,boxWidth, centerOfBox);
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



