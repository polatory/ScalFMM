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

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestParticle.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FTestParticle               ParticleClass;
    typedef FTestCell                   CellClass;
    typedef FVector<ParticleClass>      ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels<ParticleClass, CellClass, ContainerClass >         KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbParticles      = FParameters::getValue(argc,argv,"-nb", 1000);
    const int PeriodicDeep      = FParameters::getValue(argc,argv,"-per", 2);
    // choose in +x dir or -/+x dir or all dirs
    int PeriodicDirs          = (FParameters::existParameter(argc,argv,"-x")?DirMinusX:0) |
                                (FParameters::existParameter(argc,argv,"+x")?DirPlusX:0) |
                                (FParameters::existParameter(argc,argv,"-y")?DirMinusY:0) |
                                (FParameters::existParameter(argc,argv,"+y")?DirPlusY:0) |
                                (FParameters::existParameter(argc,argv,"-z")?DirMinusZ:0) |
                                (FParameters::existParameter(argc,argv,"+z")?DirPlusZ:0);
    if( PeriodicDirs == 0 ) PeriodicDirs =  AllDirs;

    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbParticles << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    FRandomLoader<ParticleClass> loader(NbParticles);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    loader.fillTree(tree);

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;
    FmmClass algo( &tree, PeriodicDeep, PeriodicDirs);
    algo.setKernel(&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    { // Check that each particle has been summed with all other
        long long int counterNbPart = 0;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            // on each leaf we should have the same number of particles
            if(octreeIterator.getCurrentCell()->getDataUp() != octreeIterator.getCurrentListSrc()->getSize() ){
                    std::cout << "Problem P2M Data up = " << octreeIterator.getCurrentCell()->getDataUp() <<
                                 " Size = " << octreeIterator.getCurrentListSrc()->getSize() << "\n";
            }
            // we also count the number of particles.
            counterNbPart += octreeIterator.getCurrentListSrc()->getSize();
        } while(octreeIterator.moveRight());

        if( counterNbPart != NbParticles){
            std::cout << "Problem global nb part, counter = " << counterNbPart << " created = " << NbParticles << std::endl;
        }
    }
    {
        const FTreeCoordinate repetitions = algo.repetitions();
        const int totalRepeatedBox = repetitions.getX() * repetitions.getY() * repetitions.getZ();
        std::cout << "The box is repeated " << repetitions.getX() <<" "<< repetitions.getY()<<" "<<
                     repetitions.getZ() << " there are " << totalRepeatedBox << " boxes in total\n";
        const long long NbParticlesEntireSystem = NbParticles * totalRepeatedBox;
        std::cout << "The total number of particles is "  << NbParticlesEntireSystem << "\n";

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            while( iter.hasNotFinished() ){
                if( NbParticlesEntireSystem - 1 != iter.data().getDataDown()){
                    std::cout << "P2P probleme, should be " << NbParticlesEntireSystem - 1 <<
                                 " iter.data().getDataDown() "<< iter.data().getDataDown() << std::endl;
                }

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



