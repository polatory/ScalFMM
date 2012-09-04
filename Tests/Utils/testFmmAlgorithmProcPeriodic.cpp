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

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestParticle.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithmThreadProcPeriodic.hpp"
#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"

#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Components/FAbstractSendable.hpp"


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

    typedef FFmmAlgorithmThreadProcPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassSeq;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbParticles      = FParameters::getValue(argc,argv,"-nb", 5);
    const int PeriodicDeep      = FParameters::getValue(argc,argv,"-per", 2);
    // choose in +x dir or -/+x dir or all dirs
    int PeriodicDirs          = (FParameters::existParameter(argc,argv,"-x")?DirMinusX:0) |
                                (FParameters::existParameter(argc,argv,"+x")?DirPlusX:0) |
                                (FParameters::existParameter(argc,argv,"-y")?DirMinusY:0) |
                                (FParameters::existParameter(argc,argv,"+y")?DirPlusY:0) |
                                (FParameters::existParameter(argc,argv,"-z")?DirMinusZ:0) |
                                (FParameters::existParameter(argc,argv,"+z")?DirPlusZ:0);
    if( PeriodicDirs == 0 ) PeriodicDirs =  AllDirs;
    if( FParameters::existParameter(argc,argv,"-nodir") ) PeriodicDirs = 0;

    FMpi app(argc, argv);

    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbParticles << " particles per boxes ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    FRandomLoader<ParticleClass> loader(NbParticles,FReal(1.0),FPoint(0,0,0), app.global().processId());
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    {
        FTestParticle*const particles = new FTestParticle[NbParticles];

        for(int idx = 0 ; idx < NbParticles ; ++idx){
            loader.fillParticle(particles[idx]);
        }

        FMpiTreeBuilder<ParticleClass>::ArrayToTree(app.global(), particles, NbParticles, loader.getCenterOfBox(),
                                                    loader.getBoxWidth(), tree);
        delete[] particles;
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;
    FmmClass algo( app.global(), &tree, PeriodicDeep, PeriodicDirs);
    algo.setKernel(&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    {
        const FTreeCoordinate repetitions = algo.repetitions();
        const int totalRepeatedBox = repetitions.getX() * repetitions.getY() * repetitions.getZ();
        std::cout << "The box is repeated " << repetitions.getX() <<" "<< repetitions.getY()<<" "<<
                     repetitions.getZ() << " there are " << totalRepeatedBox << " boxes in total\n";
        const long long NbParticlesEntireSystem = (NbParticles * app.global().processCount()) * totalRepeatedBox;
        std::cout << "The total number of particles is "  << NbParticlesEntireSystem << "\n";
        FTreeCoordinate min, max;
        algo.repetitionsIntervals(&min, &max);
        std::cout << "Min is " << min << " Max is " << max << std::endl;

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

    {
        OctreeClass treeSeq(NbLevels, SizeSubLevels, FReal(1.0), FPoint(0,0,0));
        for(int idx = 0 ; idx < app.global().processCount() ; ++idx ){
            FRandomLoader<ParticleClass> loaderSeq(NbParticles,FReal(1.0),FPoint(0,0,0), idx);
            loaderSeq.fillTree(treeSeq);
        }

        FmmClassSeq algoSeq( &treeSeq, PeriodicDeep, PeriodicDirs);
        algoSeq.setKernel(&kernels);
        algoSeq.execute();

        { // Ceck if there is number of NbPart summed at level 1
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator octreeIteratorSeq(&treeSeq);
            octreeIteratorSeq.gotoBottomLeft();

            for(int idxLevel = tree.getHeight() - 1 ; idxLevel >= 1 ; --idxLevel ){
                std::cout << "Process level " << idxLevel << std::endl;

                while( octreeIterator.getCurrentGlobalIndex() != octreeIteratorSeq.getCurrentGlobalIndex() ){
                    octreeIteratorSeq.moveRight();
                }

                do{
                    if( octreeIterator.getCurrentGlobalIndex() != octreeIteratorSeq.getCurrentGlobalIndex()){
                        std::cout << "Index problem !!!!!" << std::endl;
                    }

                    if( algo.getWorkingInterval(idxLevel).min <= octreeIteratorSeq.getCurrentGlobalIndex()){
                        if( octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorSeq.getCurrentCell()->getDataUp() ){
                            std::cout << "Up problem at " << octreeIterator.getCurrentGlobalIndex() <<
                                         " Good is " << octreeIteratorSeq.getCurrentCell()->getDataUp() <<
                                         " Bad is " << octreeIterator.getCurrentCell()->getDataUp() << std::endl;
                        }
                        if( octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorSeq.getCurrentCell()->getDataDown() ){
                            std::cout << "Down problem at " << octreeIterator.getCurrentGlobalIndex() <<
                                         " Good is " << octreeIteratorSeq.getCurrentCell()->getDataDown() <<
                                         " Bad is " << octreeIterator.getCurrentCell()->getDataDown() << std::endl;
                        }
                    }
                } while(octreeIterator.moveRight() && octreeIteratorSeq.moveRight());

                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                octreeIteratorSeq.moveUp();
                octreeIteratorSeq.gotoLeft();
            }
        }
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator octreeIteratorSeq(&treeSeq);
            octreeIteratorSeq.gotoBottomLeft();


            while( octreeIterator.getCurrentGlobalIndex() != octreeIteratorSeq.getCurrentGlobalIndex() ){
                octreeIteratorSeq.moveRight();
            }

            do{
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                typename ContainerClass::BasicIterator iterSeq(*octreeIteratorSeq.getCurrentListTargets());

                while( iter.hasNotFinished() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( iter.data().getDataDown() != iterSeq.data().getDataDown() ){
                        std::cout << "Problem L2P + P2P at " << octreeIterator.getCurrentGlobalIndex() <<
                                     " Good is " << iterSeq.data().getDataDown() <<
                                     " Bad is " << iter.data().getDataDown() << std::endl;
                    }

                    iter.gotoNext();
                    iterSeq.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }
    }
    std::cout << "Test is over...\n";

    return 0;
}



