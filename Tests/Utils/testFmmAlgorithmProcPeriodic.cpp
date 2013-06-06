// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
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

#include "../../Src/Components/FTestParticleContainer.hpp"
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
    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer      ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    typedef FFmmAlgorithmThreadProcPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    typedef FFmmAlgorithmPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassSeq;
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

    FRandomLoader loader(NbParticles,FReal(1.0),FPoint(0,0,0), app.global().processId());
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    {
        struct TestParticle{
            FPoint position;
            const FPoint& getPosition(){
                return position;
            }
        };

        TestParticle*const particles = new TestParticle[NbParticles];

        for(int idx = 0 ; idx < NbParticles ; ++idx){
            loader.fillParticle(&particles[idx].position);
        }

        FVector<TestParticle> finalParticles;
        FMpiTreeBuilder<TestParticle>::ArrayToTree(app.global(), particles, NbParticles, loader.getCenterOfBox(),
                                                    loader.getBoxWidth(), tree.getHeight(), &finalParticles);

        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            tree.insert(finalParticles[idx].position);
        }
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
            ContainerClass* container = (octreeIterator.getCurrentListTargets());
            const long long int*const dataDown = container->getDataDown();

            for(int idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                if( NbParticlesEntireSystem - 1 != dataDown[idxPart]){
                    std::cout << "P2P probleme, should be " << NbParticlesEntireSystem - 1 <<
                                 " iter.data().getDataDown() "<< dataDown[idxPart] << std::endl;
                }
            }
        } while(octreeIterator.moveRight());
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    {
        OctreeClass treeSeq(NbLevels, SizeSubLevels, FReal(1.0), FPoint(0,0,0));
        for(int idx = 0 ; idx < app.global().processCount() ; ++idx ){
            FPoint position;
            FRandomLoader loaderSeq(NbParticles,FReal(1.0),FPoint(0,0,0), idx);
            for(int idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
                loaderSeq.fillParticle(&position);
                treeSeq.insert(position);
            }
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
                ContainerClass* container = (octreeIterator.getCurrentListTargets());
                const long long int*const dataDown = container->getDataDown();

                ContainerClass* containerValide = (octreeIteratorSeq.getCurrentListTargets());
                const long long int*const dataDownValide = containerValide->getDataDown();

                for(int idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( dataDown[idxPart] != dataDownValide[idxPart]){
                        std::cout << "Problem on leaf " << octreeIterator.getCurrentGlobalIndex() <<
                                     " part " << idxPart << " valide data down " << dataDownValide[idxPart] <<
                                     " invalide " << dataDown[idxPart] << "\n";
                        std::cout << "Data down for leaf is: valide " << octreeIteratorSeq.getCurrentCell()->getDataDown()
                                  << " invalide " << octreeIterator.getCurrentCell()->getDataDown()
                                  << " size is: valide " <<  octreeIteratorSeq.getCurrentListTargets()->getNbParticles()
                                  << " invalide " << octreeIterator.getCurrentListTargets()->getNbParticles() << std::endl;
                    }
                }
            } while(octreeIterator.moveRight());
        }
    }
    std::cout << "Test is over...\n";

    return 0;
}



