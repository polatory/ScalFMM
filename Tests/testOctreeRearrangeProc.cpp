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

#include "../Src/Utils/FParameters.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"

#include "../Src/Arranger/FOctreeArrangerProc.hpp"
#include "../Src/Files/FMpiTreeBuilder.hpp"


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FTestParticle               ParticleClass;
    typedef FTestCell                   CellClass;
    typedef FVector<ParticleClass>      ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app(argc, argv);


    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const int NbPart            = FParameters::getValue(argc,argv,"-nb", 20000);
    const FReal FRandMax        = FReal(RAND_MAX);


    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal BoxWidth = 1.0;
    const FReal BoxCenter = 0.5;

    OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, F3DPosition(BoxCenter,BoxCenter,BoxCenter));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    counter.tic();


    {
        FTestParticle particles[NbPart];
        for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particles[idxPart].setPosition(
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/FReal(2.0))),
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/FReal(2.0))),
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/FReal(2.0))));
        }

        FMpiTreeBuilder<ParticleClass>::ArrayToTree(app.global(), particles, NbPart, F3DPosition(BoxCenter,BoxCenter,BoxCenter), BoxWidth, tree);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        do{
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            while( iter.hasNotFinished() ){
                if(FReal(rand())/FRandMax > 0.5){
                    iter.data().setPosition(
                                (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                                (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                                (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)));
                }
                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
    }

    counter.tac();
    std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Arrange ..." << std::endl;
    counter.tic();

    FOctreeArrangerProc<OctreeClass, ContainerClass, ParticleClass> arrange(&tree);
    arrange.rearrange(app.global());

    counter.tac();
    std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Test ..." << std::endl;
    counter.tic();

    { // Check that each particle has been put into the right leaf
        long counterPart = 0;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        MortonIndex interval[2] = {0,0};
        interval[0] = octreeIterator.getCurrentGlobalIndex();

        do{
            const MortonIndex leafIndex = octreeIterator.getCurrentGlobalIndex();
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            while( iter.hasNotFinished() ){
                const MortonIndex particleIndex = tree.getMortonFromPosition( iter.data().getPosition() );
                if( leafIndex != particleIndex){
                    std::cout << "Index problem, should be " << leafIndex <<
                                 " particleIndex "<< particleIndex << std::endl;
                }

                iter.gotoNext();
            }
            counterPart += octreeIterator.getCurrentListTargets()->getSize();
            if(octreeIterator.getCurrentListTargets()->getSize() == 0){
                std::cout << "Problem, leaf is empty at index " << leafIndex << std::endl;
            }
        } while(octreeIterator.moveRight());
        interval[1] = octreeIterator.getCurrentGlobalIndex();

        counterPart = app.global().reduceSum(counterPart);
        if(app.global().processId() == 0 && counterPart != NbPart  * app.global().processCount() ){
            std::cout <<"Wrong particles number, should be " << (NbPart  * app.global().processCount()) << " but is " << counterPart << std::endl;
        }

        MortonIndex*const allintervals = new MortonIndex[ 2 * app.global().processCount() ];
        MPI_Allgather( interval, sizeof(MortonIndex) * 2, MPI_BYTE, allintervals, sizeof(MortonIndex) * 2, MPI_BYTE, MPI_COMM_WORLD);
        if(app.global().processId() == 0){
            for(int idxProc = 1 ; idxProc < app.global().processCount() ; ++idxProc){
                if( allintervals[idxProc*2-1] > allintervals[idxProc*2] ){
                    std::cout << "Interval problem for [" << idxProc-1 << "].max = " << allintervals[idxProc*2-1] <<
                                 "  [" << idxProc << "].min = "<< allintervals[idxProc*2] << std::endl;
                }
            }
        }
        delete[] allintervals;
    }

    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator(&tree);
        OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = NbLevels - 1;
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            // for each cells
            do{
                int countChild = 0;
                CellClass** const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                    if( child[idxChild] ){
                        countChild += 1;
                    }
                }

                if(countChild == 0){
                    std::cout << "Problem at level " << idxLevel << " cell has no child " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Test = " << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



