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

#include "../Src/Arranger/FOctreeArranger.hpp"


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

    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const int NbPart           = FParameters::getValue(argc,argv,"-nb", 2000000);
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
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();


    {

        FTestParticle particleToFill;
        for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particleToFill.setPosition(
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)));
            tree.insert(particleToFill);
        }
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
                iter.data().setPosition(
                            (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                            (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                            (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)));

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

    FOctreeArranger<OctreeClass, ContainerClass, ParticleClass> arrange(&tree);
    arrange.rearrange();

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
        do{
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            const MortonIndex leafIndex = octreeIterator.getCurrentGlobalIndex();

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

        if( counterPart != NbPart ){
            std::cout <<"Wrong particles number, should be " << NbPart << " but is " << counterPart << std::endl;
        }
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



