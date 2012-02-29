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
#include "FUTester.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FBasicParticle.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Utils/FTic.hpp"

/**
  In this test we create a lot of different octree by using various height and subheigt
  then we insert particle in all leaf and we test that all leaves
  and all cells has been created.
  */


/** this class test the octree container */
class TestOctree : public FUTester<TestOctree> {
    typedef FBasicParticle               ParticleClass;
    typedef FBasicCell                   CellClass;
    typedef FVector<ParticleClass>      ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

    // test size
    void TestAll(){
        const FReal BoxWidth = 1.0;
        const FReal BoxCenter = 0.5;

        // We try for many levels
        for(int idxHeight = 2 ; idxHeight < 8 ; ++idxHeight){

            // Compute the number of leaves for a tree of this height
            const int NbSmallBoxesPerSide = (1 << (idxHeight-1));
            const FReal SmallBoxWidth = BoxWidth / FReal(NbSmallBoxesPerSide);
            const FReal SmallBoxWidthDiv2 = SmallBoxWidth / 2;

            const int NbPart = NbSmallBoxesPerSide * NbSmallBoxesPerSide * NbSmallBoxesPerSide;

            // For each level we try many sub-levels
            for(int idxSub = 1 ; idxSub < idxHeight ; ++idxSub){

                OctreeClass tree(idxHeight, idxSub, BoxWidth, F3DPosition(BoxCenter,BoxCenter,BoxCenter));

                // fill the tree
                ParticleClass particleToFill;
                for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
                    for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                        for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                            particleToFill.setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2,
                                                       FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2,
                                                       FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2);
                            tree.insert(particleToFill);
                        }
                    }
                }

                // test all cells
                OctreeClass::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                int nbCell = NbPart;
                for(int idxLevel = idxHeight - 1 ; idxLevel >= 1 ; --idxLevel ){
                    MortonIndex currentIndex = 0;
                    do{
                        // Morton index must increase one by one
                        uassert( currentIndex == octreeIterator.getCurrentGlobalIndex());
                        ++currentIndex;
                    } while(octreeIterator.moveRight());
                    // Then number of cells must be divided by 8 at each level
                    uassert(currentIndex == nbCell);
                    nbCell /= 8;

                    octreeIterator.moveUp();
                    octreeIterator.gotoLeft();
                }
            }
        }
    }

    // set test
    void SetTests(){
        AddTest(&TestOctree::TestAll,"Test Octree");
    }
};

// You must do this
TestClass(TestOctree)



