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
#include "FUTester.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FBlockAllocator.hpp"
#include "../Src/Containers/FVector.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/FPoint.hpp"

#include "../Src/Components/FBasicParticleContainer.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Utils/FTic.hpp"

/**
  In this test we create a lot of different octree by using various height and subheigt
  then we insert particle in all leaf and we test that all leaves
  and all cells has been created.
  */


/** this class test the octree container */
class TestOctree : public FUTester<TestOctree> {
    typedef FBasicCell                   CellClass;
    typedef FBasicParticleContainer<0>      ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree<CellClass, ContainerClass , LeafClass , FBasicBlockAllocator<CellClass> >  OctreeClass;

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

                OctreeClass tree(idxHeight, idxSub, BoxWidth, FPoint(BoxCenter,BoxCenter,BoxCenter));

                // fill the tree
                for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
                    for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                        for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                            const FPoint pos(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2,
                                                       FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2,
                                                       FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2);
                            tree.insert(pos);
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



