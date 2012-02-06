// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
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
            const long NbSmallBoxesPerSide = (1 << (idxHeight-1));
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
                        assert( currentIndex == octreeIterator.getCurrentGlobalIndex());
                        ++currentIndex;
                    } while(octreeIterator.moveRight());
                    // Then number of cells must be divided by 8 at each level
                    assert(currentIndex == nbCell);
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



