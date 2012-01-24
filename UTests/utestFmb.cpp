// [--License--]


#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/FComputeCell.hpp"
#include "../Src/Kernels/FElecForcesKernels.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"
#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"


typedef FmbParticle             ParticleClass;
typedef FComputeCell            CellClass;
typedef FVector<ParticleClass>  ContainerClass;

//typedef FFmbKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;
typedef FElecForcesKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

class TestFmb : public FUTester<TestFmb> {

    void TestTree(){
        const int NbLevels      = 5;
        const int SizeSubLevels = 3;
        const int DevP = 12;
        FComputeCell::Init(DevP);

        FFmaBinLoader<ParticleClass> loader("../Data/utestFmb.bin.fma");
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            return;
        }

        OctreeClass testTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        {
            ParticleClass particleToFill;
            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(particleToFill);
                testTree.insert(particleToFill);
            }
        }


        KernelClass kernels(DevP,NbLevels,loader.getBoxWidth());
        FmmClass algo(&testTree,&kernels);
        algo.execute();

        FTreeIO::Save<OctreeClass, CellClass, ParticleClass, FTreeIO::Copier<CellClass, ParticleClass> >("../UTests/data/fmb.data", testTree);

        OctreeClass goodTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FTreeIO::Load<OctreeClass, CellClass, ParticleClass, FTreeIO::Copier<CellClass, ParticleClass> >("../Data/utestFmb.data", goodTree);

        Print("Check the particles...");

        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator testOctreeIterator(&testTree);
            typename OctreeClass::Iterator goodOctreeIterator(&goodTree);

            testOctreeIterator.gotoBottomLeft();
            goodOctreeIterator.gotoBottomLeft();

            do{
                if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
                    assert(false);
                    break;
                }

                if(testOctreeIterator.getCurrentListSrc()->getSize() != goodOctreeIterator.getCurrentListSrc()->getSize()){
                    assert(false);
                    break;
                }

                typename ContainerClass::BasicIterator goodIter(*goodOctreeIterator.getCurrentListTargets());
                typename ContainerClass::BasicIterator testIter(*testOctreeIterator.getCurrentListTargets());

                while( goodIter.hasNotFinished() ){
                    assert( memcmp(&goodIter.data(), &testIter.data(), sizeof(ParticleClass)) == 0);

                    goodIter.gotoNext();
                    testIter.gotoNext();
                }


                if(!testOctreeIterator.moveRight()){
                    if(goodOctreeIterator.moveRight()){
                        assert(false);
                    }
                    break;
                }
                if(!goodOctreeIterator.moveRight()){
                    assert(false);
                    break;
                }

            } while(true);
        }
        Print("Check the leaves...");
        { // Ceck if there is number of NbPart summed at level 1
            typename OctreeClass::Iterator testOctreeIterator(&testTree);
            typename OctreeClass::Iterator goodOctreeIterator(&goodTree);

            testOctreeIterator.gotoBottomLeft();
            goodOctreeIterator.gotoBottomLeft();

            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
                do{
                    if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
                        assert(false);
                        break;
                    }

                    assert( memcmp(testOctreeIterator.getCurrentCell()->getLocal(),
                                   goodOctreeIterator.getCurrentCell()->getLocal(), CellClass::GetExp() * sizeof(FComplexe)) == 0);

                    assert( memcmp(testOctreeIterator.getCurrentCell()->getMultipole(),
                                   goodOctreeIterator.getCurrentCell()->getMultipole(),CellClass::GetExp() * sizeof(FComplexe)) == 0);

                    if(!testOctreeIterator.moveRight()){
                        if(goodOctreeIterator.moveRight()){
                            assert(false);
                        }
                        break;
                    }
                    if(!goodOctreeIterator.moveRight()){
                        assert(false);
                        break;
                    }

                } while(true);

                testOctreeIterator.moveUp();
                testOctreeIterator.gotoLeft();

                goodOctreeIterator.moveUp();
                goodOctreeIterator.gotoLeft();
            }
        }
        Print("Over...");
    }


    // set test
    void SetTests(){
        AddTest(&TestFmb::TestTree,"Test Simu and compare tree");
    }
};



// You must do this
TestClass(TestFmb)


// [--END--]
