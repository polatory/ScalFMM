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
#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Spherical/FSphericalCell.hpp"
#include "../Src/Spherical/FSphericalKernel.hpp"
#include "../Src/Spherical/FSphericalParticle.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"


/**
  * This test compare a previous FMM result with a previous simulation result.
  */


typedef FSphericalParticle       ParticleClass;
typedef FSphericalCell           CellClass;
typedef FVector<ParticleClass>  ContainerClass;

typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

/** The test class */
class TestSphericalWithPrevious : public FUTester<TestSphericalWithPrevious> {
    /** the test */
    void TestTree(){
        // Warning in make test the exec dir it Build/UTests
        const char* const DataFile = "../../Data/utestSpherical.data";
        const char* const ParticleFile = "../../Data/utestSpherical.bin.fma";

        const int NbLevels      = 5;
        const int SizeSubLevels = 3;
        const int DevP = 9;

        // Load the particles file
        FFmaBinLoader<ParticleClass> loader(ParticleFile);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            assert(false);
            return;
        }

        // Create octree
        FSphericalCell::Init(DevP);
        OctreeClass testTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        loader.fillTree(testTree);

        // Run simulation
        KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&testTree,&kernels);
        algo.execute();

        // If needed save the result
        // FTreeIO::Save<OctreeClass, CellClass, ParticleClass, ContainerClass >(DataFile, testTree);

        // Load previous result
        OctreeClass goodTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FTreeIO::Load<OctreeClass, CellClass, ParticleClass, ContainerClass >(DataFile, goodTree);

        // Compare the two simulations
        Print("Check the particles...");
        { // Check that each particle has been summed with all other
            OctreeClass::Iterator testOctreeIterator(&testTree);
            OctreeClass::Iterator goodOctreeIterator(&goodTree);

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
                    assert( FMath::LookEqual(goodIter.data().getPotential(), testIter.data().getPotential()) );
                    assert( FMath::LookEqual(goodIter.data().getPosition().getX(), testIter.data().getPosition().getX()) );
                    assert( FMath::LookEqual(goodIter.data().getPosition().getY(), testIter.data().getPosition().getY()) );
                    assert( FMath::LookEqual(goodIter.data().getPosition().getZ(), testIter.data().getPosition().getZ()) );

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
            OctreeClass::Iterator testOctreeIterator(&testTree);
            OctreeClass::Iterator goodOctreeIterator(&goodTree);

            testOctreeIterator.gotoBottomLeft();
            goodOctreeIterator.gotoBottomLeft();

            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
                do{
                    if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
                        assert(false);
                        break;
                    }

                    for(int idxLocal = 0 ; idxLocal < CellClass::GetLocalSize() ; ++idxLocal){
                        FMath::LookEqual(testOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getReal(),
                                         goodOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getReal());
                        FMath::LookEqual(testOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getImag(),
                                         goodOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getImag());
                    }

                    for(int idxPole = 0 ; idxPole < CellClass::GetPoleSize() ; ++idxPole){
                        FMath::LookEqual(testOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getReal(),
                                         goodOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getReal());
                        FMath::LookEqual(testOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getImag(),
                                         goodOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getImag());
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
        AddTest(&TestSphericalWithPrevious::TestTree,"Test Simu and compare tree");
    }
};



// You must do this
TestClass(TestSphericalWithPrevious)



