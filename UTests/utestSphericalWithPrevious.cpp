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
#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"
#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/**
  * This test compare a previous FMM result with a previous simulation result.
  */


typedef FSphericalCell           CellClass;
typedef FP2PParticleContainerIndexed<>  ContainerClass;

typedef FSphericalKernel< CellClass, ContainerClass >          KernelClass;

typedef FSimpleLeaf< ContainerClass >                     LeafClass;
typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

/** To check if a value is correct */
bool IsSimilar(const FReal good, const FReal other){
    const FReal Epsilon = FReal(0.0001);
    return (FMath::Abs(good-other)/FMath::Abs(good)) < Epsilon;
}

/** The test class */
class TestSphericalWithPrevious : public FUTester<TestSphericalWithPrevious> {
    /** the test */
    void TestTree(){
        // Warning in make test the exec dir it Build/UTests
        const char* const DataFile = (sizeof(FReal) == sizeof(float))?
                    "../../Data/utestSphericalPrevious.data.single":
                    "../../Data/utestSphericalPrevious.data.double";
        const char* const ParticleFile = (sizeof(FReal) == sizeof(float))?
                    "../../Data/utestDirect.bin.fma.single":
                    "../../Data/utestDirect.bin.fma.double";

        const int NbLevels      = 5;
        const int SizeSubLevels = 3;
        const int DevP = 9;

        // Load the particles file
        FFmaBinLoader loader(ParticleFile);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            uassert(false);
            return;
        }

        // Create octree
        FSphericalCell::Init(DevP);
        OctreeClass testTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint position;
            FReal physicalValue = 0.0;
            loader.fillParticle(&position,&physicalValue);
            // put in tree
            testTree.insert(position, idxPart, physicalValue);
        }

        // Run simulation
        KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&testTree,&kernels);
        algo.execute();

        // If needed save the result
        //FTreeIO::Save<OctreeClass, CellClass, LeafClass, ContainerClass >(DataFile, testTree);

        // Load previous result
        OctreeClass goodTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FTreeIO::Load<OctreeClass, CellClass, LeafClass, ContainerClass >(DataFile, goodTree);

        // Compare the two simulations
        Print("Check the particles...");
        { // Check that each particle has been summed with all other
            OctreeClass::Iterator testOctreeIterator(&testTree);
            OctreeClass::Iterator goodOctreeIterator(&goodTree);

            testOctreeIterator.gotoBottomLeft();
            goodOctreeIterator.gotoBottomLeft();

            do{
                if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
                    uassert(false);
                    break;
                }

                if(testOctreeIterator.getCurrentListSrc()->getNbParticles() != goodOctreeIterator.getCurrentListSrc()->getNbParticles()){
                    uassert(false);
                    break;
                }

                const ContainerClass* testLeaf = testOctreeIterator.getCurrentListSrc();
                const ContainerClass* goodLeaf = goodOctreeIterator.getCurrentListSrc();

                for(int idxPart = 0 ; idxPart < testLeaf->getNbParticles() ; ++idxPart ){
                    uassert( IsSimilar(goodLeaf->getPotentials()[idxPart], testLeaf->getPotentials()[idxPart]) );
                    uassert( IsSimilar(goodLeaf->getForcesX()[idxPart], testLeaf->getForcesX()[idxPart]) );
                    uassert( IsSimilar(goodLeaf->getForcesY()[idxPart], testLeaf->getForcesY()[idxPart]) );
                    uassert( IsSimilar(goodLeaf->getForcesZ()[idxPart], testLeaf->getForcesZ()[idxPart]) );
                }

                if(!testOctreeIterator.moveRight()){
                    if(goodOctreeIterator.moveRight()){
                        uassert(false);
                    }
                    break;
                }
                if(!goodOctreeIterator.moveRight()){
                    uassert(false);
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
                        uassert(false);
                        break;
                    }

                    for(int idxLocal = 0 ; idxLocal < CellClass::GetLocalSize() ; ++idxLocal){
                        IsSimilar(testOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getReal(),
                                         goodOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getReal());
                        IsSimilar(testOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getImag(),
                                         goodOctreeIterator.getCurrentCell()->getLocal()[idxLocal].getImag());
                    }

                    for(int idxPole = 0 ; idxPole < CellClass::GetPoleSize() ; ++idxPole){
                        IsSimilar(testOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getReal(),
                                         goodOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getReal());
                        IsSimilar(testOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getImag(),
                                         goodOctreeIterator.getCurrentCell()->getMultipole()[idxPole].getImag());
                    }

                    if(!testOctreeIterator.moveRight()){
                        if(goodOctreeIterator.moveRight()){
                            uassert(false);
                        }
                        break;
                    }
                    if(!goodOctreeIterator.moveRight()){
                        uassert(false);
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



