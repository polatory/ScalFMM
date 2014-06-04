// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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

#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Kernels/Rotation/FRotationKernel.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"


/** the test class
  *
  */
class TestRotationDirect : public FUTester<TestRotationDirect> {
    /** The test method to factorize all the test based on different kernels */
    template <class CellClass, class ContainerClass, class KernelClass, class LeafClass,
              class OctreeClass, class FmmClass>
    void RunTest(){
        // Warning in make test the exec dir it Build/UTests
        // Load particles
        const char* const filename = (sizeof(FReal) == sizeof(float))?
                                        "../Data/UTest/Direct.bin.fma.single":
                                        "../Data/UTest/Direct.bin.fma.double";
        FFmaBinLoader loader(filename);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            uassert(false);
            return;
        }
        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels      = 4;
        const int SizeSubLevels = 2;


        struct TestParticle{
            FPoint position;
            FReal forces[3];
            FReal physicalValue;
            FReal potential;
        };

        TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint position;
            FReal physicalValue;
            loader.fillParticle(&position,&physicalValue);
            // put in tree
            tree.insert(position, idxPart, physicalValue);
            // get copy
            particles[idxPart].position = position;
            particles[idxPart].physicalValue = physicalValue;
            particles[idxPart].potential = 0.0;
            particles[idxPart].forces[0] = 0.0;
            particles[idxPart].forces[1] = 0.0;
            particles[idxPart].forces[2] = 0.0;
        }


        // Run FMM
        Print("Fmm...");
        //KernelClass kernels(NbLevels,loader.getBoxWidth());
        KernelClass kernels(NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&tree,&kernels);
        algo.execute();

        // Run direct computation
        Print("Direct...");
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                      particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
                                      &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
                                      &particles[idxTarget].forces[2],&particles[idxTarget].potential,
                                particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
                                &particles[idxOther].forces[0],&particles[idxOther].forces[1],
                                &particles[idxOther].forces[2],&particles[idxOther].potential);
            }
        }

        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        { // Check that each particle has been summed with all other

            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                const FVector<int>& indexes = leaf->getTargets()->getIndexes();

                for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const int indexPartOrig = indexes[idxPart];
                    potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                    fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                    fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                    fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
                }
            });
        }

        delete[] particles;

        // Print for information
        Print("Potential diff is = ");
        Print(potentialDiff.getRelativeL2Norm());
        Print(potentialDiff.getRelativeInfNorm());
        Print("Fx diff is = ");
        Print(fx.getRelativeL2Norm());
        Print(fx.getRelativeInfNorm());
        Print("Fy diff is = ");
        Print(fy.getRelativeL2Norm());
        Print(fy.getRelativeInfNorm());
        Print("Fz diff is = ");
        Print(fz.getRelativeL2Norm());
        Print(fz.getRelativeInfNorm());

        // Assert
        const FReal MaximumDiff = FReal(0.0001);
        uassert(potentialDiff.getRelativeL2Norm() < MaximumDiff);
        uassert(potentialDiff.getRelativeInfNorm() < MaximumDiff);
        uassert(fx.getRelativeL2Norm()  < MaximumDiff);
        uassert(fx.getRelativeInfNorm() < MaximumDiff);
        uassert(fy.getRelativeL2Norm()  < MaximumDiff);
        uassert(fy.getRelativeInfNorm() < MaximumDiff);
        uassert(fz.getRelativeL2Norm()  < MaximumDiff);
        uassert(fz.getRelativeInfNorm() < MaximumDiff);
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }

    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    static const int P = 9;

    /** Rotation */
    void TestRotation(){
        typedef FRotationCell<P>              CellClass;
        typedef FP2PParticleContainerIndexed<>  ContainerClass;

        typedef FRotationKernel<CellClass, ContainerClass, P >          KernelClass;

        typedef FSimpleLeaf<ContainerClass >                     LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<CellClass, ContainerClass, KernelClass, LeafClass, OctreeClass, FmmClass>();
    }

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestRotationDirect::TestRotation,"Test Rotation Kernel");
    }
};


// You must do this
TestClass(TestRotationDirect)



