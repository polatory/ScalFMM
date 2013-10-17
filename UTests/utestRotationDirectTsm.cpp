
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

#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../Src/Components/FTypedLeaf.hpp"
#include "../Src/Extensions/FExtendCellType.hpp"
#include "../Src/Kernels/Rotation/FRotationKernel.hpp"

#include "../Src/Files/FRandomLoader.hpp"

#include "../Src/Core/FFmmAlgorithmThreadTsm.hpp"
#include "../Src/Core/FFmmAlgorithmTsm.hpp"

#include "FUTester.hpp"


/** the test class the rotation and target source model.
  *
  */
class TestRotationDirectTsm : public FUTester<TestRotationDirectTsm> {
    /** The test method to factorize all the test based on different kernels */
    template <class CellClass, class ContainerClass, class KernelClass, class LeafClass,
              class OctreeClass, class FmmClass>
    void RunTest(){
        // Warning in make test the exec dir it Build/UTests
        // Load particles
        const int nbSources = 5000;
        const int nbTargets = 5000;

        FRandomLoader loader(nbSources + nbTargets);

        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels      = 4;
        const int SizeSubLevels = 3;


        struct TestParticle{
            FPoint position;
            FReal forces[3];
            FReal physicalValue;
            FReal potential;
        };


        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

        const FReal physicalValue = 0.10;

        TestParticle* const particlesTargets = new TestParticle[nbTargets];
        for(int idxPart = 0 ; idxPart < nbTargets ; ++idxPart){
            FPoint position;
            loader.fillParticle(&position);
            // put in tree
            tree.insert(position, FParticleTypeTarget, idxPart, physicalValue);
            // get copy
            particlesTargets[idxPart].position = position;
            particlesTargets[idxPart].potential = 0.0;
            particlesTargets[idxPart].physicalValue = physicalValue;
            particlesTargets[idxPart].forces[0] = 0.0;
            particlesTargets[idxPart].forces[1] = 0.0;
            particlesTargets[idxPart].forces[2] = 0.0;
        }

        TestParticle* const particlesSources = new TestParticle[nbSources];
        for(int idxPart = 0 ; idxPart < nbSources ; ++idxPart){
            FPoint position;
            loader.fillParticle(&position);
            // put in tree
            tree.insert(position, FParticleTypeSource, idxPart, physicalValue);
            // get copy
            particlesSources[idxPart].position = position;
            particlesSources[idxPart].physicalValue = physicalValue;
        }


        // Run FMM
        Print("Fmm...");
        //KernelClass kernels(NbLevels,loader.getBoxWidth());
        KernelClass kernels(NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&tree,&kernels);
        algo.execute();

        // Run direct computation
        Print("Direct...");
        for(int idxTarget = 0 ; idxTarget < nbTargets ; ++idxTarget){
            for(int idxOther = 0 ; idxOther < nbSources ; ++idxOther){
                FP2P::NonMutualParticles(
                            particlesSources[idxOther].position.getX(), particlesSources[idxOther].position.getY(),
                            particlesSources[idxOther].position.getZ(),particlesSources[idxOther].physicalValue,
                              particlesTargets[idxTarget].position.getX(), particlesTargets[idxTarget].position.getY(),
                              particlesTargets[idxTarget].position.getZ(),particlesTargets[idxTarget].physicalValue,
                              &particlesTargets[idxTarget].forces[0],&particlesTargets[idxTarget].forces[1],
                              &particlesTargets[idxTarget].forces[2],&particlesTargets[idxTarget].potential);
            }
        }

        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        { // Check that each particle has been summed with all other

            tree.forEachLeaf([&](LeafClass* leaf){
                if( leaf->getTargets()->getNbParticles() ){
                    const FReal*const potentials = leaf->getTargets()->getPotentials();
                    const FReal*const forcesX = leaf->getTargets()->getForcesX();
                    const FReal*const forcesY = leaf->getTargets()->getForcesY();
                    const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                    const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                    const FVector<int>& indexes = leaf->getTargets()->getIndexes();

                    for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                        const int indexPartOrig = indexes[idxPart];
                        potentialDiff.add(particlesTargets[indexPartOrig].potential,potentials[idxPart]);
                        fx.add(particlesTargets[indexPartOrig].forces[0],forcesX[idxPart]);
                        fy.add(particlesTargets[indexPartOrig].forces[1],forcesY[idxPart]);
                        fz.add(particlesTargets[indexPartOrig].forces[2],forcesZ[idxPart]);
                    }
                }
            });
        }

        delete[] particlesTargets;
        delete[] particlesSources;

        // Print for information
        Print("Potential diff is = ");
        Print(potentialDiff.getL2Norm());
        Print(potentialDiff.getInfNorm());
        Print("Fx diff is = ");
        Print(fx.getL2Norm());
        Print(fx.getInfNorm());
        Print("Fy diff is = ");
        Print(fy.getL2Norm());
        Print(fy.getInfNorm());
        Print("Fz diff is = ");
        Print(fz.getL2Norm());
        Print(fz.getInfNorm());

        // Assert
        const FReal MaximumDiff = FReal(0.0001);
        uassert(potentialDiff.getL2Norm() < MaximumDiff);
        uassert(potentialDiff.getInfNorm() < MaximumDiff);
        uassert(fx.getL2Norm()  < MaximumDiff);
        uassert(fx.getInfNorm() < MaximumDiff);
        uassert(fy.getL2Norm()  < MaximumDiff);
        uassert(fy.getInfNorm() < MaximumDiff);
        uassert(fz.getL2Norm()  < MaximumDiff);
        uassert(fz.getInfNorm() < MaximumDiff);
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

    template <const int P>
    class CustomTypedRotationCell : public FRotationCell<P>, public FExtendCellType{
    };

    /** Rotation */
    void TestRotation(){
        typedef CustomTypedRotationCell<P>    CellClass;
        typedef FP2PParticleContainerIndexed  ContainerClass;

        typedef FRotationKernel<CellClass, ContainerClass, P >          KernelClass;

        typedef FTypedLeaf<ContainerClass >                     LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<CellClass, ContainerClass, KernelClass, LeafClass, OctreeClass, FmmClass>();
    }

    void TestRotationThread(){
        typedef CustomTypedRotationCell<P>    CellClass;
        typedef FP2PParticleContainerIndexed  ContainerClass;

        typedef FRotationKernel<CellClass, ContainerClass, P >          KernelClass;

        typedef FTypedLeaf<ContainerClass >                     LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<CellClass, ContainerClass, KernelClass, LeafClass, OctreeClass, FmmClass>();
    }

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestRotationDirectTsm::TestRotation,"Test Rotation Kernel TSM");
        AddTest(&TestRotationDirectTsm::TestRotationThread,"Test Rotation Kernel TSM thread");
    }
};


// You must do this
TestClass(TestRotationDirectTsm)



