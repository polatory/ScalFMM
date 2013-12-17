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

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Kernels/Rotation/FRotationOriginalKernel.hpp"

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
        const int nbParticles = 2;
        const FReal boxWidth = 1.0;
        const FPoint boxCenter(boxWidth/2.0,boxWidth/2.0,boxWidth/2.0);

        struct TestParticle{
            FPoint position;
            FReal forces[3];
            FReal physicalValue;
            FReal potential;
        };

        const int NbLevels      = 3;
        const int SizeSubLevels = 2;

        const int dimGrid = (1 << (NbLevels-1));
        const FReal dimLeaf = (boxWidth/FReal(dimGrid));
        const FReal quarterDimLeaf = (dimLeaf/4.0);

        Print("dimGrid:");
        Print(dimGrid);
        Print("dimLeaf:");
        Print(dimLeaf);
        Print("quarterDimLeaf:");
        Print(quarterDimLeaf);

        TestParticle* const particles = new TestParticle[nbParticles];
        particles[0].position = FPoint(quarterDimLeaf, quarterDimLeaf, quarterDimLeaf);
        particles[0].physicalValue = 0.50;
        particles[1].position = FPoint(2*quarterDimLeaf, quarterDimLeaf, quarterDimLeaf);
        particles[1].physicalValue = -0.10;

        Print("Number of particles:");
        Print(nbParticles);

        for(int idxLeafX = 0 ; idxLeafX < dimGrid ; ++idxLeafX){
            for(int idxLeafY = 0 ; idxLeafY < dimGrid ; ++idxLeafY){
                for(int idxLeafZ = 0 ; idxLeafZ < dimGrid ; ++idxLeafZ){

                particles[1].position = FPoint(FReal(idxLeafX)*dimLeaf + 2*quarterDimLeaf,
                                               FReal(idxLeafY)*dimLeaf + quarterDimLeaf,
                                               FReal(idxLeafZ)*dimLeaf + quarterDimLeaf);

                // Create octree
                OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, boxCenter);
                for(int idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
                    // get copy
                    particles[idxPart].potential = 0.0;
                    particles[idxPart].forces[0] = 0.0;
                    particles[idxPart].forces[1] = 0.0;
                    particles[idxPart].forces[2] = 0.0;

                    std::cout << idxPart << " " << particles[idxPart].position << std::endl;
                }


                // Run FMM
                Print("Fmm...");
                KernelClass kernels(NbLevels,boxWidth, boxCenter);
                FmmClass algo(&tree,&kernels);
                algo.execute();

                // Run direct computation
                Print("Direct...");
                for(int idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
                    for(int idxOther = idxTarget + 1 ; idxOther < nbParticles ; ++idxOther){
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

                            std::cout << indexPartOrig << " x " << particles[indexPartOrig].forces[0]<<" "<<forcesX[idxPart] << std::endl;
                            std::cout << indexPartOrig << " y " << particles[indexPartOrig].forces[1]<<" "<<forcesY[idxPart] << std::endl;
                            std::cout << indexPartOrig << " z " << particles[indexPartOrig].forces[2]<<" "<<forcesZ[idxPart] << std::endl;
                        }
                    });

                    tree.forEachCell([&](CellClass* cell){
                        std::cout << "Multipole:\n";
                        int index_j_k = 0;
                        for (int j = 0 ; j <= P ; ++j ){
                            std::cout <<"[" << j << "] " << "[0] " << cell->getMultipole()[index_j_k].getReal() << " i" << cell->getMultipole()[index_j_k].getImag() << "\t";
                            for (int k=1; k<=j ;++k, ++index_j_k){
                                std::cout << "[" << k << "] " << cell->getMultipole()[index_j_k].getReal() << " i" << cell->getMultipole()[index_j_k].getImag() << "    ";
                            }
                            std::cout << "\n";
                        }
                        std::cout << "\n";
                        std::cout << "Local:\n";
                        index_j_k = 0;
                        for (int j = 0 ; j <= P ; ++j ){
                            std::cout <<"[" << j << "] " << "[0] " << cell->getLocal()[index_j_k].getReal() << " i" << cell->getLocal()[index_j_k].getImag() << "\t";
                            for (int k=1; k<=j ;++k, ++index_j_k){
                                std::cout << "[" << k << "] " << cell->getLocal()[index_j_k].getReal() << " i" << cell->getLocal()[index_j_k].getImag() << "    ";
                            }
                            std::cout << "\n";
                        }
                        std::cout << "\n\n";
                    });
                }

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
            }
        }

        delete[] particles;
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

    static const int P = 5;

    /** Rotation */
    void TestRotation(){
        typedef FRotationCell<P>              CellClass;
        typedef FP2PParticleContainerIndexed  ContainerClass;

        typedef FRotationOriginalKernel<CellClass, ContainerClass, P >          KernelClass;

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



