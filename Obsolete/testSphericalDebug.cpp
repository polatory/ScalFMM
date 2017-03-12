// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// Keep in private GIT
// @SCALFMM_PRIVATE
// @FUSE_BLAS
#define DEBUG_SPHERICAL_M2L
#define  BLAS_SPHERICAL_COMPRESS
#define  BLAS_M2L_P
#include <iostream>

#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalRotationKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalBlockBlasKernel.hpp"

#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../UTests/FUTester.hpp"

/*
  In this test we compare the spherical fmm results and the direct results.
  */

/** the test class
  *
  */
class TestSphericalDirect : public FUTester<TestSphericalDirect> {
    /** The test method to factorize all the test based on different kernels */
    template < class FReal, class CellClass, class ContainerClass, class KernelClass, class LeafClass,
               class OctreeClass, class FmmClass>
    void RunTest(const bool isBlasKernel){
        const int DevP = 2;
        // Warning in make test the exec dir it Build/UTests
        // Load particles
        const FSize nbParticles = 2;
        const FReal boxWidth = 1.0;
        const FPoint<FReal> boxCenter(boxWidth/2.0,boxWidth/2.0,boxWidth/2.0);

        struct TestParticle{
            FPoint<FReal> position;
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

        FSphericalCell<FReal>::Init(DevP, isBlasKernel);

        TestParticle* const particles = new TestParticle[nbParticles];
        particles[0].position = FPoint<FReal>(quarterDimLeaf, quarterDimLeaf, quarterDimLeaf);
        particles[0].physicalValue = 1;
        //        particles[1].position = FPoint<FReal>(2*quarterDimLeaf, quarterDimLeaf, quarterDimLeaf);
        particles[1].physicalValue = -1;

        Print("Number of particles:");
        Print(nbParticles);
        //int idxLeafY = 0,idxLeafZ = 0 ;
        std::cout << "\n ------  Loop starts ---"<< std::endl ;
        for(int idxLeafX = 2 ; idxLeafX < dimGrid ; ++idxLeafX){
            /*for(int idxLeafY = 0 ; idxLeafY < dimGrid ; ++idxLeafY)*/{
                /*  for(int idxLeafZ = 0 ; idxLeafZ < dimGrid ; ++idxLeafZ)*/{
                    //std::cout << "Shift : " << idxLeafX << " " << idxLeafY << " " << idxLeafZ << std::endl;

                    /*  particles[1].position = FPoint<FReal>(FReal(idxLeafX)*dimLeaf + 3*quarterDimLeaf,
                                                   FReal(idxLeafY)*dimLeaf + quarterDimLeaf,
                                                   FReal(idxLeafZ)*dimLeaf + quarterDimLeaf);*/
                    //                   particles[1].position = FPoint<FReal>(FReal(idxLeafX)*dimLeaf + quarterDimLeaf,
                    //                               quarterDimLeaf,
                    //                               quarterDimLeaf);
                    particles[1].position = FPoint<FReal>(FReal(idxLeafX)*dimLeaf + 2*quarterDimLeaf,
                                                   quarterDimLeaf,
                                                   quarterDimLeaf);
                    /*     particles[1].position = FPoint<FReal>(FReal(idxLeafX)*dimLeaf + 2*quarterDimLeaf,
                               quarterDimLeaf,
                               quarterDimLeaf);*/

                    // Create octree
                    OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, boxCenter);
                    for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
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
                    //Print("Fmm...");
                    KernelClass kernels(DevP,NbLevels,boxWidth, boxCenter);
                    FmmClass algo(&tree,&kernels);
                    algo.execute();

                    // Run direct computation
                    //Print("Direct...");
                    for(int idxTarget = 0 ; idxTarget < nbParticles ; ++idxTarget){
                        for(int idxOther = idxTarget + 1 ; idxOther < nbParticles ; ++idxOther){
                            FP2PR::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
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
                    //Print("Compute Diff...");
                    FMath::FAccurater<FReal> potentialDiff;
                    FMath::FAccurater<FReal> fx, fy, fz;
                    { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                            const FReal*const potentials = leaf->getTargets()->getPotentials();
                            const FReal*const forcesX = leaf->getTargets()->getForcesX();
                            const FReal*const forcesY = leaf->getTargets()->getForcesY();
                            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                            const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                                const FSize indexPartOrig = indexes[idxPart];
                                potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                                fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                                fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                                fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);

                                std::cout <<"indx: " << indexPartOrig << " Potential  Direct " << particles[indexPartOrig].potential<<" FMM "<<potentials[idxPart] << std::endl;
                                std::cout <<"indx: " << indexPartOrig << " Fx Direct " << particles[indexPartOrig].forces[0]<<" FMM "<<forcesX[idxPart] << std::endl;
                                std::cout <<"indx: " << indexPartOrig << " Fy Direct " << particles[indexPartOrig].forces[1]<<" FMM "<<forcesY[idxPart] << std::endl;
                                std::cout <<"indx: " << indexPartOrig << " Fz Direct " << particles[indexPartOrig].forces[2]<<" FMM "<<forcesZ[idxPart] << std::endl;
                                printf("Printf : forces x %e y %e z %e\n", forcesX[idxPart],forcesY[idxPart],forcesZ[idxPart]);//TODO delete
                            }
                        });

                        tree.forEachCell([&](CellClass* cell){

                            std::cout << "cell " << cell->getMortonIndex() << "\n   Multipole: \n";
                            int index_j_k = 0;
                            for (int j = 0 ; j <= DevP ; ++j ){
                                std::cout <<"[" << j << "] ----- level\n";
                                for (int k=0; k<=j ;++k, ++index_j_k){
                                    std::cout << "[" << k << "] ( " << cell->getMultipole()[index_j_k].getReal() << " , " << cell->getMultipole()[index_j_k].getImag() << " )   ";
                                }
                                std::cout << "\n";
                            }
                            std::cout << "\n";
                            std::cout << "   Local:\n";
                            index_j_k = 0;
                            for (int j = 0 ; j <= DevP ; ++j ){
                                std::cout <<"[" << j << "] ----- level \n";
                                for (int k=0; k<=j ;++k, ++index_j_k){
                                    std::cout << "[" << k << "] ( " << cell->getLocal()[index_j_k].getReal() << " , " << cell->getLocal()[index_j_k].getImag() << " )   ";
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
                    //                   const FReal MaximumDiff = FReal(0.0001);
                    //                    if(fx.getL2Norm() > MaximumDiff || fx.getInfNorm() > MaximumDiff){
                    //                        std::cout << "Error in X " << fx.getL2Norm() << " " << fx.getInfNorm() << std::endl;
                    //                    }
                    //                    if(fy.getL2Norm() > MaximumDiff || fy.getInfNorm() > MaximumDiff){
                    //                        std::cout << "Error in Y " << fy.getL2Norm() << " " << fy.getInfNorm() << std::endl;
                    //                    }
                    //                    if(fz.getL2Norm() > MaximumDiff || fz.getInfNorm() > MaximumDiff){
                    //                        std::cout << "Error in Z " << fz.getL2Norm() << " " << fz.getInfNorm() << std::endl;
                    //                    }
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

    /** Classic */
    void TestSpherical(){
        typedef double FReal;

        typedef FSphericalCell<FReal>            CellClass;
        typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

        typedef FSphericalKernel< FReal, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
        std::cout << std::endl << std::endl << std::endl
                  << " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl
                  << " $$$$$$$$$                                                    TestSpherical                                             $$$$$$$$$$$$$$$$"<<std::endl
                  << " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;

        RunTest< FReal,  CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(false);
    }

    /** Rotation */
    void TestRotation(){
        typedef double FReal;

        typedef FSphericalCell<FReal>            CellClass;
        typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

        typedef FSphericalRotationKernel<  FReal, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest< FReal, CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(false);
    }

#ifdef SCALFMM_USE_BLAS
    /** Blas */
    void TestSphericalBlas(){
        typedef double FReal;

        typedef FSphericalCell<FReal>            CellClass;
        typedef FP2PParticleContainerIndexed<FReal>  ContainerClass;

        typedef FSphericalBlasKernel<FReal, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest< FReal,  CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(true);
    }

    /** Block blas */
    void TestSphericalBlockBlas(){
        typedef double FReal;

        typedef FSphericalCell<FReal>            CellClass;
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;

        typedef FSphericalBlockBlasKernel< FReal, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
        //
        std::cout << std::endl << std::endl << std::endl
                  << " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl
                  << " $$$$$$$$$                                                    TestSphericalBlockBlas                                $$$$$$$$$$$$$$$$"<<std::endl
                  << " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
        RunTest<  FReal, CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(true);
    }
#endif

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        //    AddTest(&TestSphericalDirect::TestSpherical,"Test Spherical Kernel");
        //AddTest(&TestSphericalDirect::TestRotation,"Test Rotation Spherical Kernel");
#ifdef SCALFMM_USE_BLAS
        AddTest(&TestSphericalDirect::TestSphericalBlas,"Test Spherical Blas Kernel");
        //        AddTest(&TestSphericalDirect::TestSphericalBlockBlas,"Test Spherical Block Blas Kernel");
#endif
    }
};


// You must do this
TestClass(TestSphericalDirect)



