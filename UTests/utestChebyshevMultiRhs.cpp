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

// ==== CMAKE =====
// @FUSE_BLAS
// ==============

#include "ScalFmmConfig.h"
#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Containers/FOctree.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"


#include "../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"
/*
  In this test we compare the spherical FMM results and the direct results.
*/


/** the test class
 *
 */
class TestChebyshevDirect : public FUTester<TestChebyshevDirect> {

    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    template <class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
                        class LeafClass, class OctreeClass, class FmmClass, const int NVals>
    void RunTest()	{
        // Warning in make test the exec dir it Build/UTests
        // Load particles
		//
		// Load particles
		//
		if(sizeof(FReal) == sizeof(float) ) {
			std::cerr << "No input data available for Float "<< std::endl;
			exit(EXIT_FAILURE);
		}
		const std::string parFile( (sizeof(FReal) == sizeof(float))?
				"Test/DirectFloat.bfma":
				"UTest/DirectDouble.bfma");
		//
		std::string filename(SCALFMMDataPath+parFile);
		//
		FFmaGenericLoader loader(filename);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            uassert(false);
            return;
        }
        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels        = 4;
        const int SizeSubLevels = 2;

        //
		FSize nbParticles = loader.getNumberOfParticles() ;
		FmaRWParticle<8,8>* const particles = new FmaRWParticle<8,8>[nbParticles];

		loader.fillParticle(particles,nbParticles);
         //
		// Create octree
		OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
		//   Insert particle in the tree
		//
		for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		    tree.insert(particles[idxPart].getPosition() , idxPart, particles[idxPart].getPhysicalValue() ,0.0,0.0,0.0);
		}
//		//
//        // Create octree
//        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
//            FPoint position;
//            FReal physicalValue;
//            loader.fillParticle(&position,&physicalValue);
//            // put in tree
//            tree.insert(position, idxPart, physicalValue);
//            // get copy
//            particles[idxPart].position = position;
//            particles[idxPart].physicalValue = physicalValue;
//            particles[idxPart].potential = 0.0;
//            particles[idxPart].forces[0] = 0.0;
//            particles[idxPart].forces[1] = 0.0;
//            particles[idxPart].forces[2] = 0.0;
//        }


        // Run FMM
        Print("Fmm...");
        KernelClass kernels(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&tree,&kernels);
        algo.execute();
//
		FReal energy= 0.0 , energyD = 0.0 ;

        // Run direct computation
        Print("Direct...");
	for(int idx = 0 ; idx <  loader.getNumberOfParticles()  ; ++idx){
	    *(particles[idx].setPotential()) *= NVals ;
	    particles[idx].getForces()[0] *= NVals ;
	    particles[idx].getForces()[1] *= NVals ;
	    particles[idx].getForces()[2] *= NVals ;
	    energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
	}
	//
        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
		{ // Check that each particle has been summed with all other

			tree.forEachLeaf([&](LeafClass* leaf){
				const FReal*const potentials        = leaf->getTargets()->getPotentials();
				const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
				const FReal*const forcesX            = leaf->getTargets()->getForcesX();
				const FReal*const forcesY            = leaf->getTargets()->getForcesY();
				const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
				const int nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
				const FVector<int>& indexes      = leaf->getTargets()->getIndexes();

				for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
					const int indexPartOrig = indexes[idxPart];
					potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
					fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
					fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
					fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
					energy   += potentials[idxPart]*physicalValues[idxPart];
				}
			});
		}

		delete[] particles;

        // Print for information
		Print("Potential diff is = ");
		printf("         Pot L2Norm     %e\n",potentialDiff.getL2Norm());
		printf("         Pot RL2Norm   %e\n",potentialDiff.getRelativeL2Norm());
		printf("         Pot RMSError   %e\n",potentialDiff.getRMSError());
		Print("Fx diff is = ");
		printf("         Fx L2Norm     %e\n",fx.getL2Norm());
		printf("         Fx RL2Norm   %e\n",fx.getRelativeL2Norm());
		printf("         Fx RMSError   %e\n",fx.getRMSError());
		Print("Fy diff is = ");
		printf("        Fy L2Norm     %e\n",fy.getL2Norm());
		printf("        Fy RL2Norm   %e\n",fy.getRelativeL2Norm());
		printf("        Fy RMSError   %e\n",fy.getRMSError());
		Print("Fz diff is = ");
		printf("        Fz L2Norm     %e\n",fz.getL2Norm());
		printf("        Fz RL2Norm   %e\n",fz.getRelativeL2Norm());
		printf("        Fz RMSError   %e\n",fz.getRMSError());
		FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
		printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
		printf("  Energy Error  =   %.12e\n",FMath::Abs(energy-energyD));
		printf("  Energy FMM    =   %.12e\n",FMath::Abs(energy));
		printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyD));

        // Assert
        const FReal MaximumDiffPotential = FReal(9e-3);
        const FReal MaximumDiffForces     = FReal(9e-2);

		Print("Test1 - Error Relative L2 norm Potential ");
		uassert(potentialDiff.getRelativeL2Norm() < MaximumDiffPotential);    //1
		Print("Test2 - Error RMS L2 norm Potential ");
		uassert(potentialDiff.getRMSError() < MaximumDiffPotential);  //2
		Print("Test3 - Error Relative L2 norm FX ");
		uassert(fx.getRelativeL2Norm()  < MaximumDiffForces);                       //3
		Print("Test4 - Error RMS L2 norm FX ");
		uassert(fx.getRMSError() < MaximumDiffForces);                      //4
		Print("Test5 - Error Relative L2 norm FY ");
		uassert(fy.getRelativeL2Norm()  < MaximumDiffForces);                       //5
		Print("Test6 - Error RMS L2 norm FY ");
		uassert(fy.getRMSError() < MaximumDiffForces);                      //6
		Print("Test7 - Error Relative L2 norm FZ ");
		uassert(fz.getRelativeL2Norm()  < MaximumDiffForces);                      //8
		Print("Test8 - Error RMS L2 norm FZ ");
		uassert(fz.getRMSError() < MaximumDiffForces);                                           //8
		Print("Test9 - Error Relative L2 norm F ");
		uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
		Print("Test10 - Relative error Energy ");
		uassert(FMath::Abs(energy-energyD) /energyD< MaximumDiffPotential);                     //10  Total Energy


        // Compute multipole local rhs diff
        FMath::FAccurater localDiff;
        FMath::FAccurater multiPoleDiff;
        tree.forEachCell([&](CellClass* cell){
            for( int idxRhs = 1 ; idxRhs < NVals ; ++idxRhs){
                localDiff.add(cell->getLocal(0), cell->getLocal(idxRhs), cell->getVectorSize());
                multiPoleDiff.add(cell->getMultipole(0), cell->getMultipole(idxRhs), cell->getVectorSize());
            }
        });
        Print("Local diff is = ");
        Print(localDiff.getL2Norm());
        Print(localDiff.getInfNorm());
        Print("Multipole diff is = ");
        Print(multiPoleDiff.getL2Norm());
        Print(multiPoleDiff.getInfNorm());

        uassert(localDiff.getL2Norm()  < 1e-10);
        uassert(localDiff.getInfNorm() < 1e-10);
        uassert(multiPoleDiff.getL2Norm()  < 1e-10);
        uassert(multiPoleDiff.getInfNorm() < 1e-10);
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                                << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated()
                                << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated()
                                << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }


    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////


    /** TestChebKernel */
    void TestChebKernel(){
        const int NVals = 4;
        const unsigned int ORDER = 6 ;
        typedef FP2PParticleContainerIndexed<> ContainerClass;
        typedef FSimpleLeaf<ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR MatrixKernelClass;
        typedef FChebCell<ORDER, 1, 1, NVals> CellClass;
        typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER, NVals> KernelClass;
        typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        // run test
        RunTest<CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass, NVals>();
    }

    /** TestChebSymKernel */
    void TestChebSymKernel(){
        const int NVals = 4;
        const unsigned int ORDER = 6;
        typedef FP2PParticleContainerIndexed<> ContainerClass;
        typedef FSimpleLeaf<ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR MatrixKernelClass;
        typedef FChebCell<ORDER, 1, 1, NVals> CellClass;
        typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER, NVals> KernelClass;
        typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        // run test
        RunTest<CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass, NVals>();
    }



    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestChebyshevDirect::TestChebKernel,"Test Chebyshev Kernel with one big SVD");
        AddTest(&TestChebyshevDirect::TestChebSymKernel,"Test Chebyshev Kernel with 16 small SVDs and symmetries");
    }
};


// You must do this
TestClass(TestChebyshevDirect)





