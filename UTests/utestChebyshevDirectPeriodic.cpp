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
// ================

#include "Utils/FGlobal.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FTreeIO.hpp"

#include "Core/FFmmAlgorithmPeriodic.hpp"

#include "FUTester.hpp"


#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/*
  In this test we compare the Chebyshev fmm results and the direct results.
 */


/** the test class
 *
 */
class TestChebyshevDirect : public FUTester<TestChebyshevDirect> {

	///////////////////////////////////////////////////////////
	// The tests!
	///////////////////////////////////////////////////////////

	template <class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
	class LeafClass, class OctreeClass, class FmmClass>
	void RunTest()	{
		// Warning in make test the exec dir it Build/UTests
		// Load particles

		const int NbLevels        = 4;
		const int SizeSubLevels = 2;
		const int PeriodicDeep  = 2;
		const int NbParticles     = 250;

		FRandomLoader<FReal> loader(NbParticles);

		Print("Number of particles:");
		Print(loader.getNumberOfParticles());

		// Create octree
		OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // interaction kernel evaluator
        typedef double FReal;
    const MatrixKernelClass MatrixKernel;

		struct TestParticle{
            FPoint<FReal> position;
			FReal forces[3];
			FReal physicalValue;
			FReal potential;
		};
		FReal coeff = -1.0, value = 0.10, sum = 0.0;
		TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
		for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> position;
			loader.fillParticle(&position);
			// put in tree
			value *= coeff ;
			sum += value ;
			// put in tree
			tree.insert(position, idxPart, value);
			// get copy
			particles[idxPart].position         = position;
			particles[idxPart].physicalValue = value;
			particles[idxPart].potential        = 0.0;
			particles[idxPart].forces[0]        = 0.0;
			particles[idxPart].forces[1]        = 0.0;
			particles[idxPart].forces[2]        = 0.0;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Run FMM computation
		/////////////////////////////////////////////////////////////////////////////////////////////////
		Print("Fmm...");
		FmmClass algo(&tree,PeriodicDeep );
		KernelClass kernels(algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter(),&MatrixKernel);
		algo.setKernel(&kernels);
		algo.execute();
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Run direct computation
		/////////////////////////////////////////////////////////////////////////////////////////////////
		Print("Direct...");

		FTreeCoordinate min, max;
		algo.repetitionsIntervals(&min, &max);
		FReal energy= 0.0 , energyD = 0.0 ;

		for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
			for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
				FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
						particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
						&particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
						&particles[idxTarget].forces[2],&particles[idxTarget].potential,
						particles[idxOther].position.getX(), particles[idxOther].position.getY(),
						particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
						&particles[idxOther].forces[0],&particles[idxOther].forces[1],
                              &particles[idxOther].forces[2],&particles[idxOther].potential,&MatrixKernel);

			}
			for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
				for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
					for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
						if(idxX ==0 && idxY == 0 && idxZ == 0) continue;
						// next lines for test

                        const FPoint<FReal> offset(loader.getBoxWidth() * FReal(idxX),
								loader.getBoxWidth() * FReal(idxY),
								loader.getBoxWidth() * FReal(idxZ));

						for(int idxSource = 0 ; idxSource < NbParticles ; ++idxSource){
							TestParticle source = particles[idxSource];
							source.position += offset;

							FP2P::NonMutualParticles(
									source.position.getX(), source.position.getY(),
									source.position.getZ(),source.physicalValue,
									particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
									particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
									&particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
									&particles[idxTarget].forces[2],&particles[idxTarget].potential,&MatrixKernel);
						}
					}
				}
			}
		}
		for(int idx = 0 ; idx <  loader.getNumberOfParticles()  ; ++idx){
			energyD +=  particles[idx].potential*particles[idx].physicalValue ;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Compare
		/////////////////////////////////////////////////////////////////////////////////////////////////

		Print("Compute Diff...");
		FMath::FAccurater<FReal> potentialDiff;
		FMath::FAccurater<FReal> fx, fy, fz;

		{ // Check that each particle has been summed with all other

			tree.forEachLeaf([&](LeafClass* leaf){
				const FReal*const potentials        = leaf->getTargets()->getPotentials();
				const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
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
					energy+=potentials[idxPart]*physicalValues[idxPart];
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


        uassert(potentialDiff.getL2Norm() < MaximumDiffPotential);    //1
        uassert(potentialDiff.getRMSError() < MaximumDiffPotential);  //2
        uassert(fx.getL2Norm()  < MaximumDiffForces);                       //3
        uassert(fx.getRMSError() < MaximumDiffForces);                      //4
        uassert(fy.getL2Norm()  < MaximumDiffForces);                       //5
        uassert(fy.getRMSError() < MaximumDiffForces);                      //6
        uassert(fz.getL2Norm()  < MaximumDiffForces);                      //8
        uassert(fz.getRMSError() < MaximumDiffForces);                                           //8
        uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
        uassert(FMath::Abs(energy-energyD) < 10*MaximumDiffPotential);                     //10  Total Energy

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
		const unsigned int ORDER = 6;
		typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
		typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
		typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
		typedef FChebCell<FReal,ORDER> CellClass;
		typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
		typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithmPeriodic<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		// run test
		std::cout <<" TEST 1  "<<std::endl;
		RunTest<CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>();
	}

	/** TestChebSymKernel */
	void TestChebSymKernel(){
		const unsigned int ORDER = 7;
		typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
		typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
		typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
		typedef FChebCell<FReal,ORDER> CellClass;
		typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
		typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithmPeriodic<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		// run test
		std::cout <<std::endl<<" TEST 2 "<<std::endl;
		RunTest<CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>();
	}



	///////////////////////////////////////////////////////////
	// Set the tests!
	///////////////////////////////////////////////////////////

	/** set test */
	void SetTests(){
		AddTest(&TestChebyshevDirect::TestChebKernel,"Test Chebyshev Kernel with one big SVD") ;
		AddTest(&TestChebyshevDirect::TestChebSymKernel,"Test Chebyshev Kernel with 16 small SVDs and symmetries");
	}
};


// You must do this
TestClass(TestChebyshevDirect)




