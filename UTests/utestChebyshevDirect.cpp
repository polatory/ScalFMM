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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"

#include "../Src/Kernels/Chebyshev/FChebParticle.hpp"
#include "../Src/Kernels/Chebyshev/FChebLeaf.hpp"
#include "../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../Src/Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

/*
  In this test we compare the spherical fmm results and the direct results.
*/

/** We need to know the position of the particle in the array */
class IndexedParticle : public FChebParticle {
	int index;
public:
	IndexedParticle(): index(-1){}

	int getIndex() const{
		return index;
	}
	void setIndex( const int inIndex ){
		index = inIndex;
	}
};

/** the test class
 *
 */
class TestChebyshevDirect : public FUTester<TestChebyshevDirect> {
	
	///////////////////////////////////////////////////////////
	// The tests!
	///////////////////////////////////////////////////////////
	
	template <class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
						class LeafClass, class OctreeClass, class FmmClass>
	void RunTest(const FReal epsilon)	{
		// Warning in make test the exec dir it Build/UTests
		// Load particles
        const char* const filename = (sizeof(FReal) == sizeof(float))?
                                        "../../Data/utestDirect.bin.fma.single":
                                        "../../Data/utestDirect.bin.fma.double";
        FFmaBinLoader<ParticleClass> loader(filename);
		if(!loader.isOpen()){
			Print("Cannot open particles file.");
			uassert(false);
			return;
		}
		Print("Number of particles:");
		Print(loader.getNumberOfParticles());

		const int NbLevels      = 4;
		const int SizeSubLevels = 2;
		//const FReal epsilon = FReal(1e-5);

		// Create octree
		OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
		ParticleClass* const particles = new ParticleClass[loader.getNumberOfParticles()];
		for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
			loader.fillParticle(particles[idxPart]);
			particles[idxPart].setIndex( idxPart );
			tree.insert(particles[idxPart]);
		}


		// Run FMM
		Print("Fmm...");
		KernelClass kernels(NbLevels, loader.getCenterOfBox(), loader.getBoxWidth(), epsilon);
		FmmClass algo(&tree,&kernels);
		algo.execute();

		// Run direct computation
		const MatrixKernelClass MatrixKernel;
		Print("Direct...");
		for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
			for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
				//kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
				const FReal wt = particles[idxTarget].getPhysicalValue();
				const FReal ws = particles[idxOther ].getPhysicalValue();
				const FReal one_over_r = MatrixKernel.evaluate(particles[idxTarget].getPosition(),
																											 particles[idxOther].getPosition());
				// potential
				particles[idxTarget].incPotential(one_over_r * ws);
				particles[idxOther ].incPotential(one_over_r * wt);
				// force
				FPoint force(particles[idxOther].getPosition() - particles[idxTarget].getPosition());
				force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
				particles[idxTarget].incForces(  force.getX(),  force.getY(),  force.getZ());
				particles[idxOther ].incForces( -force.getX(), -force.getY(), -force.getZ());
			}
		}

		// Compare
		Print("Compute Diff...");
		FMath::FAccurater potentialDiff;
		FMath::FAccurater fx, fy, fz;
		{ // Check that each particle has been summed with all other
			typename OctreeClass::Iterator octreeIterator(&tree);
			octreeIterator.gotoBottomLeft();

			do{
				typename ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

				while( leafIter.hasNotFinished() ){
					const ParticleClass& other = particles[leafIter.data().getIndex()];

					potentialDiff.add(other.getPotential(),leafIter.data().getPotential());

					fx.add(other.getForces().getX(),leafIter.data().getForces().getX());

					fy.add(other.getForces().getY(),leafIter.data().getForces().getY());

					fz.add(other.getForces().getZ(),leafIter.data().getForces().getZ());

					leafIter.gotoNext();
				}
			} while(octreeIterator.moveRight());
		}

		delete[] particles;

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
        const FReal MaximumDiffPotential = FReal(9e-5);
        const FReal MaximumDiffForces = FReal(9e-3);

		uassert(potentialDiff.getL2Norm() < MaximumDiffPotential);
		uassert(potentialDiff.getInfNorm() < MaximumDiffPotential);
		uassert(fx.getL2Norm()  < MaximumDiffForces);
		uassert(fx.getInfNorm() < MaximumDiffForces);
		uassert(fy.getL2Norm()  < MaximumDiffForces);
		uassert(fy.getInfNorm() < MaximumDiffForces);
		uassert(fz.getL2Norm()  < MaximumDiffForces);
		uassert(fz.getInfNorm() < MaximumDiffForces);
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
		const unsigned int ORDER = 5;
		const FReal epsilon = FReal(1e-5);
		typedef IndexedParticle ParticleClass;
		typedef FVector<ParticleClass> ContainerClass;
		typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
		typedef FChebMatrixKernelR MatrixKernelClass;
		typedef FChebCell<ORDER> CellClass;
		typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
		typedef FChebKernel<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithm<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		// run test
		RunTest<ParticleClass,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(epsilon);
	}

	/** TestChebSymKernel */
	void TestChebSymKernel(){
		const unsigned int ORDER = 5;
		const FReal epsilon = FReal(1e-5);
		typedef IndexedParticle ParticleClass;
		typedef FVector<ParticleClass> ContainerClass;
		typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
		typedef FChebMatrixKernelR MatrixKernelClass;
		typedef FChebCell<ORDER> CellClass;
		typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
		typedef FChebSymKernel<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithm<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		// run test
		RunTest<ParticleClass,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(epsilon);
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




