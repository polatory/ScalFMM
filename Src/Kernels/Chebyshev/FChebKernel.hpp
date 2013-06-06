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
#ifndef FCHEBKERNEL_HPP
#define FCHEBKERNEL_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FSmartPointer.hpp"

#include "./FAbstractChebKernel.hpp"
#include "./FChebM2LHandler.hpp"

class FTreeCoordinate;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Chebyshev interpolation based FMM operators. It
 * implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER>
class FChebKernel
    : public FAbstractChebKernel< CellClass, ContainerClass, MatrixKernelClass, ORDER>
{
	// private types
	typedef FChebM2LHandler<ORDER,MatrixKernelClass> M2LHandlerClass;

	// using from 
    typedef FAbstractChebKernel< CellClass, ContainerClass, MatrixKernelClass, ORDER>
	AbstractBaseClass;

	/// Needed for M2L operator
	FSmartPointer<  M2LHandlerClass,FSmartPointerMemory> M2LHandler;

public:
	/**
	 * The constructor initializes all constant attributes and it reads the
	 * precomputed and compressed M2L operators from a binary file (an
	 * runtime_error is thrown if the required file is not valid).
	 */
	FChebKernel(const int inTreeHeight,
							const FPoint& inBoxCenter,
							const FReal inBoxWidth,
							const FReal Epsilon)
        : FAbstractChebKernel< CellClass, ContainerClass, MatrixKernelClass, ORDER>(inTreeHeight,
																																															inBoxCenter,
																																															inBoxWidth),
			M2LHandler(new M2LHandlerClass(Epsilon))
	{
		// read precomputed compressed m2l operators from binary file
		M2LHandler->ReadFromBinaryFileAndSet();
		//M2LHandler->ComputeAndCompressAndSet();
	}


	void P2M(CellClass* const LeafCell,
					 const ContainerClass* const SourceParticles)
	{
		// 1) apply Sy
        const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));
		AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter,
																							AbstractBaseClass::BoxWidthLeaf,
																							LeafCell->getMultipole(),
																							SourceParticles);
		// 2) apply B
		M2LHandler->applyB(LeafCell->getMultipole(),
											 LeafCell->getMultipole() + AbstractBaseClass::nnodes);
	}


	void M2M(CellClass* const FRestrict ParentCell,
					 const CellClass*const FRestrict *const FRestrict ChildCells,
					 const int TreeLevel)
	{
		// 1) apply Sy
		FBlas::scal(AbstractBaseClass::nnodes*2, FReal(0.), ParentCell->getMultipole());
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				AbstractBaseClass::Interpolator->applyM2M(ChildIndex,
																									ChildCells[ChildIndex]->getMultipole(),
																									ParentCell->getMultipole());
		// 2) apply B
		M2LHandler->applyB(ParentCell->getMultipole(),
											 ParentCell->getMultipole() + AbstractBaseClass::nnodes);
	}


//	void M2L(CellClass* const FRestrict TargetCell,
//					 const CellClass* SourceCells[343],
//					 const int NumSourceCells,
//					 const int TreeLevel) const
//	{
//		const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
//		const FTreeCoordinate& cx = TargetCell->getCoordinate();
//		for (int idx=0; idx<NumSourceCells; ++idx) {
//			const FTreeCoordinate& cy = SourceCells[idx]->getCoordinate();
//			const int transfer[3] = {cy.getX()-cx.getX(),
//															 cy.getY()-cx.getY(),
//															 cy.getZ()-cx.getZ()};
//			M2LHandler->applyC(transfer, CellWidth,
//												SourceCells[idx]->getMultipole() + AbstractBaseClass::nnodes,
//												TargetCell->getLocal() + AbstractBaseClass::nnodes);
//		}
//	}

	void M2L(CellClass* const FRestrict TargetCell,
					 const CellClass* SourceCells[343],
					 const int NumSourceCells,
					 const int TreeLevel)
	{
		FReal *const CompressedLocalExpansion = TargetCell->getLocal() + AbstractBaseClass::nnodes;
		const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(FMath::pow(2, TreeLevel)));
		for (int idx=0; idx<343; ++idx)
			if (SourceCells[idx])
				M2LHandler->applyC(idx, CellWidth,
													 SourceCells[idx]->getMultipole() + AbstractBaseClass::nnodes,
													 CompressedLocalExpansion);
	}

//	void M2L(CellClass* const FRestrict TargetCell,
//					 const CellClass* SourceCells[343],
//					 const int NumSourceCells,
//					 const int TreeLevel) const
//	{
//		const unsigned int rank = M2LHandler.getRank();
//		FBlas::scal(343*rank, FReal(0.), MultipoleExpansion);
//		const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
//		for (int idx=0; idx<343; ++idx)
//			if (SourceCells[idx])
//				FBlas::copy(rank, const_cast<FReal *const>(SourceCells[idx]->getMultipole())+AbstractBaseClass::nnodes,
//										MultipoleExpansion+idx*rank);
//		
//		M2LHandler->applyC(CellWidth, MultipoleExpansion, TargetCell->getLocal() + AbstractBaseClass::nnodes);
//	}


	void L2L(const CellClass* const FRestrict ParentCell,
					 CellClass* FRestrict *const FRestrict ChildCells,
					 const int TreeLevel)
	{
		// 1) apply U
		M2LHandler->applyU(ParentCell->getLocal() + AbstractBaseClass::nnodes,
											 const_cast<CellClass*>(ParentCell)->getLocal());
		// 2) apply Sx
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				AbstractBaseClass::Interpolator->applyL2L(ChildIndex,
																									ParentCell->getLocal(),
																									ChildCells[ChildIndex]->getLocal());
	}

	void L2P(const CellClass* const LeafCell,
					 ContainerClass* const TargetParticles)
	{
		// 1) apply U
		M2LHandler->applyU(LeafCell->getLocal() + AbstractBaseClass::nnodes,
											 const_cast<CellClass*>(LeafCell)->getLocal());
		
        const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));

		//// 2.a) apply Sx
		//AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter,
		//																					AbstractBaseClass::BoxWidthLeaf,
		//																					LeafCell->getLocal(),
		//																					TargetParticles);
		//// 2.b) apply Px (grad Sx)
		//AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter,
		//																									AbstractBaseClass::BoxWidthLeaf,
		//																									LeafCell->getLocal(),
		//																									TargetParticles);

		// 2.c) apply Sx and Px (grad Sx)
		AbstractBaseClass::Interpolator->applyL2PTotal(LeafCellCenter,
																									 AbstractBaseClass::BoxWidthLeaf,
																									 LeafCell->getLocal(),
																									 TargetParticles);
	}

};


#endif //FCHEBKERNELS_HPP

// [--END--]
