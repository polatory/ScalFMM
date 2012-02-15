#ifndef FCHEBKERNELS_HPP
#define FCHEBKERNELS_HPP
// [--License--]

#include "../Utils/FGlobal.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FSmartPointer.hpp"

#include "./FChebInterpolator.hpp"
#include "./FChebM2LHandler.hpp"

class FTreeCoordinate;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebKernels
 * @brief
 * Please read the license
 *
 * This kernels implement the Chebyshev interpolation based FMM operators. It
 * implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam ParticleClass Type of particle
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template<class ParticleClass,
				 class CellClass,
				 class ContainerClass,
				 class MatrixKernelClass,
				 int ORDER>
class FChebKernels
{
  enum {nnodes = TensorTraits<ORDER>::nnodes};
	typedef FChebInterpolator<ORDER> InterpolatorClass;
	typedef FChebM2LHandler<ORDER,MatrixKernelClass> M2LHandlerClass;

	/// Needed for P2M, M2M, L2L and L2P operators
	FSmartPointer<InterpolatorClass,FSmartPointerMemory> Interpolator;
	/// Needed for M2L operator
	FSmartPointer<  M2LHandlerClass,FSmartPointerMemory> M2LHandler;
	/// Needed for P2P operator
	FSmartPointer<MatrixKernelClass,FSmartPointerMemory> MatrixKernel;
	/// Height of the entire oct-tree
	const unsigned int TreeHeight;
	/// Corner of oct-tree box
	const F3DPosition BoxCorner;
	/// Width of oct-tree box   
	const FReal BoxWidth;    
	/// Width of a leaf cell box 
	const FReal BoxWidthLeaf;
	/// Prescribed accuracy of the compressed M2L operators
	const FReal Epsilon;        

	/**
	 * Compute center of leaf cell from its tree coordinate.
	 * @param[in] Coordinate tree coordinate
	 * @return center of leaf cell
	 */
	const F3DPosition getLeafCellCenter(const FTreeCoordinate& Coordinate) const
	{
		return F3DPosition(BoxCorner.getX() + (FReal(Coordinate.getX()) + FReal(.5)) * BoxWidthLeaf,
											 BoxCorner.getY() + (FReal(Coordinate.getY()) + FReal(.5)) * BoxWidthLeaf,
											 BoxCorner.getZ() + (FReal(Coordinate.getZ()) + FReal(.5)) * BoxWidthLeaf);
	}

public:
	/**
	 * The constructor initializes all constant attributes and it reads the
	 * precomputed and compressed M2L operators from a binary file (an
	 * runtime_error is thrown if the required file is not valid).
	 */
	FChebKernels(const int inTreeHeight,
							 const F3DPosition& inBoxCenter,
							 const FReal inBoxWidth,
							 const FReal inEpsilon)
		: Interpolator(new InterpolatorClass()),
			M2LHandler(new M2LHandlerClass(inEpsilon)),
			MatrixKernel(new MatrixKernelClass()),
			TreeHeight(inTreeHeight),
			BoxCorner(inBoxCenter - inBoxWidth / FReal(2.)),
			BoxWidth(inBoxWidth),
			BoxWidthLeaf(BoxWidth / FReal(FMath::pow(2, inTreeHeight - 1))),
			Epsilon(inEpsilon)
	{
		// read precomputed compressed m2l operators from binary file
		M2LHandler->ReadFromBinaryFileAndSet();
	}

//	/** Default destructor */
//	~FChebKernels()	{}


	void P2M(CellClass* const LeafCell,
					 const ContainerClass* const SourceParticles) const
	{
		// 1) apply Sy
		const F3DPosition LeafCellCenter(getLeafCellCenter(LeafCell->getCoordinate()));
		Interpolator->applyP2M(LeafCellCenter,
													 BoxWidthLeaf,
													 LeafCell->getMultipole(),
													 SourceParticles);
		// 2) apply B
		M2LHandler->applyB(LeafCell->getMultipole(),
											 LeafCell->getMultipole() + nnodes);
	}


	void M2M(CellClass* const FRestrict ParentCell,
					 const CellClass*const FRestrict *const FRestrict ChildCells,
					 const int TreeLevel) const
	{
		// 1) apply Sy
		FBlas::scal(nnodes*2, FReal(0.), ParentCell->getMultipole());
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				Interpolator->applyM2M(ChildIndex,
															 ChildCells[ChildIndex]->getMultipole(),
															 ParentCell->getMultipole());
		// 2) apply B
		M2LHandler->applyB(ParentCell->getMultipole(),
											 ParentCell->getMultipole() + nnodes);
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
//												SourceCells[idx]->getMultipole() + nnodes,
//												TargetCell->getLocal() + nnodes);
//		}
//	}

	void M2L(CellClass* const FRestrict TargetCell,
					 const CellClass* SourceCells[343],
					 const int NumSourceCells,
					 const int TreeLevel) const
	{
		FReal *const CompressedLocalExpansion = TargetCell->getLocal() + nnodes;
		const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
		for (int idx=0; idx<343; ++idx)
			if (SourceCells[idx])
				M2LHandler->applyC(idx, CellWidth,
													 SourceCells[idx]->getMultipole() + nnodes,
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
//				FBlas::copy(rank, const_cast<FReal *const>(SourceCells[idx]->getMultipole())+nnodes,
//										MultipoleExpansion+idx*rank);
//		
//		M2LHandler->applyC(CellWidth, MultipoleExpansion, TargetCell->getLocal() + nnodes);
//	}


	void L2L(const CellClass* const FRestrict ParentCell,
					 CellClass* FRestrict *const FRestrict ChildCells,
					 const int TreeLevel) const
	{
		// 1) apply U
		M2LHandler->applyU(ParentCell->getLocal() + nnodes,
											const_cast<CellClass *const>(ParentCell)->getLocal());
		// 2) apply Sx
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				Interpolator->applyL2L(ChildIndex,
															 ParentCell->getLocal(),
															 ChildCells[ChildIndex]->getLocal());
	}

	void L2P(const CellClass* const LeafCell,
					 ContainerClass* const TargetParticles) const
	{
		// 1) apply U
		M2LHandler->applyU(LeafCell->getLocal() + nnodes,
											 const_cast<CellClass *const>(LeafCell)->getLocal());

		// 2) apply Sx
		const F3DPosition LeafCellCenter(getLeafCellCenter(LeafCell->getCoordinate()));
		Interpolator->applyL2P(LeafCellCenter,
													 BoxWidthLeaf,
													 LeafCell->getLocal(),
													 TargetParticles);
	}


	void P2P(const FTreeCoordinate& LeafCellCoordinate,
					 ContainerClass* const FRestrict TargetParticles,
					 const ContainerClass* const FRestrict SourceParticles,
					 const ContainerClass* const NeighborSourceParticles[27],
					 const int ) const
	{
		// loop: target particles
		typename ContainerClass::BasicIterator iTargets(*TargetParticles);
		while (iTargets.hasNotFinished()) {
			ParticleClass& Target = iTargets.data();

			{ // loop: source particles (target leaf cell == source leaf cell)
				typename ContainerClass::ConstBasicIterator iSources(*SourceParticles);
				while (iSources.hasNotFinished()) {
					const ParticleClass& Source = iSources.data();
					// only if target and source are not identical
					if (&Target != &Source)
						Target.incPotential(MatrixKernel->evaluate(Target.getPosition(), Source.getPosition())
																* Source.getPhysicalValue());
					// progress sources
					iSources.gotoNext();
				}
			}
			
			{ // loop: source particles (target leaf cell != source leaf cell)
				for (unsigned int idx=0; idx<27; ++idx) {
					if (NeighborSourceParticles[idx]) {
						typename ContainerClass::ConstBasicIterator	iSources(*NeighborSourceParticles[idx]);
						while (iSources.hasNotFinished()) {
							const ParticleClass& Source = iSources.data();
							// target and source cannot be identical
							Target.incPotential(MatrixKernel->evaluate(Target.getPosition(), Source.getPosition())
																	* Source.getPhysicalValue());
							// progress sources
							iSources.gotoNext();
						}
					}
				}
			}
			
			// progress targets
			iTargets.gotoNext();
		}
	}

};


#endif //FCHEBKERNELS_HPP

// [--END--]
