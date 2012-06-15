#ifndef FCHEBSYMKERNEL_HPP
#define FCHEBSYMKERNEL_HPP
// [--License--]

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FSmartPointer.hpp"

#include "./FAbstractChebKernel.hpp"
#include "./FChebInterpolator.hpp"
#include "./FChebSymM2LHandler.hpp"

class FTreeCoordinate;



/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebSymKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Chebyshev interpolation based FMM operators
 * exploiting the symmetries in the far-field. It implements all interfaces
 * (P2P, P2M, M2M, M2L, L2L, L2P) which are required by the FFmmAlgorithm and
 * FFmmAlgorithmThread.
 *
 * @tparam ParticleClass Type of particle
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template <class ParticleClass, class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER>
class FChebSymKernel
	: public FAbstractChebKernel<ParticleClass, CellClass, ContainerClass, MatrixKernelClass, ORDER>
{
	typedef FAbstractChebKernel<ParticleClass, CellClass, ContainerClass, MatrixKernelClass, ORDER>	AbstractBaseClass;
	typedef SymmetryHandler<ORDER, MatrixKernelClass::Type> SymmetryHandlerClass;
  enum {nnodes = AbstractBaseClass::nnodes};

	/// Needed for handling all symmetries
	const FSmartPointer<SymmetryHandlerClass,FSmartPointerMemory> SymHandler;

	// permuted local and multipole expansions
	FReal** Loc;
	FReal** Mul;
	unsigned int* countExp;


	
	/**
	 * Allocate memory for storing locally permuted mulipole and local expansions
	 */
	void allocateMemoryForPermutedExpansions()
	{
		assert(Loc==NULL && Mul==NULL && countExp==NULL);
		Loc = new FReal* [343];
		Mul = new FReal* [343];
		countExp = new unsigned int [343];

		// set all 343 to NULL
		for (unsigned int idx=0; idx<343; ++idx) {
			Mul[idx] = Loc[idx] = NULL;
		}

		// init only 16 of 343 possible translations due to symmetries
		for (int i=2; i<=3; ++i)
			for (int j=0; j<=i; ++j)
				for (int k=0; k<=j; ++k) {
					const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
					assert(Mul[idx]==NULL || Loc[idx]==NULL);
					Mul[idx] = new FReal [24 * nnodes];
					Loc[idx] = new FReal [24 * nnodes];
				}
	}
	


public:
	/**
	 * The constructor initializes all constant attributes and it reads the
	 * precomputed and compressed M2L operators from a binary file (an
	 * runtime_error is thrown if the required file is not valid).
	 */
	FChebSymKernel(const int inTreeHeight,
								 const FPoint& inBoxCenter,
								 const FReal inBoxWidth,
								 const FReal Epsilon)
		: AbstractBaseClass(inTreeHeight, inBoxCenter, inBoxWidth),
			SymHandler(new SymmetryHandlerClass(AbstractBaseClass::MatrixKernel.getPtr(), Epsilon, inBoxWidth, inTreeHeight)),
			Loc(NULL), Mul(NULL), countExp(NULL)
	{
		this->allocateMemoryForPermutedExpansions();
	}
	
	

	/** Copy constructor */
	FChebSymKernel(const FChebSymKernel& other)
		: AbstractBaseClass(other),
			SymHandler(other.SymHandler),
			Loc(NULL), Mul(NULL), countExp(NULL)
	{
		this->allocateMemoryForPermutedExpansions();
	}



	/** Destructor */
	~FChebSymKernel()
	{
		for (unsigned int t=0; t<343; ++t) {
			if (Loc[t]!=NULL) delete [] Loc[t];
			if (Mul[t]!=NULL) delete [] Mul[t];
		}
		if (Loc!=NULL)      delete [] Loc;
		if (Mul!=NULL)      delete [] Mul;
		if (countExp!=NULL) delete [] countExp;
	}


	const SymmetryHandlerClass *const getPtrToSymHandler() const
	{	return SymHandler.getPtr();	}
	


	void P2M(CellClass* const LeafCell,
					 const ContainerClass* const SourceParticles)
	{
		// apply Sy
		const FPoint LeafCellCenter(getLeafCellCenter(LeafCell->getCoordinate()));
		AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter,
																							AbstractBaseClass::BoxWidthLeaf,
																							LeafCell->getMultipole(),
																							SourceParticles);
	}



	void M2M(CellClass* const FRestrict ParentCell,
					 const CellClass*const FRestrict *const FRestrict ChildCells,
					 const int TreeLevel)
	{
		// apply Sy
		FBlas::scal(nnodes*2, FReal(0.), ParentCell->getMultipole());
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				AbstractBaseClass::Interpolator->applyM2M(ChildIndex,
																									ChildCells[ChildIndex]->getMultipole(),
																									ParentCell->getMultipole());
	}



	void M2L(CellClass* const FRestrict TargetCell,
					 const CellClass* SourceCells[343],
					 const int NumSourceCells,
					 const int TreeLevel)
	{
		// permute and copy multipole expansion
		memset(countExp, 0, sizeof(int) * 343);
		for (unsigned int idx=0; idx<343; ++idx) {
			if (SourceCells[idx]) {
				const unsigned int pidx = SymHandler->pindices[idx];
				const unsigned int count = (countExp[pidx])++;
				FReal *const mul = Mul[pidx] + count*nnodes;
				const unsigned int *const pvec = SymHandler->pvectors[idx];
				const FReal *const MultiExp = SourceCells[idx]->getMultipole();
				for (unsigned int n=0; n<nnodes; ++n)
					mul[pvec[n]] = MultiExp[n];
			}
		}

		// multiply (mat-mat-mul)
		FReal Compressed [nnodes * 24];
		const FReal scale = AbstractBaseClass::MatrixKernel->getScaleFactor(AbstractBaseClass::BoxWidth, TreeLevel);
		for (unsigned int pidx=0; pidx<343; ++pidx) {
			const unsigned int count = countExp[pidx];
			if (count) {
				const unsigned int rank = SymHandler->getLowRank(TreeLevel, pidx);
				// rank * count * (2*nnodes-1) flops
				FBlas::gemtm(nnodes, rank, count, FReal(1.),
										 const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx))+rank*nnodes,
										 nnodes, Mul[pidx], nnodes, Compressed, rank);
				// nnodes *count * (2*rank-1) flops
				FBlas::gemm( nnodes, rank, count, scale,
										 const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx)),
										 nnodes, Compressed, rank, Loc[pidx], nnodes);
			}
		}

		// permute and add contribution to local expansions
		FReal *const LocalExpansion = TargetCell->getLocal();
		memset(countExp, 0, sizeof(int) * 343);
		for (unsigned int idx=0; idx<343; ++idx) {
			if (SourceCells[idx]) {
				const unsigned int pidx = SymHandler->pindices[idx];
				const unsigned int count = (countExp[pidx])++;
				const FReal *const loc = Loc[pidx] + count*nnodes;
				const unsigned int *const pvec = SymHandler->pvectors[idx];
				// nnodes flops
				for (unsigned int n=0; n<nnodes; ++n)
					LocalExpansion[n] += loc[pvec[n]];
			}
		}

	}



	/*
	void M2L(CellClass* const FRestrict TargetCell,
					 const CellClass* SourceCells[343],
					 const int NumSourceCells,
					 const int TreeLevel)
	{
		// permute and copy multipole expansion
		memset(countExp, 0, sizeof(int) * 343);
		for (unsigned int idx=0; idx<343; ++idx) {
			if (SourceCells[idx]) {
				const unsigned int pidx = SymHandler->pindices[idx];
				const unsigned int count = (countExp[pidx])++;
				const FReal *const MultiExp = SourceCells[idx]->getMultipole();
				for (unsigned int n=0; n<nnodes; ++n)
					Mul[pidx][count*nnodes + SymHandler->pvectors[idx][n]] = MultiExp[n];
			}
		}

		// multiply (mat-mat-mul)
		FReal Compressed [nnodes * 30];
		const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(FMath::pow(2, TreeLevel)));
		const FReal scale(AbstractBaseClass::MatrixKernel->getScaleFactor(CellWidth));
		for (unsigned int pidx=0; pidx<343; ++pidx) {
			const unsigned int count = countExp[pidx];
			if (count) {
				const unsigned int rank = SymHandler->LowRank[pidx];
				FBlas::gemtm(nnodes, rank, count, FReal(1.),
										 SymHandler->K[pidx]+rank*nnodes, nnodes, Mul[pidx], nnodes, Compressed, rank);
				FBlas::gemm( nnodes, rank, count, scale,
										 SymHandler->K[pidx], nnodes, Compressed, rank, Loc[pidx], nnodes);
			}
		}

		// permute and add contribution to local expansions
		FReal *const LocalExpansion = TargetCell->getLocal();
		memset(countExp, 0, sizeof(int) * 343);
		for (unsigned int idx=0; idx<343; ++idx) {
			if (SourceCells[idx]) {
				const unsigned int pidx = SymHandler->pindices[idx];
				const unsigned int count = (countExp[pidx])++;
				for (unsigned int n=0; n<nnodes; ++n)
					LocalExpansion[n] += Loc[pidx][count*nnodes + SymHandler->pvectors[idx][n]];
			}
		}

	}
	*/

	/*
		void M2L(CellClass* const FRestrict TargetCell,
		const CellClass* SourceCells[343],
		const int NumSourceCells,
		const int TreeLevel)
		{
		FReal *const LocalExpansion = TargetCell->getLocal();
		const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
		const FReal scale(MatrixKernel->getScaleFactor(CellWidth));

		FReal PermLocalExp[nnodes];
		FReal PermMultiExp[nnodes];
		FReal Compressed[nnodes];
		for (int i=-3; i<=3; ++i) {
		for (int j=-3; j<=3; ++j) {
		for (int k=-3; k<=3; ++k) {
					
		const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
					
		if (SourceCells[idx]) {
		const FReal *const MultiExp = SourceCells[idx]->getMultipole();

		// permute
		for (unsigned int n=0; n<nnodes; ++n) PermMultiExp[pvectors[idx][n]] = MultiExp[n];

		// mat-vec-mult (svd)
		assert(K[pindices[idx]]!=NULL);
		const int rank = LowRank[pindices[idx]];
		FBlas::gemtv(nnodes, rank, FReal(1.), K[pindices[idx]]+rank*nnodes, PermMultiExp, Compressed);
		FBlas::gemv( nnodes, rank, scale, K[pindices[idx]], Compressed, PermLocalExp);

		// permute
		for (unsigned int n=0; n<nnodes; ++n) LocalExpansion[n] += PermLocalExp[pvectors[idx][n]];
		}

		}
		}
		}

		}
	*/



	void L2L(const CellClass* const FRestrict ParentCell,
					 CellClass* FRestrict *const FRestrict ChildCells,
					 const int TreeLevel)
	{
		// apply Sx
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				AbstractBaseClass::Interpolator->applyL2L(ChildIndex,
																									ParentCell->getLocal(),
																									ChildCells[ChildIndex]->getLocal());
	}



	void L2P(const CellClass* const LeafCell,
					 ContainerClass* const TargetParticles)
	{
		const FPoint LeafCellCenter(getLeafCellCenter(LeafCell->getCoordinate()));

//		// a) apply Sx
//		AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter,
//																							AbstractBaseClass::BoxWidthLeaf,
//																							LeafCell->getLocal(),
//																							TargetParticles);
//		// b) apply Px (grad Sx)
//		AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter,
//																											AbstractBaseClass::BoxWidthLeaf,
//																											LeafCell->getLocal(),
//																											TargetParticles);

		// c) apply Sx and Px (grad Sx)
		AbstractBaseClass::Interpolator->applyL2PTotal(LeafCellCenter,
																									 AbstractBaseClass::BoxWidthLeaf,
																									 LeafCell->getLocal(),
																									 TargetParticles);
	}


};








#endif //FCHEBSYMKERNELS_HPP

// [--END--]
