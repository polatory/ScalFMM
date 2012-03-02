#ifndef FCHEBSYMKERNELS_HPP
#define FCHEBSYMKERNELS_HPP
// [--License--]

#include "../Utils/FGlobal.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FSmartPointer.hpp"

#include "../Components/FAbstractKernels.hpp"

#include "./FChebInterpolator.hpp"
#include "./FChebSymmetries.hpp"

class FTreeCoordinate;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebSymKernels
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
class FChebSymKernels : public FAbstractKernels<ParticleClass, CellClass, ContainerClass>
{
  enum {nnodes = TensorTraits<ORDER>::nnodes};
	typedef FChebInterpolator<ORDER> InterpolatorClass;
	typedef FChebSymmetries<ORDER> SymmetriesClass;

	/// Needed for P2M, M2M, L2L and L2P operators
	FSmartPointer<InterpolatorClass,FSmartPointerMemory> Interpolator;
	/// Needed for M2L operator
	FSmartPointer<  SymmetriesClass,FSmartPointerMemory> Symmetries;
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

	// M2L operators
	FReal*    K[343];
	int LowRank[343];

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

	void precomputeDense()
	{
		// interpolation points of source (Y) and target (X) cell
		F3DPosition X[nnodes], Y[nnodes];
		// set roots of target cell (X)
		FChebTensor<ORDER>::setRoots(F3DPosition(0.,0.,0.), FReal(2.), X);
		unsigned int counter = 0;
		for (int i=2; i<=3; ++i) {
			for (int j=0; j<=i; ++j) {
				for (int k=0; k<=j; ++k) {
					const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
					K[idx] = new FReal [nnodes*nnodes];
					const F3DPosition cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
					FChebTensor<ORDER>::setRoots(cy, FReal(2.), Y);
					for (unsigned int n=0; n<nnodes; ++n)
						for (unsigned int m=0; m<nnodes; ++m)
							K[idx][n*nnodes + m] = MatrixKernel->evaluate(X[m], Y[n]);
					counter++;
				}
			}
		}
		std::cout << "num interactions = " << counter << std::endl;
	}

	void precomputeSVD()
	{
		// interpolation points of source (Y) and target (X) cell
		F3DPosition X[nnodes], Y[nnodes];
		// set roots of target cell (X)
		FChebTensor<ORDER>::setRoots(F3DPosition(0.,0.,0.), FReal(2.), X);
		// temporary matrix
		FReal* U = new FReal [nnodes*nnodes];

		unsigned int counter = 0;
		for (int i=2; i<=3; ++i) {
			for (int j=0; j<=i; ++j) {
				for (int k=0; k<=j; ++k) {

					// assemble matrix
					const F3DPosition cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
					FChebTensor<ORDER>::setRoots(cy, FReal(2.), Y);
					for (unsigned int n=0; n<nnodes; ++n)
						for (unsigned int m=0; m<nnodes; ++m)
							U[n*nnodes + m] = MatrixKernel->evaluate(X[m], Y[n]);

					// truncated singular value decomposition of matrix
					const unsigned int LWORK = 2 * (3*nnodes + nnodes);
					FReal *const WORK = new FReal [LWORK];
					FReal *const VT = new FReal [nnodes*nnodes];
					FReal *const S = new FReal [nnodes];
					const unsigned int info	= FBlas::gesvd(nnodes, nnodes, U, S, VT, nnodes, LWORK, WORK);
					if (info!=0) throw std::runtime_error("SVD did not converge with " + info);
					const unsigned int rank = getRank<ORDER>(S, Epsilon);

					// store 
					const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
					K[idx] = new FReal [2*rank*nnodes];
					LowRank[idx] = rank;
					for (unsigned int r=0; r<rank; ++r)
						FBlas::scal(nnodes, S[r], U + r*nnodes);
					FBlas::copy(rank*nnodes, U,  K[idx]);
					for (unsigned int r=0; r<rank; ++r)
						FBlas::copy(nnodes, VT + r, nnodes, K[idx] + rank*nnodes + r*nnodes, 1);

					std::cout << "(" << i << "," << j << "," << k << ") low rank = " << rank << std::endl;

					counter++;
				}
			}
		}
		std::cout << "num interactions = " << counter << std::endl;
		delete [] U;
	}
	
public:
	/**
	 * The constructor initializes all constant attributes and it reads the
	 * precomputed and compressed M2L operators from a binary file (an
	 * runtime_error is thrown if the required file is not valid).
	 */
	FChebSymKernels(const int inTreeHeight,
									const F3DPosition& inBoxCenter,
									const FReal inBoxWidth,
									const FReal inEpsilon)
		: Interpolator(new InterpolatorClass()),
			Symmetries(new SymmetriesClass()),
			MatrixKernel(new MatrixKernelClass()),
			TreeHeight(inTreeHeight),
			BoxCorner(inBoxCenter - inBoxWidth / FReal(2.)),
			BoxWidth(inBoxWidth),
			BoxWidthLeaf(BoxWidth / FReal(FMath::pow(2, inTreeHeight - 1))),
			Epsilon(inEpsilon)
	{
		// precompute M2L operators
		for (unsigned int t=0; t<343; ++t) {K[t] = NULL; LowRank[t] = -1;}
		
		//this -> precomputeDense();
		this -> precomputeSVD();
		

	}

	~FChebSymKernels()
	{
		for (unsigned int i=0; i<343; ++i) if (K[i]!=NULL) delete [] K[i];
	}

	void P2M(CellClass* const LeafCell,
					 const ContainerClass* const SourceParticles)
	{
		// apply Sy
		const F3DPosition LeafCellCenter(getLeafCellCenter(LeafCell->getCoordinate()));
		Interpolator->applyP2M(LeafCellCenter,
													 BoxWidthLeaf,
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
				Interpolator->applyM2M(ChildIndex,
															 ChildCells[ChildIndex]->getMultipole(),
															 ParentCell->getMultipole());
	}



	void M2L(CellClass* const FRestrict TargetCell,
					 const CellClass* SourceCells[343],
					 const int NumSourceCells,
					 const int TreeLevel)
	{
		FReal *const LocalExpansion = TargetCell->getLocal();
		const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
		const FReal scale(MatrixKernel->getScaleFactor(CellWidth));

		unsigned int permutations[nnodes];
		FReal PermLocalExp[nnodes];
		FReal PermMultiExp[nnodes];
		FReal Compressed[nnodes];
		for (int i=-3; i<=3; ++i) {
			for (int j=-3; j<=3; ++j) {
				for (int k=-3; k<=3; ++k) {
					const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
					if (SourceCells[idx]) {
						const FReal *const MultiExp = SourceCells[idx]->getMultipole();
						const unsigned int pidx = Symmetries->getPermutationArrayAndIndex(i,j,k, permutations);
						if (K[pidx]==NULL) std::cout << pidx << "does not exist" << std::endl;
						for (unsigned int n=0; n<nnodes; ++n) PermMultiExp[permutations[n]] = MultiExp[n];

						//// dense
						//FBlas::gemv(nnodes, nnodes, scale, K[pidx], PermMultiExp, PermLocalExp);
						
						// svd
						const int rank = LowRank[pidx];
						FBlas::gemtv(nnodes, rank, FReal(1.), K[pidx]+rank*nnodes, PermMultiExp, Compressed);
						FBlas::gemv( nnodes, rank, scale, K[pidx], Compressed, PermLocalExp);

						for (unsigned int n=0; n<nnodes; ++n) LocalExpansion[n] += PermLocalExp[permutations[n]];
					}
				}
			}
		}

	}



	void L2L(const CellClass* const FRestrict ParentCell,
					 CellClass* FRestrict *const FRestrict ChildCells,
					 const int TreeLevel)
	{
		// apply Sx
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])
				Interpolator->applyL2L(ChildIndex,
															 ParentCell->getLocal(),
															 ChildCells[ChildIndex]->getLocal());
	}

	void L2P(const CellClass* const LeafCell,
					 ContainerClass* const TargetParticles)
	{
		// a) apply Sx
		const F3DPosition LeafCellCenter(getLeafCellCenter(LeafCell->getCoordinate()));
		Interpolator->applyL2P(LeafCellCenter,
													 BoxWidthLeaf,
													 LeafCell->getLocal(),
													 TargetParticles);
		// b) apply Px (grad Sx)
		Interpolator->applyL2PGradient(LeafCellCenter,
																	 BoxWidthLeaf,
																	 LeafCell->getLocal(),
																	 TargetParticles);
	}


	void P2P(const FTreeCoordinate& LeafCellCoordinate,
					 ContainerClass* const FRestrict TargetParticles,
					 const ContainerClass* const FRestrict SourceParticles,
					 ContainerClass* const NeighborSourceParticles[27],
					 const int )
	{
		// loop: target particles
		typename ContainerClass::BasicIterator iTargets(*TargetParticles);
		while (iTargets.hasNotFinished()) {
			ParticleClass& Target = iTargets.data();
			const FReal wt = Target.getPhysicalValue();

			{ // loop: source particles (target leaf cell == source leaf cell)
				typename ContainerClass::ConstBasicIterator iSources(*SourceParticles);
				while (iSources.hasNotFinished()) {
					const ParticleClass& Source = iSources.data();
					// only if target and source are not identical
					if (&Target != &Source) {
						const FReal one_over_r = MatrixKernel->evaluate(Target.getPosition(), Source.getPosition());
						const FReal ws = Source.getPhysicalValue();
						// potential
						Target.incPotential(one_over_r * ws);
						// force
						F3DPosition force(Target.getPosition() - Source.getPosition());
						force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
						Target.incForces(force.getX(), force.getY(), force.getZ());
					}
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
							const FReal one_over_r = MatrixKernel->evaluate(Target.getPosition(), Source.getPosition());
							const FReal ws = Source.getPhysicalValue();
							// potential
							Target.incPotential(one_over_r * ws);
							F3DPosition force(Target.getPosition() - Source.getPosition());
							force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
							// force
							Target.incForces(force.getX(), force.getY(), force.getZ());
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


#endif //FCHEBSYMKERNELS_HPP

// [--END--]




//			for (int i=-3; i<=3; ++i) {
//				for (int j=-3; j<=3; ++j) {
//					for (int k=-3; k<=3; ++k) {
//						if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
//
//							const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
//							K[idx] = new FReal [nnodes*nnodes];
//							const F3DPosition cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
//							FChebTensor<ORDER>::setRoots(cy, FReal(2.), Y);
//							for (unsigned int n=0; n<nnodes; ++n)
//								for (unsigned int m=0; m<nnodes; ++m)
//									K[idx][n*nnodes + m] = MatrixKernel->evaluate(X[m], Y[n]);
//							counter++;
//
//						}
//					}
//				}
//			}



//		for (unsigned int idx=0; idx<343; ++idx)
//			if (SourceCells[idx])
//				FBlas::gemva(nnodes, nnodes, scale, K[idx],
//										 const_cast<CellClass*>(SourceCells[idx])->getMultipole(), LocalExpansion);
