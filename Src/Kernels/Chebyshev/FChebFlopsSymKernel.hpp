#ifndef FCHEBFLOPSSYMKERNEL_HPP
#define FCHEBFLOPSSYMKERNEL_HPP
// [--License--]

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FSmartPointer.hpp"

#include "../../Components/FAbstractKernels.hpp"

#include "./FChebInterpolator.hpp"
#include "./FChebSymmetries.hpp"

class FTreeCoordinate;



/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebFlopsSymKernel
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
template <class ParticleClass, class CellClass,	class ContainerClass, class MatrixKernelClass, int ORDER>
class FChebFlopsSymKernel
	: public FAbstractKernels<ParticleClass, CellClass, ContainerClass>
{
  enum {nnodes = TensorTraits<ORDER>::nnodes};

	/// Handler to deal with all symmetries: Stores permutation indices and
	/// vectors to reduce 343 different interactions to 16 only.
	struct SymmetryHandler;

	/// Needed for handling all symmetries
	const FSmartPointer<MatrixKernelClass,FSmartPointerMemory> MatrixKernel;
	const FSmartPointer<SymmetryHandler,  FSmartPointerMemory> SymHandler;

	// count permuted local and multipole expansions
	unsigned int* countExp;


	unsigned long long flopsP2M, flopsM2M, flopsM2L, flopsL2L, flopsL2P, flopsP2P;


	// start flop counters 
	unsigned int countFlopsP2MorL2P(const unsigned int N) const
	{	return     N        * (nnodes*(24*ORDER-26) - 1); }
	unsigned int countFlopsM2MorL2L() const
	{ return     8*nnodes * (nnodes*2             - 1); }
	unsigned int countFlopsL2PGradient(const unsigned int N) const
	{ return 3 * N        * (nnodes*(17*ORDER-21) - 1); }
	unsigned int countFlopsL2PTotal(const unsigned int N) const
	{ return     N        * (nnodes*2             -1 ) + countFlopsL2PGradient(N); }
	unsigned int countFlopsM2L(const unsigned int nexp, const unsigned int rank) const
	{ return nexp * (4*nnodes*rank - rank - nnodes); }
	unsigned int countFlopsP2P() const
	{ return 21; }
	unsigned int countFlopsP2Pmutual() const
	{ return 26; }
	// end flop counters

public:
	/**
	 * The constructor initializes all constant attributes and it reads the
	 * precomputed and compressed M2L operators from a binary file (an
	 * runtime_error is thrown if the required file is not valid).
	 */
	FChebFlopsSymKernel(const int inTreeHeight,
								 const FPoint& inBoxCenter,
								 const FReal inBoxWidth,
								 const FReal Epsilon)
		:	MatrixKernel(new MatrixKernelClass()),
			SymHandler(new SymmetryHandler(MatrixKernel.getPtr(), Epsilon)),
			flopsP2M(0), flopsM2M(0), flopsM2L(0), flopsL2L(0), flopsL2P(0), flopsP2P(0)
	{	countExp = new unsigned int [343]; }
	


	/** Copy constructor */
	FChebFlopsSymKernel(const FChebFlopsSymKernel& other)
		: SymHandler(other.SymHandler),
			flopsP2M(0), flopsM2M(0), flopsM2L(0), flopsL2L(0), flopsL2P(0), flopsP2P(0)
	{	countExp = new unsigned int [343]; }



	/** Destructor */
	~FChebFlopsSymKernel()
	{
		delete [] countExp;

		std::cout << "\n==================================================" 
							<< "\n- Flops for P2M = " << flopsP2M 
							<< "\n- Flops for M2M = " << flopsM2M
							<< "\n- Flops for M2L = " << flopsM2L
							<< "\n- Flops for L2L = " << flopsL2L
							<< "\n- Flops for L2P = " << flopsL2P
							<< "\n- Flops for P2P = " << flopsP2P
							<< "\n==================================================\n"
							<< std::endl;

	}
	
	
	
	void P2M(CellClass* const /* not needed */, const ContainerClass* const SourceParticles)
	{
		flopsP2M += countFlopsP2MorL2P(SourceParticles->getSize());
	}



	void M2M(CellClass* const FRestrict /* not needed */,
					 const CellClass*const FRestrict *const FRestrict ChildCells,
					 const int /* not needed */)
	{
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])	flopsM2M += countFlopsM2MorL2L();
	}




	void M2L(CellClass* const FRestrict /* not needed */,
					 const CellClass* SourceCells[343],
					 const int /* not needed */,
					 const int /* not needed*/)
	{
		// count how ofter each of the 16 interactions is used
		memset(countExp, 0, sizeof(int) * 343);
		for (unsigned int idx=0; idx<343; ++idx)
			if (SourceCells[idx])	countExp[SymHandler->pindices[idx]]++;
		// multiply (mat-mat-mul)
		for (unsigned int pidx=0; pidx<343; ++pidx)
			if (countExp[pidx]) flopsM2L += countFlopsM2L(countExp[pidx], SymHandler->LowRank[pidx]) + countExp[pidx]*nnodes;
	}


	void L2L(const CellClass* const FRestrict /* not needed */,
					 CellClass* FRestrict *const FRestrict ChildCells,
					 const int /* not needed */)
	{
		for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
			if (ChildCells[ChildIndex])	flopsL2L += countFlopsM2MorL2L() + nnodes;
	}



	void L2P(const CellClass* const /* not needed */,
					 ContainerClass* const TargetParticles)
	{
		//// 1.a) apply Sx
		//flopsL2P += countFlopsP2MorL2P(TargetParticlesParticles->getSize()) + TargetParticles->getSize();
		//// 1.b) apply Px (grad Sx)
		//flopsL2P += countFlopsL2PGradient(TargetParticlesParticles->getSize()) + 3 * TargetParticles->getSize();

		// or

		// 2) apply Sx and Px (grad Sx)
		flopsL2P += countFlopsL2PTotal(TargetParticles->getSize()) + 4 * TargetParticles->getSize();
	}



	void P2P(const FTreeCoordinate& /* LeafCellCoordinate */, // needed for periodic boundary conditions
					 ContainerClass* const FRestrict TargetParticles,
					 const ContainerClass* const FRestrict SourceParticles,
					 ContainerClass* const NeighborSourceParticles[27],
					 const int /* size */)
	{
		if (TargetParticles != SourceParticles) {
			flopsP2P += countFlopsP2P() * TargetParticles->getSize() * SourceParticles->getSize();
			for (unsigned int idx=0; idx<27; ++idx)
				if (NeighborSourceParticles[idx])
					flopsP2P += countFlopsP2P() * TargetParticles->getSize() * NeighborSourceParticles[idx]->getSize();
		} else {
			flopsP2P += countFlopsP2Pmutual() * ((TargetParticles->getSize()*TargetParticles->getSize()+TargetParticles->getSize()) / 2); 
			for (unsigned int idx=0; idx<=13; ++idx)
				if (NeighborSourceParticles[idx])
					flopsP2P += countFlopsP2Pmutual() * TargetParticles->getSize() * NeighborSourceParticles[idx]->getSize();
		}
	}


};










/**
 * Handler to deal with all symmetries: Stores permutation indices and vectors
 * to reduce 343 different interactions to 16 only.
 */
template <class ParticleClass,
					class CellClass,
					class ContainerClass,
					class MatrixKernelClass,
					int ORDER>
struct FChebFlopsSymKernel<ParticleClass, CellClass, ContainerClass, MatrixKernelClass, ORDER>
::SymmetryHandler
{
	// M2L operators
	FReal*    K[343];
	int LowRank[343];
		
	// permutation vectors and permutated indices
	unsigned int pvectors[343][nnodes];
	unsigned int pindices[343];


	// compute rank
	unsigned int getRank(const FReal singular_values[], const double eps)
	{
		FReal nrm2(0.);
		for (unsigned int k=0; k<nnodes; ++k)
			nrm2 += singular_values[k] * singular_values[k];
		
		FReal nrm2k(0.);
		for (unsigned int k=nnodes; k>0; --k) {
			nrm2k += singular_values[k-1] * singular_values[k-1];
			if (nrm2k > eps*eps * nrm2)	return k;
		}
		throw std::runtime_error("rank cannot be larger than nnodes");
		return 0;
	}
	

	/** Constructor */
	SymmetryHandler(const MatrixKernelClass *const MatrixKernel, const double Epsilon)
	{
		// init all 343 item to zero, because effectively only 16 exist
		for (unsigned int t=0; t<343; ++t) {
			K[t] = NULL;
			LowRank[t] = 0;
		}
			
		// set permutation vector and indices
		const FChebSymmetries<ORDER> Symmetries;
		for (int i=-3; i<=3; ++i)
			for (int j=-3; j<=3; ++j)
				for (int k=-3; k<=3; ++k)
					if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
						const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
						pindices[idx] = Symmetries.getPermutationArrayAndIndex(i,j,k, pvectors[idx]);
					}

		// precompute 16 M2L operators
		this->precomputeSVD(MatrixKernel, Epsilon);
	}



	/** Destructor */
	~SymmetryHandler()
	{
		for (unsigned int t=0; t<343; ++t)
			if (K[  t]!=NULL) delete [] K[  t];
	}



private:
	void precomputeSVD(const MatrixKernelClass *const MatrixKernel, const double Epsilon)
	{
		// interpolation points of source (Y) and target (X) cell
		FPoint X[nnodes], Y[nnodes];
		// set roots of target cell (X)
		FChebTensor<ORDER>::setRoots(FPoint(0.,0.,0.), FReal(2.), X);
		// temporary matrix
		FReal* U = new FReal [nnodes*nnodes];

		// needed for the SVD
		const unsigned int LWORK = 2 * (3*nnodes + nnodes);
		FReal *const WORK = new FReal [LWORK];
		FReal *const VT = new FReal [nnodes*nnodes];
		FReal *const S = new FReal [nnodes];
		
		unsigned int counter = 0;
		for (int i=2; i<=3; ++i) {
			for (int j=0; j<=i; ++j) {
				for (int k=0; k<=j; ++k) {

					// assemble matrix
					const FPoint cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
					FChebTensor<ORDER>::setRoots(cy, FReal(2.), Y);
					for (unsigned int n=0; n<nnodes; ++n)
						for (unsigned int m=0; m<nnodes; ++m)
							U[n*nnodes + m] = MatrixKernel->evaluate(X[m], Y[n]);

					// applying weights ////////////////////////////////////////
					FReal weights[nnodes];
					FChebTensor<ORDER>::setRootOfWeights(weights);
					for (unsigned int n=0; n<nnodes; ++n) {
						FBlas::scal(nnodes, weights[n], U + n,  nnodes); // scale rows
						FBlas::scal(nnodes, weights[n], U + n * nnodes); // scale cols
					}
					//////////////////////////////////////////////////////////		

					// truncated singular value decomposition of matrix
					const unsigned int info	= FBlas::gesvd(nnodes, nnodes, U, S, VT, nnodes, LWORK, WORK);
					if (info!=0) throw std::runtime_error("SVD did not converge with " + info);
					const unsigned int rank = this->getRank(S, Epsilon);

					// store 
					const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
					assert(K[idx]==NULL);
					K[idx] = new FReal [2*rank*nnodes];
					LowRank[idx] = rank;
					for (unsigned int r=0; r<rank; ++r)
						FBlas::scal(nnodes, S[r], U + r*nnodes);
					FBlas::copy(rank*nnodes, U,  K[idx]);
					for (unsigned int r=0; r<rank; ++r)
						FBlas::copy(nnodes, VT + r, nnodes, K[idx] + rank*nnodes + r*nnodes, 1);

					// un-weighting ////////////////////////////////////////////
					for (unsigned int n=0; n<nnodes; ++n) {
						FBlas::scal(rank, FReal(1.) / weights[n], K[idx] + n,               nnodes); // scale rows
						FBlas::scal(rank, FReal(1.) / weights[n], K[idx] + rank*nnodes + n, nnodes); // scale rows
					}
					//////////////////////////////////////////////////////////		

					std::cout << "(" << i << "," << j << "," << k << ") " << idx <<
						", low rank = " << rank << std::endl;

					counter++;
				}
			}
		}
		std::cout << "num interactions = " << counter << std::endl;
		delete [] U;
		delete [] WORK;
		delete [] VT;
		delete [] S;
	}

};











#endif //FCHEBSYMKERNELS_HPP

// [--END--]
