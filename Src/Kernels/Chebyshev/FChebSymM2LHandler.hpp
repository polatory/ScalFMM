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
#ifndef FCHEBSYMM2LHANDLER_HPP
#define FCHEBSYMM2LHANDLER_HPP

#include <climits>

#include "../../Utils/FBlas.hpp"

#include "./FChebTensor.hpp"
#include "../Interpolation/FInterpSymmetries.hpp"
#include "./FChebM2LHandler.hpp"

/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */


/*!  Choose either \a FULLY_PIVOTED_ACASVD or \a PARTIALLY_PIVOTED_ACASVD or
	\a ONLY_SVD.
*/
//#define ONLY_SVD
//#define FULLY_PIVOTED_ACASVD
#define PARTIALLY_PIVOTED_ACASVD




/*!  The fully pivoted adaptive cross approximation (fACA) compresses a
	far-field interaction as \f$K\sim UV^\top\f$. The fACA requires all entries
	to be computed beforehand, then the compression follows in
	\f$\mathcal{O}(2\ell^3k)\f$ operations based on the required accuracy
	\f$\varepsilon\f$. The matrix K will be destroyed as a result.

	@param[in] K far-field to be approximated
	@param[in] nx number of rows
	@param[in] ny number of cols
	@param[in] eps prescribed accuracy
	@param[out] U matrix containing \a k column vectors
	@param[out] V matrix containing \a k row vectors
	@param[out] k final low-rank depends on prescribed accuracy \a eps
*/
void fACA(FReal *const K,
					const unsigned int nx, const unsigned int ny,
					const double eps, FReal* &U, FReal* &V, unsigned int &k)
{
	// control vectors (true if not used, false if used)
	bool *const r = new bool[nx];
	bool *const c = new bool[ny];
	for (unsigned int i=0; i<nx; ++i) r[i] = true;
	for (unsigned int j=0; j<ny; ++j) c[j] = true;

	// compute Frobenius norm of original Matrix K
	FReal norm2K = 0;
	for (unsigned int j=0; j<ny; ++j) {
		const FReal *const colK = K + j*nx;
		norm2K += FBlas::scpr(nx, colK, colK);
	}

	// initialize rank k and UV'
	k = 0;
	const int maxk = (nx + ny) / 2;
	U = new FReal[nx * maxk];
	V = new FReal[ny * maxk];
	FBlas::setzero(nx*maxk, U);
	FBlas::setzero(ny*maxk, V);
	FReal norm2R;

	////////////////////////////////////////////////
	// start fully pivoted ACA
	do {
		
		// find max(K) and argmax(K)
		FReal maxK = 0.;
		int pi=0, pj=0;
		for (unsigned int j=0; j<ny; ++j)
			if (c[j]) {
				const FReal *const colK = K + j*nx;
				for (unsigned int i=0; i<nx; ++i)
					if (r[i] && maxK < FMath::Abs(colK[i])) {
						maxK = FMath::Abs(colK[i]);
						pi = i; 
						pj = j;
					}
			}

		// copy pivot cross into U and V
		FReal *const colU = U + k*nx;
		FReal *const colV = V + k*ny;
		const FReal pivot = K[pj*nx + pi];
		for (unsigned int i=0; i<nx; ++i) if (r[i]) colU[i] = K[pj*nx + i];
		for (unsigned int j=0; j<ny; ++j) if (c[j]) colV[j] = K[j *nx + pi] / pivot;
		
		// dont use these cols and rows anymore
		c[pj] = false;
		r[pi] = false;
		
		// subtract k-th outer product from K
		for (unsigned int j=0; j<ny; ++j)
			if (c[j]) {
				FReal *const colK = K + j*nx;
				FBlas::axpy(nx, FReal(-1. * colV[j]), colU, colK);
			}

		// compute Frobenius norm of updated K
		norm2R = 0.;
		for (unsigned int j=0; j<ny; ++j)
			if (c[j]) {
				const FReal *const colK = K + j*nx;
				norm2R += FBlas::scpr(nx, colK, colK);
			}

		// increment rank k
		k++;
		
	} while (norm2R > eps*eps * norm2K);
	////////////////////////////////////////////////

	delete [] r;
	delete [] c;
}











/*!  The partially pivoted adaptive cross approximation (pACA) compresses a
	far-field interaction as \f$K\sim UV^\top\f$. The pACA computes the matrix
	entries on the fly, as they are needed. The compression follows in
	\f$\mathcal{O}(2\ell^3k)\f$ operations based on the required accuracy
	\f$\varepsilon\f$. The matrix K will be destroyed as a result.

	@tparam ComputerType the functor type which allows to compute matrix entries
	
	@param[in] Computer the entry-computer functor
	@param[in] eps prescribed accuracy
	@param[in] nx number of rows
	@param[in] ny number of cols
	@param[out] U matrix containing \a k column vectors
	@param[out] V matrix containing \a k row vectors
	@param[out] k final low-rank depends on prescribed accuracy \a eps
*/
template <typename ComputerType>
void pACA(const ComputerType& Computer,
					const unsigned int nx, const unsigned int ny,
					const FReal eps, FReal* &U, FReal* &V, unsigned int &k)
{
	// control vectors (true if not used, false if used)
	bool *const r = new bool[nx];
	bool *const c = new bool[ny];
	for (unsigned int i=0; i<nx; ++i) r[i] = true;
	for (unsigned int j=0; j<ny; ++j) c[j] = true;
	
	// initialize rank k and UV'
	k = 0;
	const FReal eps2 = eps * eps;
	const int maxk = (nx + ny) / 2;
	U = new FReal[nx * maxk];
	V = new FReal[ny * maxk];
	
	// initialize norm
	FReal norm2S(0.);
	FReal norm2uv(0.);
	
	////////////////////////////////////////////////
	// start partially pivoted ACA
	unsigned int J = 0, I = 0;
	
	do {
		FReal *const colU = U + nx*k;
		FReal *const colV = V + ny*k;
		
		////////////////////////////////////////////
		// compute row I and its residual
		Computer(I, I+1, 0, ny, colV);
		r[I] = false;
		for (unsigned int l=0; l<k; ++l) {
			FReal *const u = U + nx*l;
			FReal *const v = V + ny*l;
			FBlas::axpy(ny, FReal(-1. * u[I]), v, colV);
		}
		
		// find max of residual and argmax
		FReal maxval = 0.;
		for (unsigned int j=0; j<ny; ++j) {
			const FReal abs_val = FMath::Abs(colV[j]);
			if (c[j] && maxval < abs_val) {
				maxval = abs_val;
				J = j;
			}
		}
		// find pivot and scale column of V
		const FReal pivot = FReal(1.) / colV[J];
		FBlas::scal(ny, pivot, colV);
		
		////////////////////////////////////////////
		// compute col J and its residual
		Computer(0, nx, J, J+1, colU);
		c[J] = false;
		for (unsigned int l=0; l<k; ++l) {
			FReal *const u = U + nx*l;
			FReal *const v = V + ny*l;
			FBlas::axpy(nx, FReal(-1. * v[J]), u, colU);
		}
		
		// find max of residual and argmax
		maxval = 0.;
		for (unsigned int i=0; i<nx; ++i) {
			const FReal abs_val = FMath::Abs(colU[i]);
			if (r[i] && maxval < abs_val) {
				maxval = abs_val;
				I = i;
			}
		}
		
		////////////////////////////////////////////
			// increment Frobenius norm: |Sk|^2 += |uk|^2 |vk|^2 + 2 sumj ukuj vjvk
		FReal normuuvv(0.);
		for (unsigned int l=0; l<k; ++l)
			normuuvv += FBlas::scpr(nx, colU, U + nx*l) * FBlas::scpr(ny, V + ny*l, colV);
		norm2uv = FBlas::scpr(nx, colU, colU) * FBlas::scpr(ny, colV, colV);
		norm2S += norm2uv + 2*normuuvv;
		
		////////////////////////////////////////////
		// increment low-rank
		k++;

	} while (norm2uv > eps2 * norm2S);
	
	delete [] r;
	delete [] c;
}



/*!  Precomputes the 16 far-field interactions (due to symmetries in their
  arrangement all 316 far-field interactions can be represented by
  permutations of the 16 we compute in this function). Depending on whether
  FACASVD is defined or not, either ACA+SVD or only SVD is used to compress
  them. */
template <int ORDER, typename MatrixKernelClass>
static void precompute(const MatrixKernelClass *const MatrixKernel, const FReal CellWidth,
											 const FReal Epsilon, FReal* K[343], int LowRank[343])
{
  //	std::cout << "\nComputing 16 far-field interactions (l=" << ORDER << ", eps=" << Epsilon
  //						<< ") for cells of width w = " << CellWidth << std::endl;

	static const unsigned int nnodes = ORDER*ORDER*ORDER;

	// interpolation points of source (Y) and target (X) cell
	FPoint X[nnodes], Y[nnodes];
	// set roots of target cell (X)
	FChebTensor<ORDER>::setRoots(FPoint(0.,0.,0.), CellWidth, X);
	// temporary matrix
	FReal* U = new FReal [nnodes*nnodes];

	// needed for the SVD
	unsigned int INFO;
	const unsigned int LWORK = 2 * (3*nnodes + nnodes);
	FReal *const WORK = new FReal [LWORK];
	FReal *const VT = new FReal [nnodes*nnodes];
	FReal *const S = new FReal [nnodes];


	// initialize timer
	FTic time;
	double overall_time(0.);
	double elapsed_time(0.);

	// initialize rank counter
	unsigned int overall_rank = 0;

	unsigned int counter = 0;
	for (int i=2; i<=3; ++i) {
		for (int j=0; j<=i; ++j) {
			for (int k=0; k<=j; ++k) {

				// assemble matrix and apply weighting matrices
				const FPoint cy(CellWidth*FReal(i), CellWidth*FReal(j), CellWidth*FReal(k));
				FChebTensor<ORDER>::setRoots(cy, CellWidth, Y);
				FReal weights[nnodes];
				FChebTensor<ORDER>::setRootOfWeights(weights);

				// now the entry-computer is responsible for weighting the matrix entries
				EntryComputer<MatrixKernelClass> Computer(nnodes, X, nnodes, Y, weights);

				// start timer
				time.tic();

#if (defined ONLY_SVD || defined FULLY_PIVOTED_ACASVD)
				Computer(0, nnodes, 0, nnodes, U);
#endif
				/*
				// applying weights ////////////////////////////////////////
				FReal weights[nnodes];
				FChebTensor<ORDER>::setRootOfWeights(weights);
				for (unsigned int n=0; n<nnodes; ++n) {
					FBlas::scal(nnodes, weights[n], U + n,  nnodes); // scale rows
					FBlas::scal(nnodes, weights[n], U + n * nnodes); // scale cols
				}
				*/

				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				// ALL PREPROC FLAGS ARE SET ON TOP OF THIS FILE !!! /////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////



				//////////////////////////////////////////////////////////////
#if (defined FULLY_PIVOTED_ACASVD || defined PARTIALLY_PIVOTED_ACASVD) ////////////
				FReal *UU, *VV;
				unsigned int rank;

#ifdef FULLY_PIVOTED_ACASVD
				fACA(U,        nnodes, nnodes, Epsilon, UU, VV, rank);
#else
				pACA(Computer, nnodes, nnodes, Epsilon, UU, VV, rank);
#endif 

				// QR decomposition
				FReal* phi = new FReal [rank*rank];
				{
					// QR of U and V
					FReal* tauU = new FReal [rank];
					INFO = FBlas::geqrf(nnodes, rank, UU, tauU, LWORK, WORK);
					assert(INFO==0);
					FReal* tauV = new FReal [rank];
					INFO = FBlas::geqrf(nnodes, rank, VV, tauV, LWORK, WORK);
					assert(INFO==0);
					// phi = Ru Rv'
					FReal* rU = new FReal [2 * rank*rank];
					FReal* rV = rU + rank*rank;
					FBlas::setzero(2 * rank*rank, rU);
					for (unsigned int l=0; l<rank; ++l) {
						FBlas::copy(l+1, UU + l*nnodes, rU + l*rank);
						FBlas::copy(l+1, VV + l*nnodes, rV + l*rank);
					}
					FBlas::gemmt(rank, rank, rank, FReal(1.), rU, rank, rV, rank, phi, rank);
					delete [] rU;
					// get Qu and Qv
					INFO = FBlas::orgqr(nnodes, rank, UU, tauU, LWORK, WORK);
					assert(INFO==0);
					INFO = FBlas::orgqr(nnodes, rank, VV, tauV, LWORK, WORK);
					assert(INFO==0);
					delete [] tauU;
					delete [] tauV;
				}
				
				const unsigned int aca_rank = rank;

				// SVD
				{
					INFO = FBlas::gesvd(aca_rank, aca_rank, phi, S, VT, aca_rank, LWORK, WORK);
					if (INFO!=0) throw std::runtime_error("SVD did not converge with " + INFO);
					rank = getRank(S, aca_rank, Epsilon);
				}					
				
				const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);

				// store
				{
					// allocate
					assert(K[idx]==NULL);
					K[idx] = new FReal [2*rank*nnodes];
					
					// set low rank
					LowRank[idx] = rank;
					
					// (U Sigma)
					for (unsigned int r=0; r<rank; ++r)
						FBlas::scal(aca_rank, S[r], phi + r*aca_rank);

					// Qu (U Sigma) 
					FBlas::gemm(nnodes, aca_rank, rank, FReal(1.), UU, nnodes, phi, aca_rank, K[idx], nnodes);
					delete [] phi;

					// Vt -> V and then Qu V
					FReal *const V = new FReal [aca_rank * rank];
					for (unsigned int r=0; r<rank; ++r)
						FBlas::copy(aca_rank, VT + r, aca_rank, V + r*aca_rank, 1);
					FBlas::gemm(nnodes, aca_rank, rank, FReal(1.), VV, nnodes, V, aca_rank, K[idx] + rank*nnodes, nnodes);
					delete [] V;
				}

				//// store recompressed UV
				//const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
				//assert(K[idx]==NULL);
				//K[idx] = new FReal [2*rank*nnodes];
				//LowRank[idx] = rank;
				//FBlas::copy(rank*nnodes, UU,  K[idx]);
				//FBlas::copy(rank*nnodes, VV,  K[idx] + rank*nnodes);
			
				delete [] UU;
				delete [] VV;

				elapsed_time = time.tacAndElapsed(); 
				overall_time += elapsed_time;
				overall_rank += rank;
				// std::cout << "(" << i << "," << j << "," << k << ") " << idx <<
				// 	", low rank = " << rank << " (" << aca_rank << ") in " << elapsed_time << "s" << std::endl;

				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				// ALL PREPROC FLAGS ARE SET ON TOP OF THIS FILE !!! /////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////

#elif defined ONLY_SVD
				// truncated singular value decomposition of matrix
				INFO = FBlas::gesvd(nnodes, nnodes, U, S, VT, nnodes, LWORK, WORK);
				if (INFO!=0) throw std::runtime_error("SVD did not converge with " + INFO);
				const unsigned int rank = getRank<ORDER>(S, Epsilon);
				
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

				elapsed_time = time.tacAndElapsed(); 
				overall_time += elapsed_time;
				overall_rank += rank;
				//				std::cout << "(" << i << "," << j << "," << k << ") " << idx <<
				//	", low rank = " << rank << " in " << elapsed_time << "s" << std::endl;
#else
#error Either fully-, partially pivoted ACA or only SVD must be defined!
#endif ///////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////

				
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				// ALL PREPROC FLAGS ARE SET ON TOP OF THIS FILE !!! /////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////


				// un-weighting ////////////////////////////////////////////
				for (unsigned int n=0; n<nnodes; ++n) {
					FBlas::scal(rank, FReal(1.) / weights[n], K[idx] + n,               nnodes); // scale rows
					FBlas::scal(rank, FReal(1.) / weights[n], K[idx] + rank*nnodes + n, nnodes); // scale rows
				}
				//////////////////////////////////////////////////////////		

				counter++;
			}
		}
	}
		std::cout << "The approximation of the " << counter
              << " far-field interactions (overall rank " << overall_rank 
              << " / " << 16*nnodes 
              << " , sizeM2L= " << 2*overall_rank*nnodes*sizeof(FReal) << ""
              << " / " << 16*nnodes*nnodes*sizeof(FReal) << " B"
              << ") took " << overall_time << "s\n" << std::endl;
	delete [] U;
	delete [] WORK;
	delete [] VT;
	delete [] S;
}









/*!  \class SymmetryHandler 

	\brief Deals with all the symmetries in the arrangement of the far-field interactions

	Stores permutation indices and permutation vectors to reduce 316 (7^3-3^3)
  different far-field interactions to 16 only. We use the number 343 (7^3)
  because it allows us to use to associate the far-field interactions based on
  the index \f$t = 7^2(i+3) + 7(j+3) + (k+3)\f$ where \f$(i,j,k)\f$ denotes
  the relative position of the source cell to the target cell. */
template <int ORDER, KERNEL_FUNCTION_TYPE TYPE> class SymmetryHandler;

/*! Specialization for homogeneous kernel functions */
template <int ORDER>
class SymmetryHandler<ORDER, HOMOGENEOUS>
{
  static const unsigned int nnodes = ORDER*ORDER*ORDER;

	// M2L operators
	FReal*    K[343];
	int LowRank[343];

public:
	
	// permutation vectors and permutated indices
	unsigned int pvectors[343][nnodes];
	unsigned int pindices[343];


	/** Constructor: with 16 small SVDs */
	template <typename MatrixKernelClass>
	SymmetryHandler(const MatrixKernelClass *const MatrixKernel, const FReal Epsilon,
									const FReal, const unsigned int)
	{
		// init all 343 item to zero, because effectively only 16 exist
		for (unsigned int t=0; t<343; ++t) {
			K[t] = NULL;
			LowRank[t] = 0;
		}
			
		// set permutation vector and indices
		const FInterpSymmetries<ORDER> Symmetries;
		for (int i=-3; i<=3; ++i)
			for (int j=-3; j<=3; ++j)
				for (int k=-3; k<=3; ++k) {
					const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
					pindices[idx] = 0;
					if (abs(i)>1 || abs(j)>1 || abs(k)>1)
						pindices[idx] = Symmetries.getPermutationArrayAndIndex(i,j,k, pvectors[idx]);
				}

		// precompute 16 M2L operators
		const FReal ReferenceCellWidth = FReal(2.);
		precompute<ORDER>(MatrixKernel, ReferenceCellWidth, Epsilon, K, LowRank);
	}



	/** Destructor */
	~SymmetryHandler()
	{
		for (unsigned int t=0; t<343; ++t) if (K[t]!=NULL) delete [] K[t];
	}


	/*! return the t-th approximated far-field interactions*/
	const FReal *const getK(const unsigned int, const unsigned int t) const
	{	return K[t]; }

	/*! return the t-th approximated far-field interactions*/
	const int getLowRank(const unsigned int, const unsigned int t) const
	{	return LowRank[t]; }

};






/*! Specialization for non-homogeneous kernel functions */
template <int ORDER>
class SymmetryHandler<ORDER, NON_HOMOGENEOUS>
{
  static const unsigned int nnodes = ORDER*ORDER*ORDER;

	// Height of octree; needed only in the case of non-homogeneous kernel functions
	const unsigned int TreeHeight;

	// M2L operators for all levels in the octree
	FReal***    K;
	int** LowRank;

public:
	
	// permutation vectors and permutated indices
	unsigned int pvectors[343][nnodes];
	unsigned int pindices[343];


	/** Constructor: with 16 small SVDs */
	template <typename MatrixKernelClass>
	SymmetryHandler(const MatrixKernelClass *const MatrixKernel, const double Epsilon,
									const FReal RootCellWidth, const unsigned int inTreeHeight)
		: TreeHeight(inTreeHeight)
	{
		// init all 343 item to zero, because effectively only 16 exist
		K       = new FReal** [TreeHeight];
		LowRank = new int*    [TreeHeight];
		K[0]       = NULL; K[1]       = NULL;
		LowRank[0] = NULL; LowRank[1] = NULL;
		for (unsigned int l=2; l<TreeHeight; ++l) {
			K[l]       = new FReal* [343];
			LowRank[l] = new int    [343];
			for (unsigned int t=0; t<343; ++t) {
				K[l][t]       = NULL;
				LowRank[l][t] = 0;
			}
		}
		

		// set permutation vector and indices
		const FInterpSymmetries<ORDER> Symmetries;
		for (int i=-3; i<=3; ++i)
			for (int j=-3; j<=3; ++j)
				for (int k=-3; k<=3; ++k) {
					const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
					pindices[idx] = 0;
					if (abs(i)>1 || abs(j)>1 || abs(k)>1)
						pindices[idx] = Symmetries.getPermutationArrayAndIndex(i,j,k, pvectors[idx]);
				}

		// precompute 16 M2L operators at all levels having far-field interactions
		FReal CellWidth = RootCellWidth / FReal(2.); // at level 1
		CellWidth /= FReal(2.);                      // at level 2
		for (unsigned int l=2; l<TreeHeight; ++l) {
			precompute<ORDER>(MatrixKernel, CellWidth, Epsilon, K[l], LowRank[l]);
			CellWidth /= FReal(2.);                    // at level l+1 
		}
	}



	/** Destructor */
	~SymmetryHandler()
	{
		for (unsigned int l=0; l<TreeHeight; ++l) {
			if (K[l]!=NULL) {
				for (unsigned int t=0; t<343; ++t) if (K[l][t]!=NULL) delete [] K[l][t];
				delete [] K[l];
			}
			if (LowRank[l]!=NULL)	delete [] LowRank[l];
		}
		delete [] K;
		delete [] LowRank;
	}

	/*! return the t-th approximated far-field interactions*/
	const FReal *const getK(const unsigned int l, const unsigned int t) const
	{	return K[l][t]; }

	/*! return the t-th approximated far-field interactions*/
	const int getLowRank(const unsigned int l, const unsigned int t) const
	{	return LowRank[l][t]; }

};








#include <fstream>
#include <sstream>


/**
 * Computes, compresses and stores the 16 M2L kernels in a binary file.
 */
template <int ORDER, typename MatrixKernelClass>
static void ComputeAndCompressAndStoreInBinaryFile(const MatrixKernelClass *const MatrixKernel, const FReal Epsilon)
{
	static const unsigned int nnodes = ORDER*ORDER*ORDER;

	// compute and compress ////////////
	FReal* K[343];
	int LowRank[343];
	for (unsigned int idx=0; idx<343; ++idx) { K[idx] = NULL; LowRank[idx] = 0;	}
	precompute<ORDER>(MatrixKernel, FReal(2.), Epsilon, K, LowRank);

	// write to binary file ////////////
	FTic time; time.tic();
	// start computing process
	const char precision = (typeid(FReal)==typeid(double) ? 'd' : 'f');
	std::stringstream sstream;
	sstream << "sym2l_" << precision << "_o" << ORDER << "_e" << Epsilon << ".bin";
	const std::string filename(sstream.str());
	std::ofstream stream(filename.c_str(),
											 std::ios::out | std::ios::binary | std::ios::trunc);
	if (stream.good()) {
		stream.seekp(0);
		for (unsigned int idx=0; idx<343; ++idx)
			if (K[idx]!=NULL) {
				// 1) write index
				stream.write(reinterpret_cast<char*>(&idx), sizeof(int));
				// 2) write low rank (int)
				int rank = LowRank[idx];
				stream.write(reinterpret_cast<char*>(&rank), sizeof(int));
				// 3) write U and V (both: rank*nnodes * FReal)
				FReal *const U = K[idx];
				FReal *const V = K[idx] + rank*nnodes;
				stream.write(reinterpret_cast<char*>(U), sizeof(FReal)*rank*nnodes);
				stream.write(reinterpret_cast<char*>(V), sizeof(FReal)*rank*nnodes);
			}
	} else throw std::runtime_error("File could not be opened to write");
	stream.close();
	// write info
	//	std::cout << "Compressed M2L operators stored in binary file " << filename
	//					<< " in " << time.tacAndElapsed() << "sec."	<< std::endl;

	// free memory /////////////////////
	for (unsigned int t=0; t<343; ++t) if (K[t]!=NULL) delete [] K[t];
}


/**
 * Reads the 16 compressed M2L kernels from the binary files and writes them
 * in K and the respective low-rank in LowRank.
 */
template <int ORDER>
void ReadFromBinaryFile(const FReal Epsilon, FReal* K[343], int LowRank[343])
{
	// compile time constants
	const unsigned int nnodes = ORDER*ORDER*ORDER;
	
	// find filename
	const char precision = (typeid(FReal)==typeid(double) ? 'd' : 'f');
	std::stringstream sstream;
	sstream << "sym2l_" << precision << "_o" << ORDER << "_e" << Epsilon << ".bin";
	const std::string filename(sstream.str());

	// read binary file
	std::ifstream istream(filename.c_str(),
												std::ios::in | std::ios::binary | std::ios::ate);
	const std::ifstream::pos_type size = istream.tellg();
	if (size<=0) throw std::runtime_error("The requested binary file does not yet exist. Exit.");
	
	if (istream.good()) {
		istream.seekg(0);
		// 1) read index (int)
		int _idx;
		istream.read(reinterpret_cast<char*>(&_idx), sizeof(int));
		// loop to find 16 compressed m2l operators
		for (int idx=0; idx<343; ++idx) {
			K[idx] = NULL;
			LowRank[idx] = 0;
			// if it exists
			if (idx == _idx) {
				// 2) read low rank (int)
				int rank;
				istream.read(reinterpret_cast<char*>(&rank), sizeof(int));
				LowRank[idx] = rank;
				// 3) read U and V (both: rank*nnodes * FReal)
				K[idx] = new FReal [2*rank*nnodes];
				FReal *const U = K[idx];
				FReal *const V = K[idx] + rank*nnodes;
				istream.read(reinterpret_cast<char*>(U), sizeof(FReal)*rank*nnodes);
				istream.read(reinterpret_cast<char*>(V), sizeof(FReal)*rank*nnodes);

				// 1) read next index
				istream.read(reinterpret_cast<char*>(&_idx), sizeof(int));
			}
		}
	}	else throw std::runtime_error("File could not be opened to read");
	istream.close();
}





#endif
