#ifndef FCHEBSYMM2LHANDLER_HPP
#define FCHEBSYMM2LHANDLER_HPP

#include <climits>

#include "../../Utils/FBlas.hpp"

#include "./FChebTensor.hpp"
#include "./FChebSymmetries.hpp"
#include "./FChebM2LHandler.hpp"

/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */


/*!
  If \a FACASVD is defined ACA+SVD will be used instead of only SVD.
*/
#define FACASVD


/*!  The fully pivoted adaptive cross approximation (fACA) compresses a
	far-field interaction as \f$K\sim UV^\top\f$. The fACA requires all entries
	to be computed beforehand, then the compression follows in
	\f$\mathcal{O}(2\ell^3k)\f$ operations based on the required accuracy
	\f$\varepsilon\f$.

	@param[in] K far-field to be approximated
	@param[in] eps prescribed accuracy
	@param[out] U matrix containing \a k column vectors
	@param[out] V matrix containing \a k row vectors
	@param[out] k final low-rank depends on prescribed accuracy \a eps
*/
template <int ORDER>
static void fACA(FReal *const K, const double eps, FReal* &U, FReal* &V, unsigned int &k)
{
	static const int nnodes = ORDER*ORDER*ORDER;

	// control vectors (true if not used, false if used)
	bool r[nnodes], c[nnodes];
	for (int i=0; i<nnodes; ++i) r[i] = c[i] = true;

	// compute Frobenius norm of original Matrix K
	FReal norm2K = 0;
	for (int j=0; j<nnodes; ++j) {
		const FReal *const colK = K + j*nnodes;
		for (int i=0; i<nnodes; ++i) norm2K += colK[i]*colK[i];
	}

	// initialize rank k and UV'
	k = 0;
	const int size = nnodes * nnodes/2;
	U = new FReal[size];
	V = new FReal[size];
	for (int i=0; i<size; ++i) U[i] = V[i] = 0.;
	FReal norm2R;

	////////////////////////////////////////////////
	// start fully pivoted ACA
	do {
		
		// find max(K) and argmax(K)
		FReal maxK = 0.;
		int pi=0, pj=0;
		for (int j=0; j<nnodes; ++j)
			if (c[j]) {
				const FReal *const colK = K + j*nnodes;
				for (int i=0; i<nnodes; ++i)
					if (r[i] && maxK < FMath::Abs(colK[i])) {
						maxK = FMath::Abs(colK[i]);
						pi = i; 
						pj = j;
					}
			}

		// copy pivot cross into U and V
		FReal *const colU = U + k*nnodes;
		FReal *const colV = V + k*nnodes;
		const FReal pivot = K[pj*nnodes + pi];
		for (int i=0; i<nnodes; ++i) {
			if (r[i]) colU[i] = K[pj*nnodes + i];
			if (c[i]) colV[i] = K[i *nnodes + pi] / pivot;
		}
		
		// dont use these cols and rows anymore
		c[pj] = false;
		r[pi] = false;
		
		// subtract k-th outer product from K
		for (int j=0; j<nnodes; ++j)
			if (c[j]) {
				FReal *const colK = K + j*nnodes;
				for (int i=0; i<nnodes; ++i)
					if (r[i])	colK[i] -= colU[i] * colV[j];
			}

		// compute Frobenius norm of updated K
		norm2R = 0.;
		for (int j=0; j<nnodes; ++j)
			if (c[j]) {
				const FReal *const colK = K + j*nnodes;
				for (int i=0; i<nnodes; ++i)
					if (r[i])	norm2R += colK[i]*colK[i];
			}

		// increment rank k
		k++;
		
	} while (norm2R > eps*eps * norm2K);
	////////////////////////////////////////////////
}


/*!  Precomputes the 16 far-field interactions (due to symmetries in their
  arrangement all 316 far-field interactions can be represented by
  permutations of the 16 we compute in this function). Depending on whether
  FACASVD is defined or not, either ACA+SVD or only SVD is used to compress
  them. */
template <int ORDER, typename MatrixKernelClass>
static void precompute(const MatrixKernelClass *const MatrixKernel, const FReal CellWidth,
											 const double Epsilon, FReal* K[343], int LowRank[343])
{
	std::cout << "\nComputing 16 far-field interactions for cells of width w = " << CellWidth << std::endl;

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
		
	unsigned int counter = 0;
	for (int i=2; i<=3; ++i) {
		for (int j=0; j<=i; ++j) {
			for (int k=0; k<=j; ++k) {

				// assemble matrix
				const FPoint cy(CellWidth*FReal(i), CellWidth*FReal(j), CellWidth*FReal(k));
				FChebTensor<ORDER>::setRoots(cy, CellWidth, Y);
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


				
				//////////////////////////////////////////////////////////////
#ifdef FACASVD ///////////////////////////////////////////////////////
				// fully pivoted ACA
				FReal *UU, *VV;
				unsigned int rank;
				fACA<ORDER>(U, Epsilon, UU, VV, rank);

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

				std::cout << "(" << i << "," << j << "," << k << ") " << idx <<
					", low rank = " << rank << " (" << aca_rank << ")" << std::endl;

#else
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

				std::cout << "(" << i << "," << j << "," << k << ") " << idx <<
					", low rank = " << rank << std::endl;

#endif ///////////////////////////////////////////////////////////////
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
	std::cout << "num interactions = " << counter << std::endl;
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
	SymmetryHandler(const MatrixKernelClass *const MatrixKernel, const double Epsilon,
									const FReal, const unsigned int)
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
		const FChebSymmetries<ORDER> Symmetries;
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
	std::cout << "Compressed M2L operators stored in binary file " << filename
						<< " in " << time.tacAndElapsed() << "sec."	<< std::endl;

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
