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


#define FACASVD


/**
 * Fully pivoted adaptive cross approximation (ACA)
 */
template <int ORDER>
static void fACA(const FReal *const R, const double eps, FReal* &U, FReal* &V, unsigned int &k)
{
	static const int nnodes = ORDER*ORDER*ORDER;

	// copy K to R
	FReal *const K = new FReal[nnodes*nnodes];
	for (int i=0; i<nnodes*nnodes; ++i) K[i] = R[i];

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

	//std::cout << "fACA rank = " << k << " and error = " << FMath::Sqrt(norm2R) << std::endl;

	delete [] K;

}


template <int ORDER, typename MatrixKernelClass>
static void precomputeSVD(const MatrixKernelClass *const MatrixKernel, const double Epsilon, FReal* K[343], int LowRank[343])
{
	static const unsigned int nnodes = ORDER*ORDER*ORDER;

	// set all to zero
	for (unsigned int t=0; t<343; ++t) { K[t] = NULL;	LowRank[t] = 0;	}

	// interpolation points of source (Y) and target (X) cell
	FPoint X[nnodes], Y[nnodes];
	// set roots of target cell (X)
	FChebTensor<ORDER>::setRoots(FPoint(0.,0.,0.), FReal(2.), X);
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
					//std::cout << " - recompression with SVD leads to rank = " << rank << std::endl;
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










/**
 * Handler to deal with all symmetries: Stores permutation indices and vectors
 * to reduce 343 different interactions to 16 only.
 */
template <int ORDER>
class SymmetryHandler
{
  static const unsigned int nnodes = ORDER*ORDER*ORDER;

public:
	// M2L operators
	FReal*    K[343];
	int LowRank[343];
	
	// permutation vectors and permutated indices
	unsigned int pvectors[343][nnodes];
	unsigned int pindices[343];


	/** Constructor: with 16 small SVDs */
	template <typename MatrixKernelClass>
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
				for (int k=-3; k<=3; ++k) {
					const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
					pindices[idx] = 0;
					if (abs(i)>1 || abs(j)>1 || abs(k)>1)
						pindices[idx] = Symmetries.getPermutationArrayAndIndex(i,j,k, pvectors[idx]);
				}

		// precompute 16 M2L operators
		precomputeSVD<ORDER>(MatrixKernel, Epsilon, K, LowRank);
	}



	/** Destructor */
	~SymmetryHandler()
	{
		for (unsigned int t=0; t<343; ++t) if (K[t]!=NULL) delete [] K[t];
	}

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
	precomputeSVD<ORDER>(MatrixKernel, Epsilon, K, LowRank);

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
