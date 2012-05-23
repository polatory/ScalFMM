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
