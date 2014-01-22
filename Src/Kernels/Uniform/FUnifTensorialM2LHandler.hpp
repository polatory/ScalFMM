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
#ifndef FUNIFTENSORIALM2LHANDLER_HPP
#define FUNIFTENSORIALM2LHANDLER_HPP

#include <numeric>
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>
#include <typeinfo>

#include "../../Utils/FBlas.hpp"
#include "../../Utils/FTic.hpp"
#include "../../Utils/FDft.hpp"

#include "../../Utils/FComplexe.hpp"


#include "./FUnifTensor.hpp"



/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FUnifTensorialM2LHandler
 * Please read the license
 *
 * This class precomputes and efficiently stores the M2L operators
 * \f$[K_1,\dots,K_{316}]\f$ for all (\f$7^3-3^3 = 316\f$ possible interacting
 * cells in the far-field) interactions for the Lagrange interpolation
 * approach. The resulting Lagrange operators have a Circulant Toeplitz 
 * structure and can be applied efficiently in Fourier Space. Hence, the
 * originally \f$K_t\f$ of size \f$\ell^3\times\ell^3\f$ times \f$316\f$ for
 * all interactions is reduced to \f$316\f$ \f$C_t\f$, each of size \f$2\ell-1\f$.
 *
 * @tparam ORDER interpolation order \f$\ell\f$
 */
template <int ORDER, class MatrixKernelClass>
class FUnifTensorialM2LHandler : FNoCopyable
{
	enum {order = ORDER,
				nnodes = TensorTraits<ORDER>::nnodes,
				ninteractions = 316, // 7^3 - 3^3 (max num cells in far-field)
        rc = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1),
        dim = MatrixKernelClass::DIM, 
        nRhs = MatrixKernelClass::NRHS,
        nLhs = MatrixKernelClass::NLHS};

	const MatrixKernelClass MatrixKernel;

	FComplexe *FC[dim];

  // for real valued kernel only n/2+1 complex values are stored 
  // after performing the DFT (the rest is deduced by conjugation)
  unsigned int opt_rc; 

  typedef FUnifTensor<ORDER> TensorType;
  unsigned int node_diff[nnodes*nnodes];

  //  FDft Dft; // Direct Discrete Fourier Transformator
  FFft<1> Dft; // Fast Discrete Fourier Transformator

	static const std::string getFileName()
	{
		const char precision_type = (typeid(FReal)==typeid(double) ? 'd' : 'f');
		std::stringstream stream;
		stream << "m2l_k"<< MatrixKernelClass::Identifier << "_" << precision_type
					 << "_o" << order << ".bin";
		return stream.str();
	}

	
public:
	FUnifTensorialM2LHandler()
		: MatrixKernel(), // PB: used only to getPosition and getScalFac in applyM2L
      opt_rc(rc/2+1),
      Dft(rc) // initialize Discrete Fourier Transformator
	{    
    // initialize array
    for (unsigned int d=0; d<dim; ++d)
      FC[d]=NULL;

    // initialize root node ids
    TensorType::setNodeIdsDiff(node_diff);
  }

	~FUnifTensorialM2LHandler()
	{
    for (unsigned int d=0; d<dim; ++d)
      if (FC[d] != NULL) delete [] FC[d];
	}

	/**
	 * Computes and sets the matrix \f$C_t\f$
	 */
	void ComputeAndSet()
	{
		// measure time
		FTic time; time.tic();

    for (unsigned int d=0; d<dim; ++d){

      // check if aready set
      if (FC[d]) throw std::runtime_error("M2L operator already set");
      // Compute matrix of interactions
      Compute(FC[d], d);

    }

    // Compute memory usage
    unsigned long sizeM2L = 343*dim*opt_rc*sizeof(FComplexe);


		// write info
		std::cout << "Compute and Set M2L operators ("<< long(sizeM2L/**1e-6*/) <<" Bytes) in "
							<< time.tacAndElapsed() << "sec."	<< std::endl;
	}

	/**
	 * Computes, writes to binary file, reads it and sets the matrices \f$Y, C_t, B\f$
	 */
	void ComputeAndStoreInBinaryFileAndReadFromFileAndSet()
	{
		FUnifTensorialM2LHandler<ORDER,MatrixKernelClass>::ComputeAndStoreInBinaryFile();
		this->ReadFromBinaryFileAndSet();
	}

	/**
	 * Computes all \f$K_t\f$.
	 *
	 * @param[out] C matrix of size \f$r\times 316 r\f$ storing \f$[C_1,\dots,C_{316}]\f$
	 */
	static void Compute(FComplexe* &FC, const unsigned int d);

	/**
	 * Computes and stores the matrix \f$C_t\f$ in a binary
	 * file
	 */
	static void ComputeAndStoreInBinaryFile();

	/**
	 * Reads the matrices \f$Y, C_t, B\f$ from the respective binary file
	 */
	void ReadFromBinaryFileAndSet();
		

  /**
	 * Expands potentials \f$x+=IDFT(X)\f$ of a target cell. This operation can be
	 * seen as part of the L2L operation.
	 *
	 * @param[in] X transformed local expansion of size \f$r\f$
	 * @param[out] x local expansion of size \f$\ell^3\f$
	 */
  void unapplyZeroPaddingAndDFT(const FComplexe *const FX, FReal *const x) const
  { 
    FReal Px[rc];
    FBlas::setzero(rc,Px);
    // Apply forward Discrete Fourier Transform
    Dft.applyIDFT(FX,Px);

    // Unapply Zero Padding
    for (unsigned int j=0; j<nnodes; ++j)
      x[j]+=Px[node_diff[nnodes-j-1]];
  }

  /**
	 * The M2L operation \f$X+=C_tY\f$ is performed in Fourier space by 
   * application of the convolution theorem (i.e. \f$FX+=FC_t:FY\f$), 
   * where \f$FY\f$ is the transformed multipole expansion and 
   * \f$FX\f$ is the transformed local expansion, both padded with zeros and
	 * of size \f$r_c=(2\times ORDER-1)^3\f$ or \f$r_c^{opt}=\frac{r_c}{2}+1\f$ 
   * for real valued kernels. The index \f$t\f$ denotes the transfer vector 
	 * of the target cell to the source cell.
	 *
	 * @param[in] transfer transfer vector
	 * @param[in] FY transformed multipole expansion
	 * @param[out] FX transformed local expansion
	 * @param[in] CellWidth needed for the scaling of the compressed M2L operators which are based on a homogeneous matrix kernel computed for the reference cell width \f$w=2\f$, ie in \f$[-1,1]^3\f$.
	 */
  void applyFC(const unsigned int idx, FReal CellWidth,
               const FComplexe *const FY, FComplexe *const FX, 
               const unsigned int idxK) const
  {
		const FReal scale(MatrixKernel.getScaleFactor(CellWidth));
    unsigned int d = MatrixKernel.getPosition(idxK);

    FComplexe tmpFX; 

    // Perform entrywise product manually
    for (unsigned int j=0; j<opt_rc; ++j){
      tmpFX=FC[d][idx*opt_rc + j];
      tmpFX*=FY[j];
      tmpFX*=scale;
      FX[j]+=tmpFX;
    }

//    // Perform entrywise product using BLAS and MKL routines
//    // PB: not necessary faster than the naive version
//    FComplexe tmpFX[rc]; 
//    FBlas::c_setzero(rc,reinterpret_cast<FReal*>(tmpFX));
//    FMkl::c_had(rc,reinterpret_cast<const FReal* const>(FC + idx*rc),
//                reinterpret_cast<const FReal* const>(FY),
//                reinterpret_cast<FReal* const>(tmpFX));
//    // Scale
//    FBlas::c_axpy(rc,&scale,reinterpret_cast<FReal* const>(tmpFX),
//                  reinterpret_cast<FReal* const>(FX));

  }


  /**
	 * Transform densities \f$Y= DFT(y)\f$ of a source cell. This operation
	 * can be seen as part of the M2M operation.
	 *
	 * @param[in] y multipole expansion of size \f$\ell^3\f$
	 * @param[out] Y transformed multipole expansion of size \f$r\f$
	 */
  void applyZeroPaddingAndDFT(FReal *const y, FComplexe *const FY) const
  {
    FReal Py[rc];
    FBlas::setzero(rc,Py);

    // Apply Zero Padding
    for (unsigned int i=0; i<nnodes; ++i)
      Py[node_diff[i*nnodes]]=y[i];

    // Apply forward Discrete Fourier Transform
    Dft.applyDFT(Py,FY);

  }


};






//////////////////////////////////////////////////////////////////////
// definition ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////






template <int ORDER, class MatrixKernelClass>
void
FUnifTensorialM2LHandler<ORDER, MatrixKernelClass>::Compute(FComplexe* &FC,
                                                            const unsigned int d)
{
	// allocate memory and store compressed M2L operators
  	if (FC) throw std::runtime_error("M2L operators are already set");

	// interpolation points of source (Y) and target (X) cell
	FPoint X[nnodes], Y[nnodes];
	// set roots of target cell (X)
	FUnifTensor<order>::setRoots(FPoint(0.,0.,0.), FReal(2.), X);
	// init matrix kernel
	const MatrixKernelClass MatrixKernel(d);

	// allocate memory and compute 316 m2l operators
	FReal *_C;
	FComplexe *_FC;

  // reduce storage from nnodes^2=order^6 to (2order-1)^3
  const unsigned int rc = (2*order-1)*(2*order-1)*(2*order-1);
	_C = new FReal [rc];
	_FC = new FComplexe [rc * ninteractions];

  // initialize root node ids pairs
  unsigned int node_ids_pairs[rc][2];
  TensorType::setNodeIdsPairs(node_ids_pairs);

  // init Discrete Fourier Transformator
//	FDft Dft(rc);
	FFft<1> Dft(rc);

  // get first column of K via permutation
  unsigned int perm[rc];
  TensorType::setStoragePermutation(perm);

	unsigned int counter = 0;
	for (int i=-3; i<=3; ++i) {
		for (int j=-3; j<=3; ++j) {
			for (int k=-3; k<=3; ++k) {
				if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
					// set roots of source cell (Y)
					const FPoint cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
					FUnifTensor<order>::setRoots(cy, FReal(2.), Y);
					// evaluate m2l operator
          unsigned int ido=0;
          for(unsigned int l=0; l<2*order-1; ++l)
            for(unsigned int m=0; m<2*order-1; ++m)
              for(unsigned int n=0; n<2*order-1; ++n){   
          
                // store value at current position in C
                // use permutation if DFT is used because 
                // the storage of the first column is required
                // i.e. C[0] C[rc-1] C[rc-2] .. C[1] < WRONG!
                // i.e. C[rc-1] C[0] C[1] .. C[rc-2] < RIGHT!
                //                _C[counter*rc + ido]
                _C[perm[ido]]
                  = MatrixKernel.evaluate(X[node_ids_pairs[ido][0]], 
                                          Y[node_ids_pairs[ido][1]]);
                ido++;
              }

          // Apply Discrete Fourier Transformation
          Dft.applyDFT(_C,_FC+counter*rc);

					// increment interaction counter
					counter++;
				}
			}
		}
	}
	if (counter != ninteractions)
		throw std::runtime_error("Number of interactions must correspond to 316");

  // Free _C
	delete [] _C;

	// store FC
	counter = 0;
  // reduce storage if real valued kernel
  const unsigned int opt_rc = rc/2+1;
  // allocate M2L
	FC = new FComplexe[343 * opt_rc];

	for (int i=-3; i<=3; ++i)
		for (int j=-3; j<=3; ++j)
			for (int k=-3; k<=3; ++k) {
				const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
				if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
					FBlas::c_copy(opt_rc, reinterpret_cast<FReal*>(_FC + counter*rc), 
                        reinterpret_cast<FReal*>(FC + idx*opt_rc));
//          for (unsigned int n=0; n<rc; ++n){
//            FC[idx*rc+n]=_FC[counter*rc+n];
//          }
					counter++;
				} else{ 
          FBlas::c_setzero(opt_rc, reinterpret_cast<FReal*>(FC + idx*opt_rc));
//          for (unsigned int n=0; n<rc; ++n){
//            FC[idx*rc+n]=FComplexe(0.0,0.0);
//          }
        }
      }

	if (counter != ninteractions)
		throw std::runtime_error("Number of interactions must correspond to 316");
  delete [] _FC;  	
}






//template <int ORDER, class MatrixKernelClass>
//void
//FUnifTensorialM2LHandler<ORDER, MatrixKernelClass>::ComputeAndStoreInBinaryFile()
//{
//	// measure time
//	FTic time; time.tic();
//	// start computing process
//	FReal *C;
//	C = NULL;
//	const unsigned int rc = Compute(C);
//	// store into binary file
//	const std::string filename(getFileName());
//	std::ofstream stream(filename.c_str(),
//											 std::ios::out | std::ios::binary | std::ios::trunc);
//	if (stream.good()) {
//		stream.seekp(0);
//		// 1) write number of interpolation points (int)
//		int _nnodes = nnodes;
//		stream.write(reinterpret_cast<char*>(&_nnodes), sizeof(int));
//		// 2) write 343 C (343 * rc * FReal)
//		stream.write(reinterpret_cast<char*>(C), sizeof(FReal)*rc*343);
//	} 	else throw std::runtime_error("File could not be opened to write");
//	stream.close();
//	// free memory
//	if (C != NULL) delete [] C;
//	// write info
//	std::cout << "M2L operators (r2="<< rc 
//            << ",nnodes2="<< nnodes*nnodes 
//            <<") stored in binary file "	<< filename
//						<< " in " << time.tacAndElapsed() << "sec."	<< std::endl;
//}
//
//
//template <int ORDER, class MatrixKernelClass>
//void
//FUnifTensorialM2LHandler<ORDER, MatrixKernelClass>::ReadFromBinaryFileAndSet()
//{
//	// measure time
//	FTic time; time.tic();
//	// start reading process
//	if (C) throw std::runtime_error("M2L operator already set");
//	const std::string filename(getFileName());
//	std::ifstream stream(filename.c_str(),
//											 std::ios::in | std::ios::binary | std::ios::ate);
//	const std::ifstream::pos_type size = stream.tellg();
//	if (size<=0) {
//		std::cout << "Info: The requested binary file " << filename
//							<< " does not yet exist. Compute it now ... " << std::endl;
//		this->ComputeAndStoreInBinaryFileAndReadFromFileAndSet();
//		return;
//	} 
//	if (stream.good()) {
//		stream.seekg(0);
//		// 1) read number of interpolation points (int)
//		int npts;
//		stream.read(reinterpret_cast<char*>(&npts), sizeof(int));
//		if (npts!=nnodes) throw std::runtime_error("nnodes and npts do not correspond");
//		// 2) write 343 C (343 * rank*rank * FReal)
//		C = new FReal [343 * rc];
//		stream.read(reinterpret_cast<char*>(C), sizeof(FReal)*rc*343);
//	}	else throw std::runtime_error("File could not be opened to read");
//	stream.close();
//	// write info
//	std::cout << "M2L operators (rc=" << rc 
//            <<",nnodes2="<< nnodes*nnodes 
//            << ") read from binary file "
//						<< filename << " in " << time.tacAndElapsed() << "sec."	<< std::endl;
//}

/*
unsigned int ReadRankFromBinaryFile(const std::string& filename)
{
	// start reading process
	std::ifstream stream(filename.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
	const std::ifstream::pos_type size = stream.tellg();
	if (size<=0) throw std::runtime_error("The requested binary file does not exist.");
	unsigned int rank = -1;
	if (stream.good()) {
		stream.seekg(0);
		// 1) read number of interpolation points (int)
		int npts;
		stream.read(reinterpret_cast<char*>(&npts), sizeof(int));
		// 2) read low rank (int)
		stream.read(reinterpret_cast<char*>(&rank), sizeof(int));
		return rank;
	}	else throw std::runtime_error("File could not be opened to read");
	stream.close();
	return rank;
}
*/






#endif // FUNIFTENSORIALM2LHANDLER_HPP
 
// [--END--]
