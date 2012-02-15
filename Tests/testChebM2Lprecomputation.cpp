// system includes
#include <iostream>
#include <stdexcept>


#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/F3DPosition.hpp"
#include "../Src/Utils/FMath.hpp"
#include "../Src/Utils/FTic.hpp"


#include "../Src/Chebyshev/FChebTensor.hpp"
#include "../Src/Chebyshev/FChebM2LHandler.hpp"
#include "../Src/Chebyshev/FChebMatrixKernel.hpp"






int main(int argc, char* argv[])
{
	// start timer /////////////////////////////////
	FTic time;
	
	// define set matrix kernel
	typedef FChebMatrixKernelR MatrixKernelClass;
	MatrixKernelClass MatrixKernel;

	// constants
  const FReal epsilon     = FReal(atof(argv[1]));
  const unsigned int order = 5;
	
	// number of interpolation points per cell
	const unsigned int nnodes = TensorTraits<order>::nnodes;

	// interpolation points of source (Y) and target (X) cell
	F3DPosition X[nnodes], Y[nnodes];
	// set roots of target cell (X)
	FChebTensor<order>::setRoots(F3DPosition(0.,0.,0.), FReal(2.), X);
	
	// allocate memory
	FReal *Qu, *C, *Qb;
	Qu = Qb = NULL;
	C = new FReal [nnodes*nnodes * 316];

	////////////////////////////////////////////////
	std::cout << "\nAssembly of 316 times "
						<< nnodes << "x" << nnodes << " M2L operators";
	time.tic();
	// compute 316 m2l operators
	unsigned int counter = 0;
	for (int i=-3; i<=3; ++i) {
		for (int j=-3; j<=3; ++j) {
			for (int k=-3; k<=3; ++k) {
				if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
					// set roots of source cell (Y)
					const F3DPosition cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
					FChebTensor<order>::setRoots(cy, FReal(2.), Y);
					// evaluate m2l operator
					for (unsigned int n=0; n<nnodes; ++n)
						for (unsigned int m=0; m<nnodes; ++m)
							C[counter*nnodes*nnodes + n*nnodes + m] = MatrixKernel.evaluate(X[m], Y[n]);
					// increment interaction counter
					counter++;
				}
			}
		}
	}
	if (counter != 316)
		std::runtime_error("Number of interactions must correspond to 316");
	time.tac();
	std::cout << " took " << time.elapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	////////////////////////////////////////////////
	std::cout << "\nSVD compression ";
	time.tic();
	const unsigned int rank = Compress<order>(epsilon, Qu, C, Qb);
	time.tac();
	std::cout << "to low rank = " << rank << " (eps = " << epsilon
						<< ") took " << time.elapsed() << " sec." << std::endl;
	////////////////////////////////////////////////

	// free memory
	delete [] C;
	if (Qu != NULL) delete Qu;
	if (Qb != NULL) delete Qb;


	return 0;
} 	
