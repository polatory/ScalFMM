// [--License--]

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FMath.hpp"
#include "../Src/Chebyshev/FChebMatrixKernel.hpp"
#include "../Src/Chebyshev/FChebM2LHandler.hpp"




/**
* In this file we show how to use octree
*/
int main(int argc, char* argv[])
{ 
	// typedefs   
	typedef FChebMatrixKernelR MatrixKernelClass;
	

	// instantiations
	FTic time;
	MatrixKernelClass MatrixKernel;

	/*
	// constants
  const FReal epsilon = FReal(atof(argv[1]));
	const unsigned int order = 4;

	// write precomputed compressed M2l operators to binary file 
	std::cout << "\nCompute compressed M2L operators of ACC("
						<< order << ", " << epsilon << ") and write them to a binary file ..."
						<< std::endl;
	time.tic();
	typedef FChebM2LHandler<order,MatrixKernel> M2LHandlerClass;
	M2LHandlerClass::ComputeAndCompressAndStoreInBinaryFile(epsilon);
	time.tac();
	std::cout << " in " << time.elapsed() << "sec." << std::endl;

	// read precomputed compressed M2l operators from binary file 
	std::cout << "\nRead compressed M2L operators of ACC(" << order << ", " << epsilon << ") from a binary file ..." << std::endl;
	time.tic();
	M2LHandlerClass M2L(epsilon);
	M2L.ReadFromBinaryFileAndSet();
	time.tac();
	std::cout << " in " << time.elapsed() << "sec." << std::endl;
	*/

	// order 2
	time.tic();
	FChebM2LHandler<2,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-1));
	FChebM2LHandler<2,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-2));
	// order 3
	FChebM2LHandler<3,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-2));
	FChebM2LHandler<3,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-3));
	// order 4
	FChebM2LHandler<4,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-3));
	FChebM2LHandler<4,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-4));
	// order 5
	FChebM2LHandler<5,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-4));
	FChebM2LHandler<5,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-5));
	// order 6
	FChebM2LHandler<6,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-5));
	FChebM2LHandler<6,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-6));
	// order 7
	FChebM2LHandler<7,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-6));
	FChebM2LHandler<7,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-7));
	// order 8
	FChebM2LHandler<8,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-7));
	FChebM2LHandler<8,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-8));
	// order 9
	FChebM2LHandler<9,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-8));
	FChebM2LHandler<9,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-9));
	// order 10
	FChebM2LHandler<10,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-9));
	FChebM2LHandler<10,MatrixKernelClass>::ComputeAndCompressAndStoreInBinaryFile(FReal(1e-10));


	return 0;
}


// [--END--]
