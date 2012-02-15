#include <iostream>
#include <stdlib.h>

#include "../Src/Utils/FBlas.hpp"


FReal FRandom() { return (FReal(rand()) / FReal(RAND_MAX)); }

/**
 * Test functionality of C - interfaced BLAS functions
 */

int main()
{
	const unsigned int m = 4, n = 4; // to be able to test both, transpose and not transpose operations
	FReal* A = new FReal [m * n]; // matrix: column major ordering
	FReal* x = new FReal [n];
	FReal* y = new FReal [m];
	
	const FReal d = FRandom();
	
	for (unsigned int j=0; j<n; ++j) {
		x[j] = FRandom();
		for (unsigned int i=0; i<m; ++i) {
			A[j*m + i] = FRandom();
		}
	}

//	std::cout << "A = " << std::endl;
//	for (unsigned int i=0; i<m; ++i) {
//		for (unsigned int j=0; j<n; ++j) 
//			std::cout << A[j*m + i] << " ";
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;	
//	
//	std::cout << "x = " << std::endl;
//	for (unsigned int i=0; i<m; ++i)
//		std::cout << x[i] << std::endl;
//	std::cout << std::endl;	


	
	FReal* z = new FReal [m];
	

	// y = d Ax ////////////////////////////////////
	// cblas
	FBlas::gemv(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		z[i] = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			z[i] += A[j*m + i] * x[j];
		z[i] *= d;
	}
	// compare
	std::cout << "\ny = d Ax (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;


	// y = d A^Tx ////////////////////////////////////
	// cblas
	FBlas::gemtv(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		z[i] = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			z[i] += A[i*m + j] * x[j];
		z[i] *= d;
	}
	// compare
	std::cout << "\ny = d A^Tx (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;


	// y += d Ax ////////////////////////////////////
	// cblas
	FBlas::gemva(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		FReal _z = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			_z += A[j*m + i] * x[j];
		z[i] += _z * d;
	}
	// compare
	std::cout << "\ny += d Ax (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;


	// y += d A^Tx ////////////////////////////////////
	// cblas
	FBlas::gemtva(m, n, d, A, x, y);
	// direct
	for (unsigned int i=0; i<m; ++i) {
		FReal _z = FReal(0.);
		for (unsigned int j=0; j<n; ++j) 
			_z += A[i*m + j] * x[j];
		z[i] += _z * d;
	}
	// compare
	std::cout << "\ny += d A^Tx (zeros are correct)" << std::endl;
	for (unsigned int i=0; i<m; ++i)
		std::cout << z[i] - y[i] << std::endl;



	delete [] A;
	delete [] x;
	delete [] y;
	delete [] z;

	
	return 0;
}

