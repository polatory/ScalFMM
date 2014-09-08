/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
//
#include "Utils/FGlobal.hpp"
#include "Utils/FPoint.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FCompareResults.hpp"

//
/// \file  compare2files.cpp
//!
//! \brief compare2files: Gives the error between
//!
//!  Driver to transform a FMA format and/or to build a visualization file<br>
//! For a description of the FMA format see FFmaGenericLoader<br>
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fmmfile1 name1:  first file name to compare (with extension .fma (ascii) or bfma (binary)
//!     \param   -fmmfile2 name2: second file name to compare (with extension .fma (ascii) or bfma (binary)
//!     \param   -ewaldfile2 name2 if name2 contains the result done by the ewald method for 1/r kernel.
//!
//! Error code (return of the main)
//!      -1   Number of points is different in the two files.
//!
//! \b examples
//!
//!
//!   compare2files -file1 unitCubeXYZQ100.fma  -file2 unitCubeXYZQ100


//
//
void genDistusage() {
	std::cout << "Driver to change the format of the input file"
			<< std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "   -help       to see the parameters    " << std::endl
			<<     "   -fmmfile1 name1  first file name to compare (with extension .fma (ascii) or bfma (binary) " <<std::endl
			<<     "   -fmmfile2 name2 second file name to compare (with extension .fma (ascii) or bfma (binary)" <<std::endl
			<<     "   -ewaldfile2 name2 if name2 contains the result done by the ewald method for 1/r kernel "<<std::endl
			<<std::endl;
}

int main(int argc, char ** argv){
	//
	if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")|| (argc < 3 )){
		genDistusage() ;
		exit(-1);
	}
	if(FParameters::existParameter(argc, argv, "-ewaldfile2")&&FParameters::existParameter(argc, argv, "-fmmfile2")){
		std::cout << "Either -file2 or -ewald have to be set"<<std::endl;
		exit(-1);
	}

	const std::string filename1(FParameters::getStr(argc,argv,"-fmmfile1",   "data.fma"));
	const std::string filename2(FParameters::getStr(argc,argv,"-fmmfile2",   "data.fma"));

	FFmaGenericLoader loader1(filename1);
	FFmaGenericLoader loader2(filename2);
	//
	// Allocation
	//
	FSize nbParticles                    = loader1.getNumberOfParticles();
	const unsigned int nbData   = loader1.getNbRecordPerline() ;
	if(nbParticles !=  loader2.getNumberOfParticles()){
		std::cerr << "Number of points is different in the two files."<<std::endl ;
		return -1 ;
	}
	if( (nbData !=  8) && (loader2.getNbRecordPerline()) ){
		std::cerr << "Wrong files only " << std::min(loader2.getNbRecordPerline(),nbData)<<  " to read."<<std::endl ;
		return -2 ;
	}
	FmaRWParticle<8,8>* const particles1 = new FmaRWParticle<8,8>[nbParticles];
	FmaRWParticle<8,8>* const particles2 = new FmaRWParticle<8,8>[nbParticles];


	//
	loader1.fillParticle(particles1,nbParticles);
	loader2.fillParticle(particles2,nbParticles);
	if(FParameters::existParameter(argc, argv, "-ewaldfile2") ) {
		FReal volume =1.0 ;
//		double volume               =  boxsize[0] *boxsize[1]*boxsize[2] ;
		 removeFirstMoment( "DLPOLY", nbParticles, particles2,  volume) ;
			FPoint FirstMoment ;
	}
	//
	int error = 0 ;
	error = compareTwoArrays("TAG", nbParticles, particles1, particles2);

	//
	delete particles1 ;
	delete particles2 ;

	//
	return error;
}