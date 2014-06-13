/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */


#include <iostream>
#include <fstream>
//#include <sstream>
#include <string>
#include "Utils/FGlobal.hpp"
#include "Utils/FPoint.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FGenerateDistribution.hpp"

//
/// \file  changeFormat.cpp
//!
//! \brief changeFormat: Driver to transform a FMA format and to build a visualization file
//!
//!
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary)
//!
//!     \param   -fout name: generic name for files (without extension) and save data
//!                  with following format in name.fma or name.bfma if -bin is set"
//!      \param   -bin save output in binary mode (name file  name.bfma
//!      \param  -visufmt format for the visu file (vtk, vtp, cvs or cosmo). vtp is the default
//!
//!
//! \b examples
//!
//!  Transform an ascii file in a binary file
//!
//!   changeFormat -filename unitCubeXYZQ100.fma  -outfilename unitCubeXYZQ100 -bin


//
//
void genDistusage() {
	std::cout << "Driver to change the format of the input file"
			<< std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "   -help       to see the parameters    " << std::endl
			<<     "   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary) " <<std::endl
			<<     "     Output " << std::endl
			<<     "             -fout name: generic name for files (without extension) and save data" <<std::endl
			<<     "                     with following format in name.fma or name.bfma if -bin is set" <<std::endl
			<<     "      -bin save output in binary mode (name file  name.bfma)" <<std::endl
			<<     "      -visufmt  vtk, vtp, cosmo or cvs format for visualization " <<std::endl;
}

int main(int argc, char ** argv){
	//
	if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")|| (argc < 3 )){
		genDistusage() ;
		exit(-1);
	}
	const std::string filename(FParameters::getStr(argc,argv,"-fin",   "data.fma"));

	FFmaGenericLoader loader(filename);
	//
	// Allocation
	//
	long int NbPoints     = loader.getNumberOfParticles();
	const unsigned int nbData   = loader.getNbRecordPerline() ;
	const unsigned int arraySize =nbData*NbPoints;
    FReal	 * particles ;
	particles = new FReal[arraySize] ;
	memset(particles,0,arraySize*sizeof(FReal));
	//
	// Read Data
	int j = 0, k=0 ;
	for(int idxPart = 0 ; idxPart < NbPoints ;++idxPart, j+=nbData){
//		//
		loader.fillParticle(&particles[j],nbData);
//		std::cout << "idxPart "<< idxPart << "  ";
//		for (int jj= 0 ; jj<nbData ; ++jj, ++k){
//			std::cout << particles[k] << "    ";
//		}
//		std::cout << std::endl;
	}
	//
	/////////////////////////////////////////////////////////////////////////
	//                                           Save data
	/////////////////////////////////////////////////////////////////////////
	//
	//  Generate file for ScalFMM Loader
	//
	bool binaryMode = false;
    std::string outfilename(FParameters::getStr(argc,argv,"-fout",   "output"));
	if(  FParameters::existParameter(argc, argv, "-bin")){
		binaryMode = true;
		outfilename += ".bfma";
	}
	else {
		outfilename += ".fma";
	}
	FFmaGenericWriter writer(outfilename,binaryMode) ;
	writer.writeHeader( loader.getCenterOfBox(), loader.getBoxWidth() , NbPoints, sizeof(FReal), nbData) ;
	writer.writeArrayOfReal(particles, nbData, NbPoints);
	//
	//
	//
	if(FParameters::existParameter(argc, argv, "-visufmt")){
		std::string visufile(""), fmt(FParameters::getStr(argc,argv,"-visufmt",   "vtp"));
		if( fmt == "vtp" ){
			visufile = outfilename + ".vtp" ;
		}
		else	if( fmt == "vtk" ){
			visufile = outfilename + ".vtk" ;
		}
		else if( fmt == "cosmo" ){
			visufile = outfilename + ".cosmo" ;
		}
		else {
			visufile =   outfilename + ".csv" ;
		}
		std::ofstream file( visufile, std::ofstream::out);
		if(!file) {
			std::cout << "Cannot open file."<< std::endl;
			exit(-1)	;
		}	//
		//
		// Export data in cvs format
		//
		if( fmt == "vtp" ){
			std::cout << "Writes in XML VTP format  (visualization) in file "<< visufile <<std::endl ;
			exportVTKxml( file, NbPoints,  particles)   ;
		}
		else		if( fmt == "vtk" ){
			std::cout << "Writes in VTK format  (visualization) in file "<< visufile <<std::endl ;
			exportVTK( file, NbPoints, particles)  ;
		}
		else if( fmt == "cosmo" ){
			std::cout << "Writes in COSMO format  (visualization) in file "<< visufile <<std::endl ;
			exportCOSMOS( file, NbPoints,  particles)   ;
		}
		else {
			std::cout << "Writes in CVS format  (visualization) in file "<<visufile<<std::endl ;
			exportCVS( file, NbPoints,  particles)  ;
		}
	}
	//
	delete particles ;

	//
	return 1;
}
