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
#include "Files/FDlpolyLoader.hpp"
#include "Utils/FParameters.hpp"
#include "Files/FGenerateDistribution.hpp"
#include "Utils/FParameterNames.hpp"
#include "Files/FExportWriter.hpp"


//
/// \file  stamp2Fma.cpp
//!
//! \brief changeFormat: Driver to transform a stamp format to FMA format and to build a visualization file
//!
//!  Driver to transform a FMA format and/or to build a visualization file<br>
//! For a description of the FMA format see FFmaGenericLoader<br>
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary)
//!       \param  -fdlpoly name  file coming from a DLpoly simulation
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
//!   changeFormat -fin unitCubeXYZQ100.fma  -fout unitCubeXYZQ100 -bin


//
//
void genDistusage() {
	std::cout << "Driver to change the format of the input file"
			<< std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "   -help       to see the parameters    " << std::endl
			<<     " Input:  only one option is allowed"  << std::endl
			<<     "   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary) " <<std::endl
			<<     "   -fdlpoly name:  file name  to convert with extension (.bin if binary file) " <<std::endl;
}

int main(int argc, char ** argv){
	//
    FHelpDescribeAndExit(argc, argv,
                         "Driver to change the format of the input file. "
                         "fdlpoly is not supported for now.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::OutputVisuFile,);

	if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")|| (argc < 3 )){
		genDistusage() ;
		exit(-1);
	}
	FSize NbPoints;
    typedef double FReal;
	FReal	 * particles = nullptr ;
    FFmaGenericLoader<FReal> * loader;
    FSize nbData;
	bool stampFile = false ;
	if (FParameters::existParameter(argc, argv, "-fstamp")) {
		stampFile  = true ;
	}
	if (FParameters::existParameter(argc, argv, "-fin")) {
		const std::string filename(FParameters::getStr(argc,argv,"-fin",   "data.fma"));
		if(stampFile) {
            loader = new FFmaGenericLoader<FReal>(filename,false);
		}
		else {
            loader = new FFmaGenericLoader<FReal>(filename) ;
		}
		//
		// Allocation
		//
		NbPoints                              = loader->getNumberOfParticles();
		nbData                                = loader->getNbRecordPerline() ;
        const FSize arraySize =nbData*NbPoints;
		//
		particles = new FReal[arraySize] ;
		std::memset(particles,0,arraySize*sizeof(FReal));
		//
		// Read Data
        FSize j = 0 ;
		for(FSize idxPart = 0 ; idxPart < NbPoints ;++idxPart, j+=nbData){
			//		//
			loader->fillParticle(&particles[j],nbData);
			//		std::cout << "idxPart "<< idxPart << "  ";
			//		for (int jj= 0 ; jj<nbData ; ++jj, ++k){
			//			std::cout << particles[k] << "    ";
			//		}
			//		std::cout << std::endl;
		}
		if(stampFile) {
			j=0 ;
			FReal D= loader->getBoxWidth();
			for(FSize idxPart = 0 ; idxPart < NbPoints ;++idxPart, j+=nbData){
				particles[j]	         *= D ;
				particles[j+1]	 *= D ;
				particles[j+2]	 *= D ;
			}
		}
		//		if(FParameters::existParameter(argc, argv, "-fdlpoly")){
		//			FDlpolyLoader  *loader = nullptr ;
		//			//			if(FParameters::existParameter(argc, argv, "-bin")){
		// //				loader  = new FDlpolyBinLoader(filenameEwaldIn.c_str());
		// //			}
		// //			else {
		// //				loader  = new FDlpolyAsciiLoader(filenameEwaldIn.c_str());
		// //			}
		//			//				NbPoints = loader->getNumberOfParticles()  ;
		//			//				particles = new FReal[arraySize] ;
		//			//				std::memset(particles,0,arraySize*sizeof(FReal));
		//			//				for(FSize idxPart = 0 ; idxPart < NbPoints ; ++idxPart){
		//				//
		//			//					int index ;
        //			//				FPoint<FReal> P ; FReal t[3];
		// ///	/				loader->fillParticle(&P, t, &physicalValue,&index);
		//			//				particles[(index-1)*]
		//				//
		// //				totalCharge += physicalValue ;
		//				}
	}
	else {
		genDistusage() ;
		return 0;
	}


	//
	/////////////////////////////////////////////////////////////////////////
	//                                           Save data
	/////////////////////////////////////////////////////////////////////////
	//
	//  Generate file for ScalFMM FMAGenericLoader
	//
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
        std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
          FFmaGenericWriter<FReal> writer(name) ;
  		writer.writeHeader( loader->getCenterOfBox(), loader->getBoxWidth() , NbPoints, sizeof(FReal), nbData) ;
  		writer.writeArrayOfReal(particles, nbData, NbPoints);
    }

	//
	//   Generate file for visualization purpose
	//
	   if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputVisuFile.options)){
	        std::string outfilename(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.vtp"));
	        driverExportData(outfilename, particles , NbPoints);
	   }
	//
    delete[] particles ;

	//
	return 1;
}
