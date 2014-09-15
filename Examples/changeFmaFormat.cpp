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
#include "Utils/FGenerateDistribution.hpp"

#include "../Src/Utils/FParameterNames.hpp"

//
/// \file  changeFmaFormat.cpp
//!
//! \brief changeFormat: Driver to transform a FMA format and to build a visualization file
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


int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Driver to change the format of the input file. "
                         "fdlpoly is not supported for now.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::OutputVisuFile, FParameterDefinitions::OutputBinFormat);

    FSize NbPoints;
    FReal	 * particles = nullptr ;
    FFmaGenericLoader * loader;
    unsigned int nbData;

    const std::string filename(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options,   "data.fma"));

    loader = new FFmaGenericLoader(filename);
    //
    // Allocation
    //
    NbPoints                              = loader->getNumberOfParticles();
    nbData   = loader->getNbRecordPerline() ;
    const unsigned int arraySize =nbData*NbPoints;
    //
    particles = new FReal[arraySize] ;
    std::memset(particles,0,arraySize*sizeof(FReal));
    //
    // Read Data
    int j = 0 ;
    for(int idxPart = 0 ; idxPart < NbPoints ;++idxPart, j+=nbData){
        loader->fillParticle(&particles[j],nbData);
    }


    /////////////////////////////////////////////////////////////////////////
    //                                           Save data
    /////////////////////////////////////////////////////////////////////////
    //
    //  Generate file for ScalFMM FMAGenericLoader
    //
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
        std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output"));
        std::string ext(".");
        if(name.find(ext) !=std::string::npos) {
            std::cout << "No file with extension permitted for output name : " << name << std::endl;
            exit(-1);
        }
        if(  FParameters::existParameter(argc, argv, FParameterDefinitions::OutputBinFormat.options)){
            name += ".bfma";
        }
        else {
            name += ".fma";
        }
        FFmaGenericWriter writer(name) ;
        writer.writeHeader( loader->getCenterOfBox(), loader->getBoxWidth() , NbPoints, sizeof(FReal), nbData) ;
        writer.writeArrayOfReal(particles, nbData, NbPoints);
    }
    //
    //   Generate file for visualization purpose
    //
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputVisuFile.options)){
        std::string outfilename(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output"));
        std::string visufile(""), fmt(FParameters::getStr(argc,argv,FParameterDefinitions::OutputVisuFile.options,   "vtp"));
        if( fmt == "vtp" ){
            visufile = outfilename + ".vtp" ;
        }
        else	if( fmt == "vtk" ){
            visufile = outfilename + ".vtk" ;
        }
        else if( fmt == "cosmo" ){
            if(nbData !=4) {
                std::cerr << "Cosmos export accept only 4 data per particles. here: "<<nbData<<std::endl;
                std::exit(EXIT_FAILURE);
            }
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
            if(nbData==4){
                exportVTKxml( file,  particles, NbPoints)   ;
            }
            else {
                exportVTKxml( file,  particles, NbPoints,nbData)   ;
            }
        }
        else		if( fmt == "vtk" ){
            std::cout << "Writes in VTK format  (visualization) in file "<< visufile <<std::endl ;
            exportVTK( file, particles, NbPoints,nbData)  ;
        }
        else if( fmt == "cosmo" ){
            std::cout << "Writes in COSMO format  (visualization) in file "<< visufile <<std::endl ;
            exportCOSMOS( file,  particles, NbPoints)   ;
        }
        else {
            std::cout << "Writes in CVS format  (visualization) in file "<<visufile<<std::endl ;
            exportCVS( file,  particles, NbPoints,nbData)  ;
        }
    }
    //
    delete particles ;

    //
    return 1;
}
