// ===================================================================================
// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
//

//
//
//
// Tests/Release/removeMoment  -Ewald2FMM -fin ../Data/forceNacl_2000_dlpolyPer.bin -dlpoly -fout ../Data/forceNacl_2000.fma


#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include  "ScalFmmConfig.h"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"


#include "Files/FDlpolyLoader.hpp"
#include "Files/FFmaGenericLoader.hpp"


/** DLpoly particle is used in the gadget program
 * here we try to make the same simulation
 */
struct MDParticle {
	FPoint position;
	FReal  forces[3];
	FReal  physicalValue;
	FReal  potential;
	int      index;
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
	//
	FReal scaleEnergy, scaleForce  ;
	const FReal r4pie0 = FReal(138935.4835);

	///////////////////////What we do/////////////////////////////
	if( FParameters::existParameter(argc, argv, "-help") || FParameters::existParameter(argc, argv, "-h") ){
		std::cout << ">> Remove or Add First moment This executable has to be used to compute direct interaction either for periodic or non periodic system.\n";
		std::cout << " Options " << std::endl
				<< "      -fin  fileIn.{bin, txt) input file name "<< std::endl
				<< "      -fout fileout.{bfma, fma)  output file name "<< std::endl
				<< "      -Ewald2FMM to add the Electric Moment in the potential and the force " << std::endl
				<< "      -FMM2Ewald  to remove the  Electric Moment in the potential and the force " << std::endl ;
		std::cout << "     -verbose : print index x y z fx fy fy Q and V" << std::endl;
		exit(-1);
	}
	if( FParameters::existParameter(argc, argv, "-Ewald2FMM") &&FParameters::existParameter(argc, argv, "-FMM2Ewald") ){
		std::cerr << " Only -Ewald2FMM or FDMM2Ewald have to be set " <<std::endl;
		exit(EXIT_FAILURE);
	}
	///
	// ---------------------------------------------------
	//        DL_POLY CONSTANT
	//  ---------------------------------------------------
	scaleEnergy =  1.0;
	scaleForce   = 1.0 ;
	//	bool scale = true ;
	if(FParameters::existParameter(argc, argv, "-dlpoly")){
		//		scale = false ;
		scaleEnergy =  r4pie0 / 418.4 ;   // kcal mol^{-1}
		scaleForce  = -r4pie0 ;           // 10 J mol^{-1} A^{-1}
	}
	std::cout <<"Scaling factor " <<std::endl
			<<      "Energy factor " << scaleEnergy<<std::endl
			<<      "  Force factor " << scaleForce<<std::endl ;
	//////////////////////////////////////////////////////////////
	//
	const std::string    filenameIn (FParameters::getStr(argc,argv,"-fin", "../Data/forceNacl_128_dlpolyPer.bin"));
	const std::string    filenameOut (FParameters::getStr(argc,argv,"-fout", "result.fma"));
	//  file for -saveError option

	FTic counter;

	// -----------------------------------------------------
	//  LOADER EWALD PARTICLES FROM DLPOLY
	//  -----------------------------------------------------
	std::cout << "Opening : " << filenameIn << "\n";
	FDlpolyLoader  *loader = nullptr ;
	std::string ext(".bin");
	// open particle file
	if(filenameIn.find(ext) != std::string::npos) {
		loader  = new FDlpolyBinLoader(filenameIn.c_str());

	}
	else if(filenameIn.find(".txt")!=std::string::npos ) {
		loader  = new FDlpolyAsciiLoader(filenameIn.c_str());
	}
	else  {
		std::cout << "Input file not allowed only .bin or .txt extensions" <<std::endl;
		std::exit ( EXIT_FAILURE) ;
	}

	if(! loader->isOpen()){
		std::cout << "Loader Error, " << filenameIn << " is missing\n";
		exit(EXIT_FAILURE);
	}
	// ---------------------------------------------------------------------------------
	//  Read particles
	// ---------------------------------------------------------------------------------
	double boxsize[3] ;
	boxsize[0] = boxsize[1]= boxsize[2]=loader->getBoxWidth() ;

	FSize numberofParticles = loader->getNumberOfParticles() ;
	std::cout << "Creating & Inserting " << numberofParticles<< " particles ..." << std::endl;
	std::cout << "\tWidth : " << loader->getBoxWidth() << " \t center x : " << loader->getCenterOfBox().getX()
	    											<< " y : " << loader->getCenterOfBox().getY() << " z : " << loader->getCenterOfBox().getZ() << std::endl;


	counter.tic();
	FPoint electricMoment(0.0,0.0,0.0) ;
	// const --> then shared
	MDParticle * const particles = new MDParticle[numberofParticles];
	memset(particles, 0, sizeof(MDParticle) * numberofParticles) ;
	//
	double totalCharge = 0.0;
	//
	for(int idxPart = 0 ; idxPart < numberofParticles; ++idxPart){
		//
		loader->fillParticle(&particles[idxPart].position, particles[idxPart].forces,
				&particles[idxPart].physicalValue,&particles[idxPart].index);
		//
		totalCharge += particles[idxPart].physicalValue ;
		//
		electricMoment.incX(particles[idxPart].physicalValue*particles[idxPart].position.getX() );
		electricMoment.incY(particles[idxPart].physicalValue*particles[idxPart].position.getY() );
		electricMoment.incZ(particles[idxPart].physicalValue*particles[idxPart].position.getZ() );
	}
	counter.tac();

	std::cout << std::endl;
	std::cout << "Total Charge               = "<< totalCharge <<std::endl;
	std::cout << "Electric Moment          = "<< electricMoment <<std::endl;
	std::cout << "Electric Moment norm = "<< electricMoment.norm2()  <<std::endl;
	std::cout << "----------------------------------------------------"<<std::endl;
	std::cout << std::endl;

	std::cout << "Done  " << "(@Reading Ewald file  = " << counter.elapsed() << " s)." << std::endl;
	//
	//
	// ---------------------------------------------------------------------------------
	//   READ DIRECT COMPUTATION
	// ----------------------------------------------------------------
	//  Save direct computation in binary format
	// write binary output file
	//
	// ----------------------------------------------------------------------------
	//	 Remove or Add Electric Moment
	//-----------------------------------------------------------------------------
	//  remove polarization component
	//
	double volume               =  boxsize[0] *boxsize[1]*boxsize[2] ;
	double coeffCorrection  = 2.0*FMath::FPi/volume/3.0 ;
	double preScaleEnergy = 1.0,  postScaleEnergy = 1.0, preScaleForce = 1.0,  postScaleForce = 1.0 ;
	if( FParameters::existParameter(argc, argv, "-Ewald2FMM") ){
		preScaleEnergy = 1.0/scaleEnergy ;  postScaleEnergy = 1.0 ; preScaleForce = 1.0/scaleForce;  postScaleForce = 1.0 ;
	}
	else if( FParameters::existParameter(argc, argv, "-FMM2Ewald") ){
		preScaleEnergy  = 1.0 ;  postScaleEnergy = scaleEnergy ; preScaleForce = 1.0/scaleForce;  postScaleForce = scaleForce ;
		coeffCorrection  = -coeffCorrection;
	}
	else {
		std::cout << " -Ewald2FMM ou -FMM2Ewald should be set"<<std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "coeffCorrection: "<<coeffCorrection<<std::endl;
	std::cout << "preScaleEnergy:  "<<preScaleEnergy<<std::endl;
	std::cout << "postScaleEnergy: "<<postScaleEnergy<<std::endl;
	std::cout << "preScaleForce:     "<<preScaleForce<<std::endl;
	std::cout << "postScaleForce:   "<<postScaleForce<<std::endl;
	//
	double tmp, newEnergy =0.0;
	for(int idx = 0 ; idx < numberofParticles ; ++idx){
		//		std::cout << " Pos " << particles[idx].position.getX()<< "  "<< particles[idx].position.getY()<< "  "<< particles[idx].position.getZ()<< std::endl
		//				      << " F  ori " << particles[idx].forces[0]<< "  "<< particles[idx].forces[1]<< "  "<< particles[idx].forces[2]<< std::endl;

		tmp = particles[idx].position.getX()*electricMoment.getX()  + particles[idx].position.getY()*electricMoment.getY()
								+ particles[idx].position.getZ()*electricMoment.getZ()  ;
		//
		particles[idx].potential += 2.0*tmp*coeffCorrection;
		//
		particles[idx].forces[0] = preScaleForce*particles[idx].forces[0] + 2.0*particles[idx].physicalValue*coeffCorrection*electricMoment.getX() ;
		particles[idx].forces[1] = preScaleForce*particles[idx].forces[1] + 2.0*particles[idx].physicalValue*coeffCorrection*electricMoment.getY() ;
		particles[idx].forces[2] = preScaleForce*particles[idx].forces[2] + 2.0*particles[idx].physicalValue*coeffCorrection*electricMoment.getZ() ;
		//
		particles[idx].forces[0] *= postScaleForce;
		particles[idx].forces[1] *= postScaleForce;
		particles[idx].forces[2] *= postScaleForce;
		//	      std::cout <<  " F  new " << particles[idx].forces[0]<< "  "<< particles[idx].forces[1]<< "  "<< particles[idx].forces[2]<< std::endl<< std::endl;

		//
	}
	std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
	std::cout << std::scientific;
	std::cout.precision(10) ;
	std::cout << " " <<  coeffCorrection*electricMoment.norm2()  <<std::endl;
	std::cout << " " <<  preScaleEnergy*loader->getEnergy()  <<std::endl;
	newEnergy = preScaleEnergy*loader->getEnergy() + coeffCorrection*electricMoment.norm2() ;
	newEnergy     *= postScaleEnergy ;

	//
	std::cout.precision(10) ;
	std::cout << "Energy EWALD  = "<< loader->getEnergy() <<std::endl ;
	std::cout << "Energy New      = "<<newEnergy<<std::endl ;
	//
	// Save data in FMAFormat
	//FSize nbParticles = loader.getNumberOfParticles() ;
	FmaRWParticle<8,8>* const particlesOut = new FmaRWParticle<8,8>[numberofParticles];
	// remove index
	memset(particlesOut, 0, sizeof(FmaRWParticle<8,8>) * numberofParticles) ;
	for(int idx = 0 ; idx < numberofParticles ; ++idx){
		particlesOut[idx].setPosition(particles[idx].position) ;
		particlesOut[idx].setPhysicalValue(particles[idx].physicalValue) ;
		particlesOut[idx].setPotential (particles[idx].potential) ;
		particlesOut[idx].setForces(particles[idx].forces[0],particles[idx].forces[0],particles[idx].forces[2]) ;
	}
	//
	// ----------------------------------------------------------------
	//  Save  computation in binary format
	//
	//

	std::cout << "Generate " << filenameOut <<"  for output file" << std::endl;
	//
//	std::cout << " numberofParticles: " << nbParticles <<"  " << sizeof(numberofParticles) <<std::endl;
//	std::cout << " denergy: " << newEnergy <<"  " << sizeof(newEnergy) <<std::endl;
//	std::cout << " Box size: " << loader.getBoxWidth() << "  " << sizeof(loader.getBoxWidth())<<std::endl;
	//
	FFmaGenericWriter writer(filenameOut) ;
	writer.writeHeader(loader->getCenterOfBox(), loader->getBoxWidth() , numberofParticles,*particlesOut) ;
	writer.writeArrayOfParticles(particlesOut, numberofParticles);

	//
	delete[] particles;

	return 0;
}



