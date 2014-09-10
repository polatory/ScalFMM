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

void genDistusage() {
	std::cout << ">> Remove or Add First moment This executable has to be used to compute direct interaction either for periodic or non periodic system.\n";
	std::cout << " Options " << std::endl
			<< "      -fin  fileIn.{bin, txt) input file name "<< std::endl
			<< "      -fout fileout.{bfma, fma)  output file name "<< std::endl
			<< "      -Ewald2FMM to add the Electric Moment in the potential and the force " << std::endl
			<< "      -FMM2Ewald  to remove the  Electric Moment in the potential and the force " << std::endl ;
	std::cout << "     -verbose : print index x y z fx fy fy Q and V" << std::endl;
	exit(-1);
}
// Simply create particles and try the kernels
int main(int argc, char ** argv){
	//
	FReal scaleEnergy, scaleForce  ;
	const FReal r4pie0 = FReal(138935.4835);
	const double q0  = 1.6021892e-19;
	const double e0  = 8.854187e-12;
	const double  ang = 1.0e-10;
	const double Unsur4pie0 =8.98755179e+9; // 1.0/(4*FMath::FPi*e0); //8.98755179e+9 ;

	std::cout << "Unsur4pie0: " <<   Unsur4pie0 << " 8.98755179e+9  Diff " << Unsur4pie0-8.98755179e+9 <<std::endl;

	if( argc == 1 ){
		genDistusage() ;
	}
	///////////////////////What we do/////////////////////////////
	if( FParameters::existParameter(argc, argv, "-help") || FParameters::existParameter(argc, argv, "-h") ){
		genDistusage() ;
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
	else 	if(FParameters::existParameter(argc, argv, "-stamp")){
		scaleEnergy =  q0*q0*Unsur4pie0 /ang;
		scaleForce  = -scaleEnergy/ang;
	}
	std::cout <<"Scaling factor " <<std::endl
			<<      "         Energy factor " << scaleEnergy<<std::endl
			<<      "         Force factor   " << scaleForce<<std::endl ;
	//////////////////////////////////////////////////////////////
	//
	const std::string    filenameIn (FParameters::getStr(argc,argv,"-fin", "../Data/forceNacl_128_dlpolyPer.bin"));
	const std::string    filenameOut (FParameters::getStr(argc,argv,"-fout", "result.fma"));
	//  file for -saveError option

	FTic counter;

	std::cout << "Opening : " << filenameIn << "\n";
	FDlpolyLoader          *loaderDlpoly = nullptr ;
	FFmaGenericLoader  *loaderFMA    = nullptr ;
	bool useDLPOLY=false,useFMA=false;
	std::string ext(".bin");
	std::string ext1(".fma"),ext2(".bfma");

	// open particle file
	if(filenameIn.find(ext) != std::string::npos) {
		loaderDlpoly  = new FDlpolyBinLoader(filenameIn.c_str());
		useDLPOLY = true ;
	}
	else if(filenameIn.find(".txt")!=std::string::npos ) {
		loaderDlpoly  = new FDlpolyAsciiLoader(filenameIn.c_str());
		useDLPOLY = true ;
	}
	else if(filenameIn.find(ext1)!=std::string::npos ) {
		loaderFMA  = new FFmaGenericLoader(filenameIn.c_str());
		useFMA = true ;
	}  else if(filenameIn.find(ext2)!=std::string::npos ) {
		loaderFMA  = new FFmaGenericLoader(filenameIn.c_str());
		useFMA    = true ;
	}  else
	{
		std::cout << "Input file not allowed only .bin or .txt extensions" <<std::endl;
		std::exit ( EXIT_FAILURE) ;
	}
	double boxsize[3] ;
	FPoint centre ;
	FSize numberofParticles ;
	FmaRWParticle<8,8>*  particlesIn  ; //= new FmaRWParticle<8,8>[numberofParticles];

	if (useDLPOLY) {
		// -----------------------------------------------------
		//  LOADER EWALD PARTICLES FROM DLPOLY
		//  -----------------------------------------------------
		if(! loaderDlpoly->isOpen()){
			std::cout << "Loader Error, " << filenameIn << " is missing\n";
			exit(EXIT_FAILURE);
		}
		boxsize[0] = boxsize[1]= boxsize[2]=loaderDlpoly->getBoxWidth() ;
		centre = 	loaderDlpoly->getCenterOfBox();
		numberofParticles = loaderDlpoly->getNumberOfParticles() ;
		particlesIn  = new FmaRWParticle<8,8>[numberofParticles];
		memset(particlesIn, 0, sizeof( FmaRWParticle<8,8>) * numberofParticles) ;

		for(int idxPart = 0 ; idxPart < numberofParticles; ++idxPart){
			//
			FPoint position;
			FReal  forces[3];
			FReal  physicalValue;
			FReal  potential=0.0;
			int index ;
			loaderDlpoly->fillParticle(&position, forces,&physicalValue,&index);
			particlesIn[index].setPosition(position) ;
			particlesIn[index].setPhysicalValue(physicalValue) ;
			particlesIn[index].setPotential(potential) ;
			particlesIn[index].setForces(forces[0],forces[1],forces[2]) ;
		}
	}
	else if (useFMA) {
		// -----------------------------------------------------
		//  LOADER EWALD PARTICLES FROM STAMP
		//  -----------------------------------------------------
		if(! loaderFMA->isOpen()){
			std::cout << "Loader Error, " << filenameIn << " is missing\n";
			exit(EXIT_FAILURE);
		}
		boxsize[0] = boxsize[1]= boxsize[2] = loaderFMA->getBoxWidth()/ang;

		//		scaleEnergy /=    boxsize[0]  ;
		//		scaleForce     = -scaleEnergy/boxsize[0] ;
		//		boxsize[0] = boxsize[1]= boxsize[2] = 1.0 ;


		numberofParticles                    = loaderFMA->getNumberOfParticles() ;
		particlesIn  = new FmaRWParticle<8,8>[numberofParticles];
		memset(particlesIn, 0, sizeof( FmaRWParticle<8,8>) * numberofParticles) ;

		loaderFMA->fillParticle(particlesIn,numberofParticles);
		centre = loaderFMA->getCenterOfBox();
		centre *= boxsize[0];
		for(int idxPart = 0 ; idxPart < numberofParticles; ++idxPart){
			//
			FPoint PP = particlesIn[idxPart].getPosition() ;
			PP *=  boxsize[0];
			particlesIn[idxPart].setPosition(PP) ;
			//
		}
	}


	counter.tic();
	FPoint electricMoment(0.0,0.0,0.0) ;
	// const --> then shared
	//
	double totalCharge = 0.0;
	//
	for(int idxPart = 0 ; idxPart < numberofParticles; ++idxPart){
		//
		totalCharge += particlesIn[idxPart].getPhysicalValue() ;
		//
		electricMoment.incX(particlesIn[idxPart].getPhysicalValue()*particlesIn[idxPart].getPosition().getX() );
		electricMoment.incY(particlesIn[idxPart].getPhysicalValue()*particlesIn[idxPart].getPosition().getY() );
		electricMoment.incZ(particlesIn[idxPart].getPhysicalValue()*particlesIn[idxPart].getPosition().getZ() );
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
	// ----------------------------------------------------------------------------
	//	 Remove or Add Electric Moment
	//-----------------------------------------------------------------------------
	//  remove polarization component
	//
	double volume               =  boxsize[0] *boxsize[1]*boxsize[2] ;
	double coeffCorrection   = 4.0*FMath::FPi/volume/3.0 ;
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
			tmp = particlesIn[idx].getPosition().getX()*electricMoment.getX()  + particlesIn[idx].getPosition().getY()*electricMoment.getY()
					+ particlesIn[idx].getPosition().getZ()*electricMoment.getZ()  ;
		//
		double P = particlesIn[idx].getPotential();
		particlesIn[idx].setPotential(  P + 2.0*tmp*coeffCorrection);
		//
		double * forces ;
		forces  = particlesIn[idx].getForces() ;
		forces[0] = preScaleForce*forces[0] - coeffCorrection*particlesIn[idx].getPhysicalValue()*electricMoment.getX() ;
		forces[1] = preScaleForce*forces[1] - coeffCorrection*particlesIn[idx].getPhysicalValue()*electricMoment.getY() ;
		forces[2] = preScaleForce*forces[2] - coeffCorrection*particlesIn[idx].getPhysicalValue()*electricMoment.getZ() ;
		//
		forces[0] *= postScaleForce;
		forces[1] *= postScaleForce;
		forces[2] *= postScaleForce;
		particlesIn[idx].setForces(forces[0],forces[1],forces[2]) ;
		//
	}
	std::cout.precision(10) ;

	std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
	std::cout << std::scientific;
	std::cout.precision(10) ;
	if (useDLPOLY){
		std::cout << " " <<  coeffCorrection*electricMoment.norm2()  <<std::endl;
		std::cout << " " <<  preScaleEnergy*loaderDlpoly->getEnergy()  <<std::endl;
		newEnergy = preScaleEnergy*loaderDlpoly->getEnergy() + coeffCorrection*electricMoment.norm2() ;
		newEnergy     *= postScaleEnergy ;


		//
		std::cout << "Energy EWALD  = "<< loaderDlpoly->getEnergy() <<std::endl ;
	}
	std::cout << "Energy New      = "<<newEnergy<<std::endl ;
	//
	// Save data in FMAFormat
	//
	// ----------------------------------------------------------------
	//  Save  computation in binary format
	//
	//

	std::cout << "Generate " << filenameOut <<"  for output file" << std::endl;
	//
	FFmaGenericWriter writer(filenameOut) ;
	writer.writeHeader(centre, boxsize[0] , numberofParticles,*particlesIn) ;
	writer.writeArrayOfParticles(particlesIn, numberofParticles);

	//
	delete[] particlesIn;
	return 0;
}



