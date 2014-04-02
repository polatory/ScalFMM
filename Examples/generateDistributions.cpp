/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Utils/FMath.hpp"
#include  "Utils/FGenerateDistribution.hpp"
/// \file
//! \example  generateDistributions.cpp
//! \brief generateDistributions: Driver to generate N points (non)uniformly distributed on a geometry

//! Options
//!       -h     to see the parameters available in this driver
//!      -N     The number of points in the distribution
//!    Geometry
//!      \arg     \b -unitSphere uniform distribution on unit sphere
//!      \arg  -sphere  uniform distribution on  sphere of radius given by
//!          \arg         -radius  R - default value for R is 2.0
//!        \arg   -ellipsoid non uniform distribution on  an ellipsoid of aspect ratio given by
//!                  -ar a:b:c   with a, b and c > 0
//!
//!        -prolate ellipsoid with aspect ratio a:a:c
//!                  -ar a:a:c   with  c > a > 0
//!          -plummer (Highly non unuiform) plummer distrinution (astrophysics)"
//!                 -radius  R - default value 10.0"
//!
//!        -filename name: generic name for files (without extension) and save data
//!                 with following format in name.xxx or name.bin in -bin is set
//!        -visu save output in name.txt
//!
//!      \b Physical values
//!        -charge generate physical values between -1 and 1 otherwise generate between 0 and 1
		void genDistusage() {
			std::cout << "Driver to generate N points (non)uniformly distributed on a geometry"
 << std::endl;
           std::cout <<	 "Options  "<< std::endl
		                 <<     "   -h     to see the parameters    " << std::endl
		                 <<     "   -N     The number of points in the distribution    " << std::endl
		                 <<    std::endl
			             <<     "    Distribution    " << std::endl
			             <<     "        Uniform on    " << std::endl
			             <<     "             -unitSphere  uniform distribution on unit sphere" <<std::endl
			             <<     "             -sphere  uniform distribution on  sphere of radius given by" <<std::endl
			             <<     "                     -radius  R - default value for R is 2.0" <<std::endl
			             <<     "             -prolate ellipsoid with aspect ratio a:a:c" <<std::endl
			             <<     "                     -ar a:a:c   with  c > a > 0" <<std::endl<<std::endl
			             <<     "        Non Uniform on    " << std::endl
			             <<     "             -ellipsoid non uniform distribution on  an ellipsoid of aspect ratio given by" <<std::endl
			             <<     "                     -ar a:b:c   with a, b and c > 0" <<std::endl
			             <<     "             -plummer (Highly non unuiform) plummer distrinution (astrophysics)"<<std::endl
			             <<     "                     -radius  R - default value 10.0" <<std::endl
			             <<     "    Physical values" <<std::endl
			             <<     "             -charge generate physical values between -1 and 1 otherwise generate between 0 and 1		" <<std::endl<<std::endl
			             <<     "     Output " << std::endl
			             <<     "             -filename name: generic name for files (without extension) and save data" <<std::endl
			             <<     "                     with following format in name.xxx or name.bin in -bin is set" <<std::endl
			             <<     "             -visu save output in name.txt" <<std::endl;
		}

int main(int argc, char ** argv){
	//
	if(FParameters::existParameter(argc, argv, "-h") ||FParameters::existParameter(argc, argv, "-help") ){
		genDistusage() ;
		exit(-1);
	}
	const FReal       eps = 0.001 ;
	const int NbPoints  = FParameters::getValue(argc,argv,"-N",   20);
	const std::string genericFileName(FParameters::getStr(argc,argv,"-filename",   "unifPointDist"));
	FReal BoxWith ;
	//
	// Allocation
	//
	FReal * particles ;
	particles = new FReal[4*NbPoints] ;
	memset(particles,0,4*NbPoints*sizeof(FReal));
	//
	// Generate physical values
	//

	FReal phyVal, sum,a,b ;
	if(FParameters::existParameter(argc, argv, "-charge")){
		a= 2.0 ; b = -1.0 ;
	}
	else {
		a= 1.0 ; b = 0.0 ;
	}
	sum = 0.0 ;
	int j = 3 ;
	for(int i = 0 ; i< NbPoints; ++i, j+=4){
		phyVal            = a*getRandom() +b  ;
		sum              += phyVal ;
		particles[j]       = phyVal ;
	}
	std::cout << "Sum physical value "<< sum << "   Mean Value " << sum/NbPoints<<std::endl ;
	//
	// Point  generation
	//
	if(FParameters::existParameter(argc, argv, "-unitSphere")){
		unifRandonPointsOnUnitSphere(NbPoints, particles) ;
		BoxWith = 2.0 ;
	}
	else if(FParameters::existParameter(argc, argv, "-sphere")){
		const FReal Radius  = FParameters::getValue(argc,argv,"-radius",  2.0);
		unifRandonPointsOnSphere(NbPoints, Radius,particles) ;
		BoxWith = 2.0*Radius ;
	}
	else if(FParameters::existParameter(argc, argv, "-prolate")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
		if(A != B){
			std::cerr << " A /= B in prolate sllipsoide A =B. Your aspect ratio: "<< aspectRatio<<std::endl;
		}
			std::cout << "A: "<<A<<" B "<< B << " C: " << C<<std::endl;
		unifRandonPointsOnProlate(NbPoints,A,C,particles);
		BoxWith =  C;
	}
	else if(FParameters::existParameter(argc, argv, "-ellipsoid")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
		//		std::cout << "A: "<<A<<" B "<< B << " C: " << C<<std::endl;
		nonunifRandonPointsOnElipsoid(NbPoints,A,B,C,particles);
		BoxWith =  FMath::Max( A,FMath::Max( B,C)) ;
	}
		else if(FParameters::existParameter(argc, argv, "-plummer")){
			const FReal Radius  = FParameters::getValue(argc,argv,"-radius",  10.0);
			unifRandonPlummer(NbPoints, Radius, sum, particles) ;
			BoxWith = 2.0*Radius ;
		}

		else {
			std::cout << "Bad geometry option"<< std::endl;
			exit(-1) ;
		}
		if(FParameters::existParameter(argc, argv, "-visu")){
			std::ofstream file( genericFileName + ".txt", std::ofstream::out);
			if(!file) {
				std::cout << "Cannot open file."<< std::endl;
				exit(-1)	;
			}	//
			//
			// Export data in cvs format
			//
			std::cout << "Writes in CVS format  (visualization) in file "<< genericFileName + ".txt" <<std::endl ;
			exportCVS( file, NbPoints, particles)  ;
			//
			// Export data in vtk format
			//
		}
		//
		//  Generate file for ScalFMM Loader
		//
		std::ofstream outfile( genericFileName + ".fma", std::ofstream::out);
		if(!outfile) {
			std::cout << "Cannot open file."<< std::endl;
			exit(-1)	 ;
		}
		BoxWith += eps ;
		std::cout << "Writes in FMA format  in file "<< genericFileName + ".fma" <<std::endl ;
		std::cout << " Points are in a cube of size  "<< BoxWith << "  Centered in the Origin"<<std::endl;
		//
		outfile << 	NbPoints << "  " << BoxWith << "   0.0   0.0  0.0 " << std::endl;
		j=0;
		for(int i = 0 ; i< NbPoints; ++i, j+=4){
			outfile <<    particles[j]    << "       "    <<   particles[j+1]    << "       "   <<   particles[j+2]    << "       "   <<   particles[j+3]   <<std::endl;
		}
		//
		delete particles ;

		//
		return 1;

	}
