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
#include "Utils/FPoint.hpp"
#include  "Utils/FGenerateDistribution.hpp"
//
/// \file  generateDistributions.cpp
//!
//! \brief generateDistributions: Driver to generate N points (non)uniformly distributed on a given geometry
//!
//! The goal of this driver is to generate uniform or non uniform points on the following geometries
//!
//!   Uniform : cube, cuboid, sphere, prolate,
//!
//!   Non uniform : ellipsoid, prolate
//!
//!  You can set two kind of physical values depending of your problem. By default all values are between 0 and 1.
//!   If you select the argument -charge (see bellow) the values are between -1 and 1.
//!  The arguments available are
//!
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param  -N     The number of points in the distribution (default 20000)
//!     \param   -filename name: generic name for files (without extension) and save data
//!                  with following format in name.fma or name.bfma in -bin is set"
//!      \param  -visufmt format for the visu file (vtk, vtp, cvs or cosmo). vtp is the default
//!      \param -extraLength   value    extra length to add to the boxWidth (default 0.0)
//!
//!  <b> Geometry arguments:</b>
//!      \param  -unitCube uniform distribution on unit cube
//!      \param  -cube uniform distribution on a cube
//!          \arg         -length  R - default value for R is 2.0
//!      \param  -unitSphere uniform distribution on unit sphere
//!      \param  -sphere  uniform distribution on  sphere of radius given by
//!          \arg         -radius  R - default value for R is 2.0
//!        \param   -ellipsoid non uniform distribution on  an ellipsoid of aspect ratio given by
//!              \arg          -ar a:b:c   with a, b and c > 0
//!         \param  -prolate ellipsoid with aspect ratio a:a:c  given by
//!                \arg             -ar a:a:c   with  c > a > 0
//!          \param   -plummer (Highly non uniform) plummer distribution (astrophysics)
//!                   \arg         -radius  R - default value 10.0"
//!
//!
//!  <b> Physical values argument:</b>
//!         \param -charge generate physical values between -1 and 1 otherwise generate between 0 and 1
//!         \param -zeromean  the average of the physical values is zero
//!
//!
//! \b examples
//!
//!   generateDistributions -prolate -ar 2:2:4   -N 20000 -filename prolate
//!
//! or
//!
//!   generateDistributions  -cuboid 2:2:4 -filename cuboid  -visufmt vtp -charge  -zeromean
//!


//
//
void genDistusage() {
	std::cout << "Driver to generate N points (non)uniformly distributed on a given geometry"
			<< std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "   -help       to see the parameters    " << std::endl
			<<     "   -N       The number of points in the distribution    " << std::endl
			<<     "   -extraLength   value    extra length to add to the boxWidth"<< std::endl
			<<    std::endl
			<<     "    Distributions   " << std::endl
			<<     "        Uniform on    " << std::endl
			<<     "             -unitCube  uniform distribution on unit cube" <<std::endl
			<<     "             -cuboid  uniform distribution on rectangular cuboid of size  a:a:c" <<std::endl
			<<     "                     -lengths   a:a:c - default value for R is 1.0:1.0:2.0" <<std::endl
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
			<<     "             -zeromean  the average of the physical values is zero		" <<std::endl<<std::endl
			<<     "     Output " << std::endl
			<<     "             -filename name: generic name for files (without extension) and save data" <<std::endl
			<<     "                     with following format in name.fma or name.bfma in -bin is set" <<std::endl
					<<     "      -bin save output in binary file name.txt" <<std::endl
			<<     "             -visufmt  vtk, vtp, cosmo or cvs format " <<std::endl;
}


int main(int argc, char ** argv){
	//
	if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")|| (argc < 3 )){
		genDistusage() ;
		exit(-1);
	}
	 FReal       extraRadius = 0.000 ;
	const int NbPoints  = FParameters::getValue(argc,argv,"-N",   20000);
	const std::string genericFileName(FParameters::getStr(argc,argv,"-filename",   "unifPointDist"));
	FReal BoxWith = 0.0;
	FPoint Centre(0.0, 0.0,0.0);
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
	if(FParameters::existParameter(argc, argv, "-zeromean")){
		FReal  rm = sum/NbPoints ; sum = 0.0 ;
		j = 3 ;
		for(int i = 0 ; i< NbPoints; ++i, j+=4){
			particles[j]    -= rm ;
			sum              += particles[j]   ;
		}
	}
	std::cout << "Sum physical value "<< sum << "   Mean Value " << sum/NbPoints<<std::endl ;
	//
	// Point  generation
	//
	if(FParameters::existParameter(argc, argv, "-unitCube")){
		unifRandonPointsOnUnitCube(NbPoints, particles) ;
		Centre.setPosition(0.5,0.5,0.5);
		BoxWith = 0.5 ;
	}
	else if(FParameters::existParameter(argc, argv, "-cuboid")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-lengths",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
		unifRandonPointsOnCube(NbPoints, A,B,C,particles) ;
		BoxWith = 0.5*FMath::Max(A,FMath::Max(B,C) );
		FReal halfBW = BoxWith;
		Centre.setPosition(halfBW,halfBW,halfBW);
	}
	else if(FParameters::existParameter(argc, argv, "-unitSphere")){
		unifRandonPointsOnUnitSphere(NbPoints, particles) ;
		BoxWith = 1.0 ;
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
		BoxWith =  2.0*C;
	}
	else if(FParameters::existParameter(argc, argv, "-ellipsoid")){
		std::string  dd(":"),aspectRatio  = FParameters::getStr(argc,argv,"-ar",  "1:1:2");
		FReal A,B,C ;
		size_t pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		pos = aspectRatio.find(":");		aspectRatio.replace(pos,1," ");
		std::stringstream ss(aspectRatio); ss >>A >> B >> C ;
		//		std::cout << "A: "<<A<<" B "<< B << " C: " << C<<std::endl;
		nonunifRandonPointsOnElipsoid(NbPoints,A,B,C,particles);
		BoxWith =  2.0*FMath::Max( A,FMath::Max( B,C)) ;
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
    /////////////////////////////////////////////////////////////////////////
	//                                           Save data
    /////////////////////////////////////////////////////////////////////////


	//
	//  Generate file for ScalFMM Loader
	//
	if(FParameters::existParameter(argc, argv, "-extraLength")){
		extraRadius  = FParameters::getValue(argc,argv,"-extraLength",  0.0);
		BoxWith += 2*extraRadius ;
	}
	if( ! FParameters::existParameter(argc, argv, "-bin")){
		std::ofstream outfileA( genericFileName + ".fma", std::ofstream::out);
		if(!outfileA) {
			std::cout << "Cannot open file."<< std::endl;
			exit(-1)	 ;
		}
		std::cout << "Writes in ascii FMA format  in file "<< genericFileName + ".fma" <<std::endl ;
		std::cout << " Points are in a cube of size  "<< BoxWith << "  Centered in the Origin"<<std::endl;
		//
		outfileA << 	NbPoints << "  " << BoxWith << "  " << Centre.getX() << "  " << Centre.getY() << "  " << Centre.getZ() << std::endl;
		j=0;
		for(int i = 0 ; i< NbPoints; ++i, j+=4){
			outfileA <<    particles[j]    << "       "    <<   particles[j+1]    << "       "   <<   particles[j+2]    << "       "   <<   particles[j+3]   <<std::endl;
		}
	}
	else{
		std::cout << "Writes in binary FMA format  in file "<< genericFileName + ".bfma" <<std::endl ;
		std::cout << " Points are in a cube of size  "<< BoxWith << "  Centered in the Origin"<<std::endl;
		//
		std::fstream outfile(genericFileName + ".bfma",std::ifstream::out| std::ios::binary| std::ios::trunc);
		if(!outfile) {
			std::cout << "Cannot open file."<< std::endl;
			return 1;
		}
		int typeFReal = sizeof( particles[0]) ;
		outfile.write((char* )&typeFReal,sizeof(int));
		outfile.write((char* )const_cast<int*>(&NbPoints),sizeof(NbPoints));
		outfile.write((char*)&BoxWith,sizeof(BoxWith));
		outfile.write((char*)Centre.getDataValue(),sizeof(FReal)*3);
		//
		outfile.write ((char* )&particles[0], 4*sizeof(FReal)*NbPoints);
		outfile.flush();
		//
	}
	//
	//
	//
	if(FParameters::existParameter(argc, argv, "-visufmt")){
		std::string visufile(""), fmt(FParameters::getStr(argc,argv,"-visufmt",   "vtp"));
		if( fmt == "vtp" ){
			visufile = genericFileName + ".vtp" ;
		}
		else	if( fmt == "vtk" ){
			visufile = genericFileName + ".vtk" ;
		}
		else if( fmt == "cosmo" ){
			visufile = genericFileName + ".cosmo" ;
		}
		else {
			visufile =   genericFileName + ".csv" ;
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
			exportVTKxml( file, NbPoints, particles)  ;
		}
		else		if( fmt == "vtk" ){
			std::cout << "Writes in VTK format  (visualization) in file "<< visufile <<std::endl ;
			exportVTK( file, NbPoints, particles)  ;
		}
		else if( fmt == "cosmo" ){
			std::cout << "Writes in COSMO format  (visualization) in file "<< visufile <<std::endl ;
			exportCOSMOS( file, NbPoints, particles)  ;
		}
		else {
			std::cout << "Writes in CVS format  (visualization) in file "<<visufile<<std::endl ;
			exportCVS( file, NbPoints, particles)  ;
		}
	}
	//
	delete particles ;

	//
	return 1;
}
