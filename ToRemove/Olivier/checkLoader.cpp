// ===================================================================================
// ===================================================================================
// Copyright ScalFmm 2014 INRIA
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
// @author O. Coulaud
//
//
//
// Tests/Release/DirectEwaldcmp  -bin -fewald ../Data/forceNacl_128_dlpolyPer.bin -fdirect ../Bench/NaCL/Nacl-128-per=6.bin

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include  "ScalFmmConfig.h"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"


#include "Files/FFmaGenericLoader.hpp"

// Simply create particles and try the kernels
int main(int argc, char ** argv){
	//
	///////////////////////What we do/////////////////////////////
	if( FParameters::existParameter(argc, argv, "-help") || FParameters::existParameter(argc, argv, "-h") || argc < 3){
		std::cout << ">> This executable check theFMAGenericLoader.\n";
		std::cout << ">> options are -file filename  \n";
		std::cout << ">> Recommended files : ../Data/EwalTest_Periodic.run ../Data/EwalTest_NoPeriodic.run\n";
		std::cout << " Options " << std::endl;
		std::cout << "      -file filename {with extension }  \n";
		std::cout << "     -verbose : print index x y z Q  " << std::endl;
		exit(-1);

	}
	//////////////////////////////////////////////////////////////

	const std::string    filenameIn (FParameters::getStr(argc,argv,"-file", "unitCubeXYZQ100.bfma"));

	FTic counter,counter1;

	// -----------------------------------------------------
	//  LOADER  PARTICLES
	//  -----------------------------------------------------
	std::cout << "Opening : " << filenameIn << "\n";
	FFmaGenericLoader  loader(filenameIn);

	// ---------------------------------------------------------------------------------
	//  Read particles
	// ---------------------------------------------------------------------------------

	int nbPart =  loader.getNumberOfParticles()  ;
	std::cout << "Reading" <<nbPart << " particles ..." << std::endl;
	std::cout << "\tWidth : " << loader.getBoxWidth() << " \t center x : " << loader.getCenterOfBox().getX()
	    													<< " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;


	//	// const --> then shared
	//	FmaBasicParticle * const particles = new FmaBasicParticle[nbPart];
	FmaParticle *  particles = new FmaParticle[nbPart];
	if(  FParameters::existParameter(argc, argv, "-d")){

		counter.tic();
		loader.fillParticle(particles,nbPart);
		counter.tac();
		//
	}	//
	else {
		counter.tic();

		for(int idxPart = 0 ; idxPart < nbPart ;++idxPart){
			loader.fillParticle(particles[idxPart]);
		}
		counter.tac();
	}
	std::cout << "Done  " << "(@Reading  file one bloc = " << counter.elapsed() << " s)." << std::endl;

	//	std::cout << "i "<<idxPart<<" "<< particles[idxPart].position << " "<< particles[idxPart].physicalValue <<std::endl;

	std::cout << std::endl;
	std::cout << "----------------------------------------------------"<<std::endl;
	std::cout << std::endl;

	std::cout << "Done  " << "(@Reading  file  = " << counter.elapsed() << " s)." << std::endl;
	//

	delete[] particles;
	return 0;
}



