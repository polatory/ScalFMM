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
// Tests/Release/DirectEwaldcmp  -bin -fewald ../Data/forceNacl_128_dlpolyPer.bin -fdirect ../Bench/NaCL/Nacl-128-per=6.bin

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include  "ScalFmmConfig.h"
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"


#include "../../Src/Files/FDlpolyLoader.hpp"


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
	if( FParameters::existParameter(argc, argv, "-help")){
		std::cout << ">> This executable has to be used to compute direct interaction either for periodic or non periodic system.\n";
		std::cout << ">> options are -fewald filenameIN (-bin)  -fdirect  filenameOUT \n";
		std::cout << ">> Recommended files : ../Data/EwalTest_Periodic.run ../Data/EwalTest_NoPeriodic.run\n";
		std::cout << " Options " << std::endl;
		std::cout << "     -verbose : print index x y z fx fy fy Q and V" << std::endl;
		exit(-1);

	}
	if(FParameters::existParameter(argc, argv, "-per") &&FParameters::existParameter(argc, argv, "-noper")  ){
		std::cerr <<" Error -per X and -noper are forbidden together " << std::endl;
		exit(-1);
	}
	///
	// ---------------------------------------------------
	//        DL_POLY CONSTANT
	//  ---------------------------------------------------
//	bool scale = true ;
	if(FParameters::existParameter(argc, argv, "-noscale")){
//		scale = false ;
		scaleEnergy =  1.0;   // kcal mol^{-1}
		scaleForce   = 1.0 ;           // 10 J mol^{-1} A^{-1}
	}
	else {
		scaleEnergy =  r4pie0 / 418.4 ;   // kcal mol^{-1}
		scaleForce  = -r4pie0 ;           // 10 J mol^{-1} A^{-1}
	}

	//////////////////////////////////////////////////////////////

	const std::string    filenameEwaldIn (FParameters::getStr(argc,argv,"-fewald", "../Data/forceNacl_128_dlpolyPer.bin"));
	const std::string    filenameDirectIn (FParameters::getStr(argc,argv,"-fdirect", "periodicDirect.bin"));
	//  file for -saveError option

	FTic counter;

	// -----------------------------------------------------
	//  LOADER EWALD PARTICLES
	//  -----------------------------------------------------
	std::cout << "Opening : " << filenameEwaldIn << "\n";
	FDlpolyLoader  *loader = nullptr ;
	if(FParameters::existParameter(argc, argv, "-bin")){
		loader  = new FDlpolyBinLoader(filenameEwaldIn.c_str());
	}
	else {
		loader  = new FDlpolyAsciiLoader(filenameEwaldIn.c_str());
	}

	if(! loader->isOpen()){
		std::cout << "Loader Error, " << filenameEwaldIn << " is missing\n";
		return 1;
	}
	// ---------------------------------------------------------------------------------
	//  Read particles
	// ---------------------------------------------------------------------------------
	std::cout << "Creating & Inserting " << loader->getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tWidth : " << loader->getBoxWidth() << " \t center x : " << loader->getCenterOfBox().getX()
	    					<< " y : " << loader->getCenterOfBox().getY() << " z : " << loader->getCenterOfBox().getZ() << std::endl;


	counter.tic();
	FPoint electricMoment(0.0,0.0,0.0) ;
	// const --> then shared
	MDParticle * const particles = new MDParticle[loader->getNumberOfParticles()];
	memset(particles, 0, sizeof(MDParticle) * loader->getNumberOfParticles()) ;
	//
	MDParticle* particlesDirect = nullptr;
	particlesDirect = new MDParticle[loader->getNumberOfParticles()];
	memset(particlesDirect, 0, sizeof(MDParticle) * loader->getNumberOfParticles()) ;
	//
	double totalCharge = 0.0;
	//
	for(int idxPart = 0 ; idxPart < loader->getNumberOfParticles() ; ++idxPart){
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
	std::cout << "Generate " << filenameDirectIn <<"  from input file" << std::endl;
	std::ifstream fileIn(filenameDirectIn, std::ifstream::in| std::ios::binary);
	if(!fileIn) {
		std::cout << "Cannot open file."<< std::endl;
		return 1;
	}
	int nbParticles;
	//
	counter.tic();
	fileIn.seekg (std::ios::beg);
	fileIn.read( (char*)(&nbParticles), sizeof(int));
	if (loader->getNumberOfParticles() != nbParticles){
		std::cerr << "Not the same number of particles. reads "<<nbParticles << " expected "<<loader->getNumberOfParticles() <<std::endl;
		exit(-1);
	}
	double boxsize[3] ; boxsize[0] = boxsize[1]= boxsize[2]=0.0;
	int PER[4] ;
	double denergy ;

	fileIn.read( (char*)(&boxsize[0]), sizeof(double)*3);
	fileIn.read((char*)&PER[0],sizeof(int)*4);
	fileIn.read( (char*)(&denergy), sizeof(denergy));
	fileIn.read((char*)&particlesDirect[0], sizeof(MDParticle)*nbParticles);
	//
	totalCharge  =  0.0 ; electricMoment.setX(0.0) ;electricMoment.setY(0.0) ;electricMoment.setZ(0.0) ;
	for(int idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
		//
		totalCharge += particlesDirect[idxPart].physicalValue ;
		electricMoment.incX(particlesDirect[idxPart].physicalValue*particlesDirect[idxPart].position.getX() );
		electricMoment.incY(particlesDirect[idxPart].physicalValue*particlesDirect[idxPart].position.getY() );
		electricMoment.incZ(particlesDirect[idxPart].physicalValue*particlesDirect[idxPart].position.getZ() );
	}
	counter.tac();
	std::cout << " Box size: " << boxsize[0] <<"x"<<boxsize[1]<<"x"<<boxsize[2] <<std::endl;
	std::cout << " Periodicity: " << PER[0] <<"   deep "<< PER[1] <<"x"<<PER[2]<<"x"<<PER[3] <<std::endl;

	std::cout << std::endl;
	std::cout << "Total Charge                 = "<< totalCharge <<std::endl;
	std::cout << "Electric Moment            = "<< electricMoment <<std::endl;
	std::cout << "Electric Moment norm^2 = "<< electricMoment.norm2()  <<std::endl;
	std::cout << "----------------------------------------------------"<<std::endl;
	std::cout << std::endl;
	std::cout << "Done  " << "(@Reading direct file  = " << counter.elapsed() << " s)." << std::endl;

	//
	std::cout << " nbParticles: " << nbParticles <<"  " << sizeof(nbParticles) <<std::endl;
	std::cout << " denergy: " << denergy <<"  " << sizeof(denergy) <<std::endl;
	std::cout << " Box size: " << loader->getBoxWidth() << "  " << sizeof(loader->getBoxWidth())<<std::endl;
	//
	// end generate
	// -----------------------------------------------------
	//	  Correct energy, potential and forces
	//-----------------------------------------------------------------------------
	// remove polarization component
	//
	double volume               =  boxsize[0] *boxsize[1]*boxsize[2] ;
	double coeffCorrection  = 2.0*FMath::FPi/volume/3.0 ;
	//
	double tmp, newEnergy =0.0;
	for(int idx = 0 ; idx < nbParticles ; ++idx){
		tmp = particlesDirect[idx].position.getX()*electricMoment.getX()  + particlesDirect[idx].position.getY()*electricMoment.getY()
						+ particlesDirect[idx].position.getZ()*electricMoment.getZ()  ;
		//
		particlesDirect[idx].potential -= 2.0*tmp*coeffCorrection;
		//
		particlesDirect[idx].forces[0] -= 2.0*particlesDirect[idx].physicalValue*coeffCorrection*electricMoment.getX() ;
		particlesDirect[idx].forces[1] -= 2.0*particlesDirect[idx].physicalValue*coeffCorrection*electricMoment.getY() ;
		particlesDirect[idx].forces[2] -= 2.0*particlesDirect[idx].physicalValue*coeffCorrection*electricMoment.getZ() ;
		//
		particlesDirect[idx].forces[0] *= scaleForce;
		particlesDirect[idx].forces[1] *= scaleForce;
		particlesDirect[idx].forces[2] *= scaleForce;
		//
		newEnergy += particlesDirect[idx].physicalValue*particlesDirect[idx].potential  ;
	}
	denergy -= coeffCorrection*electricMoment.norm2() ;
	//
	denergy     *= scaleEnergy ;
	newEnergy *= 0.5*scaleEnergy ;
	//
	//  Comparisons
	//
	FMath::FAccurater fx, fy, fz ;
	for(int idx = 0 ; idx < nbParticles ; ++idx){
		fx.add(particles[idx].forces[0],particlesDirect[idx].forces[0] );
		fy.add(particles[idx].forces[1],particlesDirect[idx].forces[1] );
		fz.add(particles[idx].forces[2],particlesDirect[idx].forces[2] );

		if(FParameters::existParameter(argc, argv, "-verbose")){
			std::cout << ">> index " << particlesDirect[idx].index << std::endl;
			std::cout << " Dx   " << particlesDirect[idx].position.getX() << " y  " << particlesDirect[idx].position.getY() << " z  " << particlesDirect[idx].position.getZ() << std::endl;
			std::cout << " Ex   " << particles[idx].position.getX() << " y  " << particles[idx].position.getY() << " z  " << particles[idx].position.getZ() << std::endl;
			std::cout << " Dfx  " << particlesDirect[idx].forces[0]         << " fy " << particlesDirect[idx].forces[1]       << " fz " << particlesDirect[idx].forces[2] << std::endl;
			std::cout << " Efx  " << particles[idx].forces[0]         << " fy " << particles[idx].forces[1]       << " fz " << particles[idx].forces[2] << std::endl;
			std::cout << " Q   " << particlesDirect[idx].physicalValue   << " V  " << particlesDirect[idx].potential << std::endl;
			std::cout << "\n";
		}
	}
	//
	std::cout << std::scientific;
	printf("Difference between FMM and DiRECT:\n");
	printf("Fx diff is = \n");
	printf("   RMS       %e\n",fx.getRMSError());
	printf("   L2Norm  %e\n",fx.getRelativeL2Norm());
	printf("   InfNorm %e\n",fx.getRelativeInfNorm());
	printf(" Fy diff is = \n");
	printf("   RMS       %e\n",fy.getRMSError());
	printf("   L2Norm  %e\n",fy.getRelativeL2Norm());
	printf("   InfNorm %e\n",fy.getRelativeInfNorm());
	printf(" Fz diff is = \n");
	printf("   RMS       %e\n",fz.getRMSError());
	printf("   L2Norm  %e\n",fz.getRelativeL2Norm());
	printf("   InfNorm %e\n",fz.getRelativeInfNorm());
	double L2error = (fx.getL2Norm()*fx.getL2Norm() + fy.getL2Norm()*fy.getL2Norm()   + fz.getL2Norm()*fz.getL2Norm()  );
	printf("RMS Forcel Error= %e\n",FMath::Sqrt(L2error/nbParticles)) ;;

	std::cout.precision(10) ;
	std::cout << "Energy EWALD  = "<< loader->getEnergy() <<std::endl ;
	std::cout << "Energy Direct   = " << denergy<<std::endl ;
	std::cout << "Energy New      = "<<newEnergy<<std::endl ;

	delete[] particles;
	delete[] particlesDirect;

	return 0;
}



