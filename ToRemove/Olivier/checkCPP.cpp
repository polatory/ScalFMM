/*
 * checkCPP.cpp
 *
 *  Created on: 11 juin 2014
 *      Author: coulaud
 */
#include <iostream>
#include <cstdlib>

#include  "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"



class FmaBAsicParticle {
public:
	FReal  x,y,z;
	FReal  physicalValue;

	FReal * getPtrFirstData()
		{return &x ;}
	const  FReal * getPtrFirstData() const
		{return &x ;}
	int getReadDataNumber()
		{ return 4;}
	int getWriteDataNumber() const
		{ return 4;}
	int getWriteDataSize() const
	{ return sizeof(FmaBAsicParticle);}
};


int main(){
	FmaBAsicParticle part, *ppart = nullptr;
	int NbPoints = 4 ;

	FReal * particles ;
	particles = new FReal[4*NbPoints] ;
	memset(particles,0,4*NbPoints*sizeof(FReal));

	particles[4] = 0.1 ; particles[5] = 0.5 ; particles[6] = 0.8 ; particles[7] = 0.15 ;
	ppart = (FmaBAsicParticle*)(&particles[0]);

	std::cout << "Part "<< ppart[1].x << " "<< ppart[1].y << " "<< ppart[1].z << " "<< ppart[1].physicalValue << " "<< std::endl;



	return 1;
}
