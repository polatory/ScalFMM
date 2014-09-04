#ifndef FP2PRAVX_HPP
#define FP2PRAVX_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FMath.hpp"

#include "../../Utils/FAvx.hpp"


namespace FP2PR{

#ifdef ScalFMM_USE_DOUBLE_PRECISION
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			   const int limiteNeighbors){

	const int nbParticlesTargets = inTargets->getNbParticles();
	const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
	const FReal*const targetsX = inTargets->getPositions()[0];
	const FReal*const targetsY = inTargets->getPositions()[1];
	const FReal*const targetsZ = inTargets->getPositions()[2];
	FReal*const targetsForcesX = inTargets->getForcesX();
	FReal*const targetsForcesY = inTargets->getForcesY();
	FReal*const targetsForcesZ = inTargets->getForcesZ();
	FReal*const targetsPotentials = inTargets->getPotentials();


	const __m256d mOne = _mm256_set1_pd(1.0);
	//P2P between target cells and others ones
	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+3)/4;
		const __m256d*const sourcesPhysicalValues = (const __m256d*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m256d*const sourcesX = (const __m256d*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m256d*const sourcesY = (const __m256d*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m256d*const sourcesZ = (const __m256d*)inNeighbors[idxNeighbors]->getPositions()[2];
		__m256d*const sourcesForcesX = (__m256d*)inNeighbors[idxNeighbors]->getForcesX();
		__m256d*const sourcesForcesY = (__m256d*)inNeighbors[idxNeighbors]->getForcesY();
		__m256d*const sourcesForcesZ = (__m256d*)inNeighbors[idxNeighbors]->getForcesZ();
		__m256d*const sourcesPotentials = (__m256d*)inNeighbors[idxNeighbors]->getPotentials();

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m256d tx = _mm256_broadcast_sd(&targetsX[idxTarget]);
		    const __m256d ty = _mm256_broadcast_sd(&targetsY[idxTarget]);
		    const __m256d tz = _mm256_broadcast_sd(&targetsZ[idxTarget]);
		    const __m256d tv = _mm256_broadcast_sd(&targetsPhysicalValues[idxTarget]);
		    __m256d  tfx = _mm256_setzero_pd();
		    __m256d  tfy = _mm256_setzero_pd();
		    __m256d  tfz = _mm256_setzero_pd();
		    __m256d  tpo = _mm256_setzero_pd();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m256d dx = sourcesX[idxSource] - tx;
			__m256d dy = sourcesY[idxSource] - ty;
			__m256d dz = sourcesZ[idxSource] - tz;

			__m256d inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m256d inv_distance = _mm256_sqrt_pd(inv_square_distance);

			inv_square_distance *= inv_distance;
			inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

			dx *= inv_square_distance;
			dy *= inv_square_distance;
			dz *= inv_square_distance;

			tfx += dx;
			tfy += dy;
			tfz += dz;
			tpo += inv_distance * sourcesPhysicalValues[idxSource];

			sourcesForcesX[idxSource] -= dx;
			sourcesForcesY[idxSource] -= dy;
			sourcesForcesZ[idxSource] -= dz;
			sourcesPotentials[idxSource] += inv_distance * tv;
		    }

		    __attribute__((aligned(32))) double buffer[4];

		    _mm256_store_pd(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm256_store_pd(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm256_store_pd(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm256_store_pd(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];
		}
	    }
	}
	//P2P mutual (between targets cells and its own parts 4 by 4)

	{
	    const int nbParticlesSources = (nbParticlesTargets+3)/4;
	    const __m256d*const sourcesPhysicalValues = (const __m256d*)targetsPhysicalValues;
	    const __m256d*const sourcesX = (const __m256d*)targetsX;
	    const __m256d*const sourcesY = (const __m256d*)targetsY;
	    const __m256d*const sourcesZ = (const __m256d*)targetsZ;
	    __m256d*const sourcesForcesX = (__m256d*)targetsForcesX;
	    __m256d*const sourcesForcesY = (__m256d*)targetsForcesY;
	    __m256d*const sourcesForcesZ = (__m256d*)targetsForcesZ;
	    __m256d*const sourcesPotentials = (__m256d*)targetsPotentials;

	    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		const __m256d tx = _mm256_broadcast_sd(&targetsX[idxTarget]);
		const __m256d ty = _mm256_broadcast_sd(&targetsY[idxTarget]);
		const __m256d tz = _mm256_broadcast_sd(&targetsZ[idxTarget]);
		const __m256d tv = _mm256_broadcast_sd(&targetsPhysicalValues[idxTarget]);
		__m256d  tfx = _mm256_setzero_pd();
		__m256d  tfy = _mm256_setzero_pd();
		__m256d  tfz = _mm256_setzero_pd();
		__m256d  tpo = _mm256_setzero_pd();

		for(int idxSource = (idxTarget+4)/4 ; idxSource < nbParticlesSources ; ++idxSource){
		    __m256d dx = sourcesX[idxSource] - tx;
		    __m256d dy = sourcesY[idxSource] - ty;
		    __m256d dz = sourcesZ[idxSource] - tz;

		    __m256d inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
		    const __m256d inv_distance = _mm256_sqrt_pd(inv_square_distance);

		    inv_square_distance *= inv_distance;
		    inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

		    dx *= inv_square_distance;
		    dy *= inv_square_distance;
		    dz *= inv_square_distance;

		    tfx += dx;
		    tfy += dy;
		    tfz += dz;
		    tpo += inv_distance * sourcesPhysicalValues[idxSource];

		    sourcesForcesX[idxSource] -= dx;
		    sourcesForcesY[idxSource] -= dy;
		    sourcesForcesZ[idxSource] -= dz;
		    sourcesPotentials[idxSource] += inv_distance * tv;
		}

		__attribute__((aligned(32))) double buffer[4];

		_mm256_store_pd(buffer, tfx);
		targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		_mm256_store_pd(buffer, tfy);
		targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		_mm256_store_pd(buffer, tfz);
		targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		_mm256_store_pd(buffer, tpo);
		targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

	    }
	}


	//P2P with the last ones
	for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; idxTarget++){
	    for(int idxS = 1 ; idxS < 4-(idxTarget%4) ; ++idxS){
		int idxSource = idxTarget + idxS;
		FReal dx = targetsX[idxSource] - targetsX[idxTarget];
		FReal dy = targetsY[idxSource] - targetsY[idxTarget];
		FReal dz = targetsZ[idxSource] - targetsZ[idxTarget];

		FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
		const FReal inv_distance = FMath::Sqrt(inv_square_distance);

		inv_square_distance *= inv_distance;
		inv_square_distance *= targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource];

		dx *= inv_square_distance;
		dy *= inv_square_distance;
		dz *= inv_square_distance;

		targetsForcesX[idxTarget] += dx;
		targetsForcesY[idxTarget] += dy;
		targetsForcesZ[idxTarget] += dz;
		targetsPotentials[idxTarget] += inv_distance * targetsPhysicalValues[idxSource];

		targetsForcesX[idxSource] -= dx;
		targetsForcesY[idxSource] -= dy;
		targetsForcesZ[idxSource] -= dz;
		targetsPotentials[idxSource] += inv_distance * targetsPhysicalValues[idxTarget];
	    }
	}
    }


    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			   const int limiteNeighbors){
	const int nbParticlesTargets = inTargets->getNbParticles();
	const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
	const FReal*const targetsX = inTargets->getPositions()[0];
	const FReal*const targetsY = inTargets->getPositions()[1];
	const FReal*const targetsZ = inTargets->getPositions()[2];
	FReal*const targetsForcesX = inTargets->getForcesX();
	FReal*const targetsForcesY = inTargets->getForcesY();
	FReal*const targetsForcesZ = inTargets->getForcesZ();
	FReal*const targetsPotentials = inTargets->getPotentials();

	const __m256d mOne = _mm256_set1_pd(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+3)/4;
		const __m256d*const sourcesPhysicalValues = (const __m256d*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m256d*const sourcesX = (const __m256d*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m256d*const sourcesY = (const __m256d*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m256d*const sourcesZ = (const __m256d*)inNeighbors[idxNeighbors]->getPositions()[2];

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m256d tx = _mm256_broadcast_sd(&targetsX[idxTarget]);
		    const __m256d ty = _mm256_broadcast_sd(&targetsY[idxTarget]);
		    const __m256d tz = _mm256_broadcast_sd(&targetsZ[idxTarget]);
		    const __m256d tv = _mm256_broadcast_sd(&targetsPhysicalValues[idxTarget]);
		    __m256d  tfx = _mm256_setzero_pd();
		    __m256d  tfy = _mm256_setzero_pd();
		    __m256d  tfz = _mm256_setzero_pd();
		    __m256d  tpo = _mm256_setzero_pd();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m256d dx = sourcesX[idxSource] - tx;
			__m256d dy = sourcesY[idxSource] - ty;
			__m256d dz = sourcesZ[idxSource] - tz;

			__m256d inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m256d inv_distance = _mm256_sqrt_pd(inv_square_distance);

			inv_square_distance *= inv_distance;
			inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

			dx *= inv_square_distance;
			dy *= inv_square_distance;
			dz *= inv_square_distance;

			tfx += dx;
			tfy += dy;
			tfz += dz;
			tpo += inv_distance * sourcesPhysicalValues[idxSource];
		    }

		    __attribute__((aligned(32))) double buffer[4];

		    _mm256_store_pd(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm256_store_pd(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm256_store_pd(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm256_store_pd(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];
		}
	    }
	}
    }

#else //Float, ScalFMM_USE_DOUBLE_PRECISION not set
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			   const int limiteNeighbors){

	const int nbParticlesTargets = inTargets->getNbParticles();
	const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
	const FReal*const targetsX = inTargets->getPositions()[0];
	const FReal*const targetsY = inTargets->getPositions()[1];
	const FReal*const targetsZ = inTargets->getPositions()[2];
	FReal*const targetsForcesX = inTargets->getForcesX();
	FReal*const targetsForcesY = inTargets->getForcesY();
	FReal*const targetsForcesZ = inTargets->getForcesZ();
	FReal*const targetsPotentials = inTargets->getPotentials();

	const __m256 mOne = _mm256_set1_ps(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+7)/8;
		const __m256*const sourcesPhysicalValues = (const __m256*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m256*const sourcesX = (const __m256*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m256*const sourcesY = (const __m256*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m256*const sourcesZ = (const __m256*)inNeighbors[idxNeighbors]->getPositions()[2];
		__m256*const sourcesForcesX = (__m256*)inNeighbors[idxNeighbors]->getForcesX();
		__m256*const sourcesForcesY = (__m256*)inNeighbors[idxNeighbors]->getForcesY();
		__m256*const sourcesForcesZ = (__m256*)inNeighbors[idxNeighbors]->getForcesZ();
		__m256*const sourcesPotentials = (__m256*)inNeighbors[idxNeighbors]->getPotentials();

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m256 tx = _mm256_broadcast_ss(&targetsX[idxTarget]);
		    const __m256 ty = _mm256_broadcast_ss(&targetsY[idxTarget]);
		    const __m256 tz = _mm256_broadcast_ss(&targetsZ[idxTarget]);
		    const __m256 tv = _mm256_broadcast_ss(&targetsPhysicalValues[idxTarget]);
		    __m256  tfx = _mm256_setzero_ps();
		    __m256  tfy = _mm256_setzero_ps();
		    __m256  tfz = _mm256_setzero_ps();
		    __m256  tpo = _mm256_setzero_ps();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m256 dx = sourcesX[idxSource] - tx;
			__m256 dy = sourcesY[idxSource] - ty;
			__m256 dz = sourcesZ[idxSource] - tz;

			__m256 inv_square_distance = mOne/(dx*dx + dy*dy + dz*dz);
			const __m256 inv_distance = _mm256_rsqrt_ps(dx*dx + dy*dy + dz*dz);

			inv_square_distance *= inv_distance;
			inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

			dx *= inv_square_distance;
			dy *= inv_square_distance;
			dz *= inv_square_distance;

			tfx += dx;
			tfy += dy;
			tfz += dz;
			tpo += inv_distance * sourcesPhysicalValues[idxSource];

			sourcesForcesX[idxSource] -= dx;
			sourcesForcesY[idxSource] -= dy;
			sourcesForcesZ[idxSource] -= dz;
			sourcesPotentials[idxSource] += inv_distance * tv;
		    }

		    __attribute__((aligned(32))) float buffer[8];

		    _mm256_store_ps(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		    _mm256_store_ps(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		    _mm256_store_ps(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		    _mm256_store_ps(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];
		}
	    }
	}

	{
	    const int nbParticlesSources = (nbParticlesTargets+7)/8;
	    const __m256*const sourcesPhysicalValues = (const __m256*)targetsPhysicalValues;
	    const __m256*const sourcesX = (const __m256*)targetsX;
	    const __m256*const sourcesY = (const __m256*)targetsY;
	    const __m256*const sourcesZ = (const __m256*)targetsZ;
	    __m256*const sourcesForcesX = (__m256*)targetsForcesX;
	    __m256*const sourcesForcesY = (__m256*)targetsForcesY;
	    __m256*const sourcesForcesZ = (__m256*)targetsForcesZ;
	    __m256*const sourcesPotentials = (__m256*)targetsPotentials;

	    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		const __m256 tx = _mm256_broadcast_ss(&targetsX[idxTarget]);
		const __m256 ty = _mm256_broadcast_ss(&targetsY[idxTarget]);
		const __m256 tz = _mm256_broadcast_ss(&targetsZ[idxTarget]);
		const __m256 tv = _mm256_broadcast_ss(&targetsPhysicalValues[idxTarget]);
		__m256  tfx = _mm256_setzero_ps();
		__m256  tfy = _mm256_setzero_ps();
		__m256  tfz = _mm256_setzero_ps();
		__m256  tpo = _mm256_setzero_ps();

		for(int idxSource = (idxTarget+8)/8 ; idxSource < nbParticlesSources ; ++idxSource){
		    __m256 dx = sourcesX[idxSource] - tx;
		    __m256 dy = sourcesY[idxSource] - ty;
		    __m256 dz = sourcesZ[idxSource] - tz;

		    __m256 square_distance     = (dx*dx + dy*dy + dz*dz);
		    __m256 inv_square_distance = mOne/square_distance;
		    const __m256 inv_distance  = _mm256_rsqrt_ps(square_distance);

		    inv_square_distance *= inv_distance;
		    inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

		    dx *= inv_square_distance;
		    dy *= inv_square_distance;
		    dz *= inv_square_distance;

		    tfx += dx;
		    tfy += dy;
		    tfz += dz;
		    tpo += inv_distance * sourcesPhysicalValues[idxSource];

		    sourcesForcesX[idxSource] -= dx;
		    sourcesForcesY[idxSource] -= dy;
		    sourcesForcesZ[idxSource] -= dz;
		    sourcesPotentials[idxSource] += inv_distance * tv;
		}

		__attribute__((aligned(32))) float buffer[8];

		_mm256_store_ps(buffer, tfx);
		targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		_mm256_store_ps(buffer, tfy);
		targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		_mm256_store_ps(buffer, tfz);
		targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		_mm256_store_ps(buffer, tpo);
		targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];
	    }
	}

	for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
	    for(int idxS = 1 ; idxS < 8-(idxTarget%8) ; ++idxS){
		const int idxSource = idxTarget + idxS;
		FReal dx = targetsX[idxSource] - targetsX[idxTarget];
		FReal dy = targetsY[idxSource] - targetsY[idxTarget];
		FReal dz = targetsZ[idxSource] - targetsZ[idxTarget];

		FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
		const FReal inv_distance = FMath::Sqrt(inv_square_distance);

		inv_square_distance *= inv_distance;
		inv_square_distance *= targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource];

		dx *= inv_square_distance;
		dy *= inv_square_distance;
		dz *= inv_square_distance;

		targetsForcesX[idxTarget] += dx;
		targetsForcesY[idxTarget] += dy;
		targetsForcesZ[idxTarget] += dz;
		targetsPotentials[idxTarget] += inv_distance * targetsPhysicalValues[idxSource];

		targetsForcesX[idxSource] -= dx;
		targetsForcesY[idxSource] -= dy;
		targetsForcesZ[idxSource] -= dz;
		targetsPotentials[idxSource] += inv_distance * targetsPhysicalValues[idxTarget];
	    }
	}

    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			   const int limiteNeighbors){
	const int nbParticlesTargets = inTargets->getNbParticles();
	const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
	const FReal*const targetsX = inTargets->getPositions()[0];
	const FReal*const targetsY = inTargets->getPositions()[1];
	const FReal*const targetsZ = inTargets->getPositions()[2];
	FReal*const targetsForcesX = inTargets->getForcesX();
	FReal*const targetsForcesY = inTargets->getForcesY();
	FReal*const targetsForcesZ = inTargets->getForcesZ();
	FReal*const targetsPotentials = inTargets->getPotentials();

	const __m256 mOne = _mm256_set1_ps(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+7)/8;
		const __m256*const sourcesPhysicalValues = (const __m256*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m256*const sourcesX = (const __m256*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m256*const sourcesY = (const __m256*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m256*const sourcesZ = (const __m256*)inNeighbors[idxNeighbors]->getPositions()[2];

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m256 tx = _mm256_broadcast_ss(&targetsX[idxTarget]);
		    const __m256 ty = _mm256_broadcast_ss(&targetsY[idxTarget]);
		    const __m256 tz = _mm256_broadcast_ss(&targetsZ[idxTarget]);
		    const __m256 tv = _mm256_broadcast_ss(&targetsPhysicalValues[idxTarget]);
		    __m256  tfx = _mm256_setzero_ps();
		    __m256  tfy = _mm256_setzero_ps();
		    __m256  tfz = _mm256_setzero_ps();
		    __m256  tpo = _mm256_setzero_ps();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m256 dx = sourcesX[idxSource] - tx;
			__m256 dy = sourcesY[idxSource] - ty;
			__m256 dz = sourcesZ[idxSource] - tz;

			__m256 inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m256 inv_distance = _mm256_rsqrt_ps(dx*dx + dy*dy + dz*dz);

			inv_square_distance *= inv_distance;
			inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

			dx *= inv_square_distance;
			dy *= inv_square_distance;
			dz *= inv_square_distance;

			tfx += dx;
			tfy += dy;
			tfz += dz;
			tpo += inv_distance * sourcesPhysicalValues[idxSource];
		    }

		    __attribute__((aligned(32))) float buffer[8];

		    _mm256_store_ps(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		    _mm256_store_ps(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		    _mm256_store_ps(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];

		    _mm256_store_ps(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];
		}
	    }
	}
    }

#endif
}
#endif //FP2PRAVX_HPP