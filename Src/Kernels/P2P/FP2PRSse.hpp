#ifndef FP2PRSSE_HPP
#define FP2PRSSE_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FMath.hpp"

#include "../../Utils/FSse.hpp"


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

	const __m128d mOne = _mm_set1_pd(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+1)/2;
		const __m128d*const sourcesPhysicalValues = (const __m128d*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m128d*const sourcesX = (const __m128d*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m128d*const sourcesY = (const __m128d*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m128d*const sourcesZ = (const __m128d*)inNeighbors[idxNeighbors]->getPositions()[2];
		__m128d*const sourcesForcesX = (__m128d*)inNeighbors[idxNeighbors]->getForcesX();
		__m128d*const sourcesForcesY = (__m128d*)inNeighbors[idxNeighbors]->getForcesY();
		__m128d*const sourcesForcesZ = (__m128d*)inNeighbors[idxNeighbors]->getForcesZ();
		__m128d*const sourcesPotentials = (__m128d*)inNeighbors[idxNeighbors]->getPotentials();

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m128d tx = _mm_load1_pd(&targetsX[idxTarget]);
		    const __m128d ty = _mm_load1_pd(&targetsY[idxTarget]);
		    const __m128d tz = _mm_load1_pd(&targetsZ[idxTarget]);
		    const __m128d tv = _mm_load1_pd(&targetsPhysicalValues[idxTarget]);
		    __m128d  tfx = _mm_setzero_pd();
		    __m128d  tfy = _mm_setzero_pd();
		    __m128d  tfz = _mm_setzero_pd();
		    __m128d  tpo = _mm_setzero_pd();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m128d dx = sourcesX[idxSource] - tx;
			__m128d dy = sourcesY[idxSource] - ty;
			__m128d dz = sourcesZ[idxSource] - tz;

			__m128d inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m128d inv_distance = _mm_sqrt_pd(inv_square_distance);

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

		    __attribute__((aligned(16))) double buffer[2];

		    _mm_store_pd(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1];

		    _mm_store_pd(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1];

		    _mm_store_pd(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1];

		    _mm_store_pd(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1];
		}
	    }
	}

	{//In this part, we compute (vectorially) the interaction
	 //within the target leaf.

	    const int nbParticlesSources = (nbParticlesTargets+1)/2;
	    const __m128d*const sourcesPhysicalValues = (const __m128d*)targetsPhysicalValues;
	    const __m128d*const sourcesX = (const __m128d*)targetsX;
	    const __m128d*const sourcesY = (const __m128d*)targetsY;
	    const __m128d*const sourcesZ = (const __m128d*)targetsZ;
	    __m128d*const sourcesForcesX = (__m128d*)targetsForcesX;
	    __m128d*const sourcesForcesY = (__m128d*)targetsForcesY;
	    __m128d*const sourcesForcesZ = (__m128d*)targetsForcesZ;
	    __m128d*const sourcesPotentials = (__m128d*)targetsPotentials;

	    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		const __m128d tx = _mm_load1_pd(&targetsX[idxTarget]);
		const __m128d ty = _mm_load1_pd(&targetsY[idxTarget]);
		const __m128d tz = _mm_load1_pd(&targetsZ[idxTarget]);
		const __m128d tv = _mm_load1_pd(&targetsPhysicalValues[idxTarget]);
		__m128d  tfx = _mm_setzero_pd();
		__m128d  tfy = _mm_setzero_pd();
		__m128d  tfz = _mm_setzero_pd();
		__m128d  tpo = _mm_setzero_pd();

		for(int idxSource = (idxTarget+2)/2 ; idxSource < nbParticlesSources ; ++idxSource){

		    __m128d dx = sourcesX[idxSource] - tx;
		    __m128d dy = sourcesY[idxSource] - ty;
		    __m128d dz = sourcesZ[idxSource] - tz;
		    __m128d inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
		    const __m128d inv_distance = _mm_sqrt_pd(inv_square_distance);

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

		__attribute__((aligned(16))) double buffer[2];

		_mm_store_pd(buffer, tfx);
		targetsForcesX[idxTarget] += buffer[0] + buffer[1];

		_mm_store_pd(buffer, tfy);
		targetsForcesY[idxTarget] += buffer[0] + buffer[1];

		_mm_store_pd(buffer, tfz);
		targetsForcesZ[idxTarget] += buffer[0] + buffer[1];

		_mm_store_pd(buffer, tpo);
		targetsPotentials[idxTarget] += buffer[0] + buffer[1];
	    }
	}

	for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; idxTarget += 2){
	    const int idxSource = idxTarget + 1;
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

	const __m128d mOne = _mm_set1_pd(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+1)/2;
		const __m128d*const sourcesPhysicalValues = (const __m128d*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m128d*const sourcesX = (const __m128d*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m128d*const sourcesY = (const __m128d*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m128d*const sourcesZ = (const __m128d*)inNeighbors[idxNeighbors]->getPositions()[2];

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m128d tx = _mm_load1_pd(&targetsX[idxTarget]);
		    const __m128d ty = _mm_load1_pd(&targetsY[idxTarget]);
		    const __m128d tz = _mm_load1_pd(&targetsZ[idxTarget]);
		    const __m128d tv = _mm_load1_pd(&targetsPhysicalValues[idxTarget]);
		    __m128d  tfx = _mm_setzero_pd();
		    __m128d  tfy = _mm_setzero_pd();
		    __m128d  tfz = _mm_setzero_pd();
		    __m128d  tpo = _mm_setzero_pd();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m128d dx = sourcesX[idxSource] - tx;
			__m128d dy = sourcesY[idxSource] - ty;
			__m128d dz = sourcesZ[idxSource] - tz;

			__m128d inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m128d inv_distance = _mm_sqrt_pd(inv_square_distance);

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

		    __attribute__((aligned(16))) double buffer[2];

		    _mm_store_pd(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1];

		    _mm_store_pd(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1];

		    _mm_store_pd(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1];

		    _mm_store_pd(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1];
		}
	    }
	}
    }

#else // Float : ScalFMM_USE_DOUBLE_PRECISION not set
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

	const __m128 mOne = _mm_set1_ps(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+3)/4;
		const __m128*const sourcesPhysicalValues = (const __m128*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m128*const sourcesX = (const __m128*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m128*const sourcesY = (const __m128*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m128*const sourcesZ = (const __m128*)inNeighbors[idxNeighbors]->getPositions()[2];
		__m128*const sourcesForcesX = (__m128*)inNeighbors[idxNeighbors]->getForcesX();
		__m128*const sourcesForcesY = (__m128*)inNeighbors[idxNeighbors]->getForcesY();
		__m128*const sourcesForcesZ = (__m128*)inNeighbors[idxNeighbors]->getForcesZ();
		__m128*const sourcesPotentials = (__m128*)inNeighbors[idxNeighbors]->getPotentials();

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m128 tx = _mm_load1_ps(&targetsX[idxTarget]);
		    const __m128 ty = _mm_load1_ps(&targetsY[idxTarget]);
		    const __m128 tz = _mm_load1_ps(&targetsZ[idxTarget]);
		    const __m128 tv = _mm_load1_ps(&targetsPhysicalValues[idxTarget]);
		    __m128  tfx = _mm_setzero_ps();
		    __m128  tfy = _mm_setzero_ps();
		    __m128  tfz = _mm_setzero_ps();
		    __m128  tpo = _mm_setzero_ps();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m128 dx = sourcesX[idxSource] - tx;
			__m128 dy = sourcesY[idxSource] - ty;
			__m128 dz = sourcesZ[idxSource] - tz;

			__m128 inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m128 inv_distance = _mm_sqrt_ps(inv_square_distance);

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

		    __attribute__((aligned(16))) float buffer[4];

		    _mm_store_ps(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm_store_ps(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm_store_ps(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm_store_ps(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];
		}
	    }
	}

	{
	    const int nbParticlesSources = (nbParticlesTargets+3)/4;
	    const __m128*const sourcesPhysicalValues = (const __m128*)targetsPhysicalValues;
	    const __m128*const sourcesX = (const __m128*)targetsX;
	    const __m128*const sourcesY = (const __m128*)targetsY;
	    const __m128*const sourcesZ = (const __m128*)targetsZ;
	    __m128*const sourcesForcesX = (__m128*)targetsForcesX;
	    __m128*const sourcesForcesY = (__m128*)targetsForcesY;
	    __m128*const sourcesForcesZ = (__m128*)targetsForcesZ;
	    __m128*const sourcesPotentials = (__m128*)targetsPotentials;

	    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		const __m128 tx = _mm_load1_ps(&targetsX[idxTarget]);
		const __m128 ty = _mm_load1_ps(&targetsY[idxTarget]);
		const __m128 tz = _mm_load1_ps(&targetsZ[idxTarget]);
		const __m128 tv = _mm_load1_ps(&targetsPhysicalValues[idxTarget]);
		__m128  tfx = _mm_setzero_ps();
		__m128  tfy = _mm_setzero_ps();
		__m128  tfz = _mm_setzero_ps();
		__m128  tpo = _mm_setzero_ps();

		for(int idxSource = (idxTarget+4)/4 ; idxSource < nbParticlesSources ; ++idxSource){
		    __m128 dx = sourcesX[idxSource] - tx;
		    __m128 dy = sourcesY[idxSource] - ty;
		    __m128 dz = sourcesZ[idxSource] - tz;

		    __m128 inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
		    const __m128 inv_distance = _mm_sqrt_ps(inv_square_distance);

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

		__attribute__((aligned(16))) float buffer[4];

		_mm_store_ps(buffer, tfx);
		targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		_mm_store_ps(buffer, tfy);
		targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		_mm_store_ps(buffer, tfz);
		targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		_mm_store_ps(buffer, tpo);
		targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];
	    }
	}

	for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
	    for(int idxS = 1 ; idxS < 4-(idxTarget%4) ; ++idxS){
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

	const __m128 mOne = _mm_set1_ps(1.0);

	for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
	    if( inNeighbors[idxNeighbors] ){
		const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+3)/4;
		const __m128*const sourcesPhysicalValues = (const __m128*)inNeighbors[idxNeighbors]->getPhysicalValues();
		const __m128*const sourcesX = (const __m128*)inNeighbors[idxNeighbors]->getPositions()[0];
		const __m128*const sourcesY = (const __m128*)inNeighbors[idxNeighbors]->getPositions()[1];
		const __m128*const sourcesZ = (const __m128*)inNeighbors[idxNeighbors]->getPositions()[2];

		for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
		    const __m128 tx = _mm_load1_ps(&targetsX[idxTarget]);
		    const __m128 ty = _mm_load1_ps(&targetsY[idxTarget]);
		    const __m128 tz = _mm_load1_ps(&targetsZ[idxTarget]);
		    const __m128 tv = _mm_load1_ps(&targetsPhysicalValues[idxTarget]);
		    __m128  tfx = _mm_setzero_ps();
		    __m128  tfy = _mm_setzero_ps();
		    __m128  tfz = _mm_setzero_ps();
		    __m128  tpo = _mm_setzero_ps();

		    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
			__m128 dx = sourcesX[idxSource] - tx;
			__m128 dy = sourcesY[idxSource] - ty;
			__m128 dz = sourcesZ[idxSource] - tz;

			__m128 inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
			const __m128 inv_distance = _mm_sqrt_ps(inv_square_distance);

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

		    __attribute__((aligned(16))) float buffer[4];

		    _mm_store_ps(buffer, tfx);
		    targetsForcesX[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm_store_ps(buffer, tfy);
		    targetsForcesY[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm_store_ps(buffer, tfz);
		    targetsForcesZ[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];

		    _mm_store_ps(buffer, tpo);
		    targetsPotentials[idxTarget] += buffer[0] + buffer[1] + buffer[2] + buffer[3];
		}
	    }
	}
    }
#endif
}
#endif //FP2PRSSE_HPP
