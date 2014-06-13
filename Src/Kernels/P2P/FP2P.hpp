#ifndef FP2P_HPP
#define FP2P_HPP

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"


/**
 * @brief The FP2P namespace
 */
namespace FP2P {

  /** P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = q_1 * q_2 / r^2
   * P_1 = q_2 / r ; P_2 = q_1 / r
   * \f]
   * In details :
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } = \Delta_x * F \f$
   */
  inline void MutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
			      FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
			      const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
			      FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential
			      ){
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
    FReal inv_distance = FMath::Sqrt(inv_square_distance);

    inv_square_distance *= inv_distance;
    inv_square_distance *= targetPhysicalValue * sourcePhysicalValue;

    dx *= inv_square_distance;
    dy *= inv_square_distance;
    dz *= inv_square_distance;

    *targetForceX += dx;
    *targetForceY += dy;
    *targetForceZ += dz;
    *targetPotential += ( inv_distance * sourcePhysicalValue );

    *sourceForceX -= dx;
    *sourceForceY -= dy;
    *sourceForceZ -= dz;
    *sourcePotential += ( inv_distance * targetPhysicalValue );
  }

  /** P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = q_1 * q_2 / r^2
   * P_1 = q_2 / r ; P_2 = q_1 / r
   * \f]
   * In details :
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } = \Delta_x * F \f$
   */
  inline void NonMutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
				 const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
				 FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential){
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
    FReal inv_distance = FMath::Sqrt(inv_square_distance);

    inv_square_distance *= inv_distance;
    inv_square_distance *= targetPhysicalValue * sourcePhysicalValue;

    dx *= inv_square_distance;
    dy *= inv_square_distance;
    dz *= inv_square_distance;

    *targetForceX += dx;
    *targetForceY += dy;
    *targetForceZ += dz;
    *targetPotential += ( inv_distance * sourcePhysicalValue );
  }

  template <class ContainerClass>
  inline void FullMutualLJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
	if( inNeighbors[idxNeighbors] ){
	  const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
	  const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
	  const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
	  const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
	  const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];
	  FReal*const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX();
	  FReal*const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY();
	  FReal*const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ();
	  FReal*const sourcesPotentials = inNeighbors[idxNeighbors]->getPotentials();

	  for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
	    FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
	    FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
	    FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

	    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
	    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
	    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;
	    FReal inv_distance_pow6 = inv_distance_pow3 * inv_distance_pow3;
	    FReal inv_distance_pow8 = inv_distance_pow6 * inv_distance_pow2;

	    FReal coef = ((targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource])
			  * (FReal(12.0)*inv_distance_pow6*inv_distance_pow8 - FReal(6.0)*inv_distance_pow8));
	    FReal potentialCoef = (inv_distance_pow6*inv_distance_pow6-inv_distance_pow6);

	    dx *= coef;
	    dy *= coef;
	    dz *= coef;

	    targetsForcesX[idxTarget] += dx;
	    targetsForcesY[idxTarget] += dy;
	    targetsForcesZ[idxTarget] += dz;
	    targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

	    sourcesForcesX[idxSource] -= dx;
	    sourcesForcesY[idxSource] -= dy;
	    sourcesForcesZ[idxSource] -= dz;
	    sourcesPotentials[idxSource] += potentialCoef * targetsPhysicalValues[idxTarget];
	  }
	}
      }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
      for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){
	FReal dx = targetsX[idxSource] - targetsX[idxTarget];
	FReal dy = targetsY[idxSource] - targetsY[idxTarget];
	FReal dz = targetsZ[idxSource] - targetsZ[idxTarget];

	FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
	FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
	FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;
	FReal inv_distance_pow6 = inv_distance_pow3 * inv_distance_pow3;
	FReal inv_distance_pow8 = inv_distance_pow6 * inv_distance_pow2;

	FReal coef = ((targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource])
		      * (FReal(12.0)*inv_distance_pow6*inv_distance_pow8 - FReal(6.0)*inv_distance_pow8));
	FReal potentialCoef = (inv_distance_pow6*inv_distance_pow6-inv_distance_pow6);

	dx *= coef;
	dy *= coef;
	dz *= coef;

	targetsForcesX[idxTarget] += dx;
	targetsForcesY[idxTarget] += dy;
	targetsForcesZ[idxTarget] += dz;
	targetsPotentials[idxTarget] += ( potentialCoef * targetsPhysicalValues[idxSource] );

	targetsForcesX[idxSource] -= dx;
	targetsForcesY[idxSource] -= dy;
	targetsForcesZ[idxSource] -= dz;
	targetsPotentials[idxSource] += potentialCoef * targetsPhysicalValues[idxTarget];
      }
    }
  }

  /**
   *
   */
  template <class ContainerClass>
  inline void FullRemoteLJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
	if( inNeighbors[idxNeighbors] ){
	  const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
	  const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
	  const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
	  const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
	  const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

	  for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
	    // lenard-jones potential
	    FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
	    FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
	    FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

	    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
	    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
	    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;
	    FReal inv_distance_pow6 = inv_distance_pow3 * inv_distance_pow3;
	    FReal inv_distance_pow8 = inv_distance_pow6 * inv_distance_pow2;

	    FReal coef = ((targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource])
			  * (FReal(12.0)*inv_distance_pow6*inv_distance_pow8 - FReal(6.0)*inv_distance_pow8));

	    dx *= coef;
	    dy *= coef;
	    dz *= coef;

	    targetsForcesX[idxTarget] += dx;
	    targetsForcesY[idxTarget] += dy;
	    targetsForcesZ[idxTarget] += dz;
	    targetsPotentials[idxTarget] += ( (inv_distance_pow6*inv_distance_pow6-inv_distance_pow6) * sourcesPhysicalValues[idxSource] );
	  }
	}
      }
    }
  }


  /**
   * @brief NonMutualParticlesLJ
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   */
  inline void NonMutualParticlesLJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
				   const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
				   FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential){
    // lenard-jones potential
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;
    FReal inv_distance_pow6 = inv_distance_pow3 * inv_distance_pow3;
    FReal inv_distance_pow8 = inv_distance_pow6 * inv_distance_pow2;

    FReal coef = ((targetPhysicalValue * sourcePhysicalValue) * (FReal(12.0)*inv_distance_pow6*inv_distance_pow8
								 - FReal(6.0)*inv_distance_pow8));

    dx *= coef;
    dy *= coef;
    dz *= coef;

    (*targetForceX) += dx;
    (*targetForceY) += dy;
    (*targetForceZ) += dz;
    (*targetPotential) += ( (inv_distance_pow6*inv_distance_pow6-inv_distance_pow6) * sourcePhysicalValue );
  }

  /**
   * @brief MutualParticlesLJ
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param sourceForceX
   * @param sourceForceY
   * @param sourceForceZ
   * @param sourcePotential
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   */
  inline void MutualParticlesLJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
				FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
				const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
				FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential
				){
    // lenard-jones potential
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;
    FReal inv_distance_pow6 = inv_distance_pow3 * inv_distance_pow3;
    FReal inv_distance_pow8 = inv_distance_pow6 * inv_distance_pow2;

    FReal coef = ((targetPhysicalValue * sourcePhysicalValue) * (FReal(12.0)*inv_distance_pow6*inv_distance_pow8
								 - FReal(6.0)*inv_distance_pow8));
    FReal potentialCoef = (inv_distance_pow6*inv_distance_pow6-inv_distance_pow6);

    dx *= coef;
    dy *= coef;
    dz *= coef;

    (*targetForceX) += dx;
    (*targetForceY) += dy;
    (*targetForceZ) += dz;
    (*targetPotential) += ( potentialCoef * sourcePhysicalValue );

    (*sourceForceX) -= dx;
    (*sourceForceY) -= dy;
    (*sourceForceZ) -= dz;
    (*sourcePotential) += ( potentialCoef * targetPhysicalValue );
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // R_IJ
  // PB: The following functions take the MatrixKernel as input arguments for more generic implementation 
  // Besides this MatrixKernel is already updated with any extra parameter (e.g. core width).
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief FullMutualRIJ
   */
  template <class ContainerClass, typename MatrixKernelClass>
  inline void FullMutualRIJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			    const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const double CoreWidth2 = MatrixKernel->getCoreWidth2(); //PB: TODO directly call evaluateBlock

    const int nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
	if( inNeighbors[idxNeighbors] ){
	  const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
	  const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
	  const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
	  const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

	  for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
	    FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
	    FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
	    FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

	    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
	    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
	    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

	    FReal r[3]={dx,dy,dz};

	    for(unsigned int i = 0 ; i < 3 ; ++i){
	      FReal*const targetsPotentials = inTargets->getPotentials(i);
	      FReal*const targetsForcesX = inTargets->getForcesX(i);
	      FReal*const targetsForcesY = inTargets->getForcesY(i);
	      FReal*const targetsForcesZ = inTargets->getForcesZ(i);
	      FReal*const sourcesPotentials = inNeighbors[idxNeighbors]->getPotentials(i);
	      FReal*const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX(i);
	      FReal*const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY(i);
	      FReal*const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ(i);

	      FReal ri2=r[i]*r[i];

	      for(unsigned int j = 0 ; j < 3 ; ++j){
		const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
		const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

		// potentials
		FReal potentialCoef;
		if(i==j)
		  potentialCoef = inv_distance - ri2 * inv_distance_pow3;
		else
		  potentialCoef = - r[i] * r[j] * inv_distance_pow3;

		// forces
		FReal rj2=r[j]*r[j];

		FReal coef[3]; 
		for(unsigned int k = 0 ; k < 3 ; ++k)
		  coef[k]= -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

		// Grad of RIJ kernel is RIJK kernel => use same expression as in FInterpMatrixKernel
		for(unsigned int k = 0 ; k < 3 ; ++k){
		  if(i==j){
		    if(j==k) //i=j=k
		      coef[k] *= FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
		    else //i=j!=k
		      coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
		  }
		  else{ //(i!=j)
		    if(i==k) //i=k!=j
		      coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
		    else if(j==k) //i!=k=j
		      coef[k] *= ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
		    else //i!=k!=j
		      coef[k] *= FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
		  }
		}// k

		targetsForcesX[idxTarget] += coef[0];
		targetsForcesY[idxTarget] += coef[1];
		targetsForcesZ[idxTarget] += coef[2];
		targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

		sourcesForcesX[idxSource] -= coef[0];
		sourcesForcesY[idxSource] -= coef[1];
		sourcesForcesZ[idxSource] -= coef[2];
		sourcesPotentials[idxSource] += potentialCoef * targetsPhysicalValues[idxTarget];

                             
	      }// j
	    }// i
	  }
	}
      }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
      for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){
	FReal dx = targetsX[idxSource] - targetsX[idxTarget];
	FReal dy = targetsY[idxSource] - targetsY[idxTarget];
	FReal dz = targetsZ[idxSource] - targetsZ[idxTarget];

	FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
	FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
	FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

	FReal r[3]={dx,dy,dz};

	for(unsigned int i = 0 ; i < 3 ; ++i){
	  FReal*const targetsPotentials = inTargets->getPotentials(i);
	  FReal*const targetsForcesX = inTargets->getForcesX(i);
	  FReal*const targetsForcesY = inTargets->getForcesY(i);
	  FReal*const targetsForcesZ = inTargets->getForcesZ(i);
	  FReal ri2=r[i]*r[i];

	  for(unsigned int j = 0 ; j < 3 ; ++j){
	    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);

	    // potentials
	    FReal potentialCoef;
	    if(i==j)
	      potentialCoef = inv_distance - ri2 * inv_distance_pow3;
	    else
	      potentialCoef = - r[i] * r[j] * inv_distance_pow3;

	    // forces
	    FReal rj2=r[j]*r[j];

	    FReal coef[3]; 
	    for(unsigned int k = 0 ; k < 3 ; ++k)
	      coef[k]= -(targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);

	    // Grad of RIJ kernel is RIJK kernel => use same expression as in FInterpMatrixKernel
	    for(unsigned int k = 0 ; k < 3 ; ++k){
	      if(i==j){
		if(j==k) //i=j=k
		  coef[k] *= FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
		else //i=j!=k
		  coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
	      }
	      else{ //(i!=j)
		if(i==k) //i=k!=j
		  coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
		else if(j==k) //i!=k=j
		  coef[k] *= ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
		else //i!=k!=j
		  coef[k] *= FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
	      }
	    }// k


	    targetsForcesX[idxTarget] += coef[0];
	    targetsForcesY[idxTarget] += coef[1];
	    targetsForcesZ[idxTarget] += coef[2];
	    targetsPotentials[idxTarget] += ( potentialCoef * targetsPhysicalValues[idxSource] );

	    targetsForcesX[idxSource] -= coef[0];
	    targetsForcesY[idxSource] -= coef[1];
	    targetsForcesZ[idxSource] -= coef[2];
	    targetsPotentials[idxSource] += potentialCoef * targetsPhysicalValues[idxTarget];
	  }// j
	}// i

      }
    }
  }

  /**
   * @brief FullRemoteRIJ
   */
  template <class ContainerClass, typename MatrixKernelClass>
  inline void FullRemoteRIJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			    const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const double CoreWidth2 = MatrixKernel->getCoreWidth2(); //PB: TODO directly call evaluateBlock

    const int nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
	if( inNeighbors[idxNeighbors] ){
	  const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
	  const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
	  const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
	  const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

	  for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
	    // potential
	    FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
	    FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
	    FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

	    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
	    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
	    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

	    FReal r[3]={dx,dy,dz};

	    for(unsigned int i = 0 ; i < 3 ; ++i){
	      FReal*const targetsPotentials = inTargets->getPotentials(i);
	      FReal*const targetsForcesX = inTargets->getForcesX(i);
	      FReal*const targetsForcesY = inTargets->getForcesY(i);
	      FReal*const targetsForcesZ = inTargets->getForcesZ(i);
	      FReal ri2=r[i]*r[i];

	      for(unsigned int j = 0 ; j < 3 ; ++j){
		const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
		const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

		// potentials
		FReal potentialCoef;
		if(i==j)
		  potentialCoef = inv_distance - ri2 * inv_distance_pow3;
		else
		  potentialCoef = - r[i] * r[j] * inv_distance_pow3;

		// forces
		FReal rj2=r[j]*r[j];

		FReal coef[3]; 
		for(unsigned int k = 0 ; k < 3 ; ++k)
		  coef[k]= -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

		// Grad of RIJ kernel is RIJK kernel => use same expression as in FInterpMatrixKernel
		for(unsigned int k = 0 ; k < 3 ; ++k){
		  if(i==j){
		    if(j==k) //i=j=k
		      coef[k] *= FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
		    else //i=j!=k
		      coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
		  }
		  else{ //(i!=j)
		    if(i==k) //i=k!=j
		      coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
		    else if(j==k) //i!=k=j
		      coef[k] *= ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
		    else //i!=k!=j
		      coef[k] *= FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
		  }
		}// k

		targetsForcesX[idxTarget] += coef[0];
		targetsForcesY[idxTarget] += coef[1];
		targetsForcesZ[idxTarget] += coef[2];
		targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

	      }// j
	    }// i

	  }
	}
      }
    }
  }


  /**
   * @brief MutualParticlesRIJ
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param sourceForceX
   * @param sourceForceY
   * @param sourceForceZ
   * @param sourcePotential
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   */
  template<typename MatrixKernelClass>
  inline void MutualParticlesRIJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal* sourcePhysicalValue,
				 FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
				 const FReal targetX,const FReal targetY,const FReal targetZ, const FReal* targetPhysicalValue,
				 FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
				 const MatrixKernelClass *const MatrixKernel){

    const double CoreWidth2 = MatrixKernel->getCoreWidth2(); //PB: TODO directly call evaluateBlock

    // GradGradR potential
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

    FReal r[3]={dx,dy,dz};

    for(unsigned int i = 0 ; i < 3 ; ++i){
      FReal ri2=r[i]*r[i];
      for(unsigned int j = 0 ; j < 3 ; ++j){

	// potentials
	FReal potentialCoef;
	if(i==j)
	  potentialCoef = inv_distance - ri2 * inv_distance_pow3;
	else
	  potentialCoef = - r[i] * r[j] * inv_distance_pow3;

	// forces
	FReal rj2=r[j]*r[j];

	FReal coef[3]; 
	for(unsigned int k = 0 ; k < 3 ; ++k)
	  coef[k]= -(targetPhysicalValue[j] * sourcePhysicalValue[j]);

	// Grad of RIJ kernel is RIJK kernel => use same expression as in FInterpMatrixKernel
	for(unsigned int k = 0 ; k < 3 ; ++k){
	  if(i==j){
	    if(j==k) //i=j=k
	      coef[k] *= FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
	    else //i=j!=k
	      coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
	  }
	  else{ //(i!=j)
	    if(i==k) //i=k!=j
	      coef[k] *= ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
	    else if(j==k) //i!=k=j
	      coef[k] *= ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
	    else //i!=k!=j
	      coef[k] *= FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
	  }
	}// k

	targetForceX[i] += coef[0];
	targetForceY[i] += coef[1];
	targetForceZ[i] += coef[2];
	targetPotential[i] += ( potentialCoef * sourcePhysicalValue[j] );

	sourceForceX[i] -= coef[0];
	sourceForceY[i] -= coef[1];
	sourceForceZ[i] -= coef[2];
	sourcePotential[i] += ( potentialCoef * targetPhysicalValue[j] );

      }// j
    }// i

  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // R_IJK
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief FullMutualRIJK
   */
  template <class ContainerClass, typename MatrixKernelClass>
  inline void FullMutualRIJK(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			     const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const double CoreWidth2 = MatrixKernel->getCoreWidth2(); //PB: TODO directly call evaluateBlock

    const int nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        if( inNeighbors[idxNeighbors] ){
          const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
          const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
          const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
          const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

          for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
            FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
            FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
            FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

            FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
            FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
            FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

            FReal r[3]={dx,dy,dz};

            for(unsigned int i = 0 ; i < 3 ; ++i){
              FReal*const targetsPotentials = inTargets->getPotentials(i);
              FReal*const targetsForcesX = inTargets->getForcesX(i);
              FReal*const targetsForcesY = inTargets->getForcesY(i);
              FReal*const targetsForcesZ = inTargets->getForcesZ(i);
              FReal*const sourcesPotentials = inNeighbors[idxNeighbors]->getPotentials(i);
              FReal*const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX(i);
              FReal*const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY(i);
              FReal*const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ(i);

              const FReal ri2=r[i]*r[i];

              for(unsigned int j = 0 ; j < 3 ; ++j){
                FReal rj2=r[j]*r[j];

                for(unsigned int k = 0 ; k < 3 ; ++k){
                  const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j*3+k);
                  const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j*3+k);

                  // potentials
                  FReal potentialCoef;
                  if(i==j){
                    if(j==k) //i=j=k
                      potentialCoef = FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
                    else //i=j!=k
                      potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
                  }
                  else{ //(i!=j)
                    if(i==k) //i=k!=j
                      potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
                    else if(j==k) //i!=k=j
                      potentialCoef = ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
                    else //i!=k!=j
                      potentialCoef = FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
                  }

                  // forces
                  FReal coef[3]; 
                  for(unsigned int l = 0 ; l < 3 ; ++l){
                    coef[l]= -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

		    //                    // Grad of RIJK kernel is RIJKL kernel => TODO Implement
		    //                    coef[l] *= ;

                  }// l

                  targetsForcesX[idxTarget] += coef[0];
                  targetsForcesY[idxTarget] += coef[1];
                  targetsForcesZ[idxTarget] += coef[2];
                  targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

                  sourcesForcesX[idxSource] -= coef[0];
                  sourcesForcesY[idxSource] -= coef[1];
                  sourcesForcesZ[idxSource] -= coef[2];
                  sourcesPotentials[idxSource] += -potentialCoef * targetsPhysicalValues[idxTarget];

                }// k
              }// j
            }// i
          }
        }
      }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
      for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){
        FReal dx = targetsX[idxSource] - targetsX[idxTarget];
        FReal dy = targetsY[idxSource] - targetsY[idxTarget];
        FReal dz = targetsZ[idxSource] - targetsZ[idxTarget];

        FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
        FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
        FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

        FReal r[3]={dx,dy,dz};

        for(unsigned int i = 0 ; i < 3 ; ++i){
          FReal*const targetsPotentials = inTargets->getPotentials(i);
          FReal*const targetsForcesX = inTargets->getForcesX(i);
          FReal*const targetsForcesY = inTargets->getForcesY(i);
          FReal*const targetsForcesZ = inTargets->getForcesZ(i);
          FReal ri2=r[i]*r[i];

          for(unsigned int j = 0 ; j < 3 ; ++j){
            FReal rj2=r[j]*r[j];

            for(unsigned int k = 0 ; k < 3 ; ++k){
              const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j*3+k);

              // potentials
              FReal potentialCoef;
              if(i==j){
                if(j==k) //i=j=k
                  potentialCoef = FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
                else //i=j!=k
                  potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
              }
              else{ //(i!=j)
                if(i==k) //i=k!=j
                  potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
                else if(j==k) //i!=k=j
                  potentialCoef = ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
                else //i!=k!=j
                  potentialCoef = FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
              }

              // forces
              FReal coef[3]; 
              for(unsigned int l = 0 ; l < 3 ; ++l){
                coef[l]= -(targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);

		//                    // Grad of RIJK kernel is RIJKL kernel => TODO Implement
		//                    coef[l] *= ;

              }// l

              targetsForcesX[idxTarget] += coef[0];
              targetsForcesY[idxTarget] += coef[1];
              targetsForcesZ[idxTarget] += coef[2];
              targetsPotentials[idxTarget] += ( potentialCoef * targetsPhysicalValues[idxSource] );

              targetsForcesX[idxSource] -= coef[0];
              targetsForcesY[idxSource] -= coef[1];
              targetsForcesZ[idxSource] -= coef[2];
              targetsPotentials[idxSource] += -potentialCoef * targetsPhysicalValues[idxTarget];
            }// k
          }// j
        }// i

      }
    }
  }

  /**
   * @brief FullRemoteRIJK
   */
  template <class ContainerClass, typename MatrixKernelClass>
  inline void FullRemoteRIJK(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
			     const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const double CoreWidth2 = MatrixKernel->getCoreWidth2(); //PB: TODO directly call evaluateBlock

    const int nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        if( inNeighbors[idxNeighbors] ){
          const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
          const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
          const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
          const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

          for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
            // potential
            FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
            FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
            FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

            FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
            FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
            FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

            FReal r[3]={dx,dy,dz};

            for(unsigned int i = 0 ; i < 3 ; ++i){
              FReal*const targetsPotentials = inTargets->getPotentials(i);
              FReal*const targetsForcesX = inTargets->getForcesX(i);
              FReal*const targetsForcesY = inTargets->getForcesY(i);
              FReal*const targetsForcesZ = inTargets->getForcesZ(i);
              FReal ri2=r[i]*r[i];

              for(unsigned int j = 0 ; j < 3 ; ++j){
                FReal rj2=r[j]*r[j];

                for(unsigned int k = 0 ; k < 3 ; ++k){

                  const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j*3+k);
                  const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j*3+k);

                  // potentials
                  FReal potentialCoef;
                  if(i==j){
                    if(j==k) //i=j=k
                      potentialCoef = FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
                    else //i=j!=k
                      potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
                  }
                  else{ //(i!=j)
                    if(i==k) //i=k!=j
                      potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
                    else if(j==k) //i!=k=j
                      potentialCoef = ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
                    else //i!=k!=j
                      potentialCoef = FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
                  }

                  // forces
                  FReal coef[3]; 
                  for(unsigned int l = 0 ; l < 3 ; ++l){
                    coef[l]= -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                    //                    // Grad of RIJK kernel is RIJKL kernel => TODO Implement
                    //                    coef[l] *= ;

                  }// l

                  targetsForcesX[idxTarget] += coef[0];
                  targetsForcesY[idxTarget] += coef[1];
                  targetsForcesZ[idxTarget] += coef[2];
                  targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

                }// k
              }// j
            }// i

          }
        }
      }
    }
  }


  /**
   * @brief MutualParticlesRIJK
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param sourceForceX
   * @param sourceForceY
   * @param sourceForceZ
   * @param sourcePotential
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   */
  template<typename MatrixKernelClass>
  inline void MutualParticlesRIJK(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal* sourcePhysicalValue,
                                  FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                                  const FReal targetX,const FReal targetY,const FReal targetZ, const FReal* targetPhysicalValue,
                                  FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                                  const MatrixKernelClass *const MatrixKernel){

    const double CoreWidth2 = MatrixKernel->getCoreWidth2(); //PB: TODO directly call evaluateBlock

    // GradGradR potential
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz + CoreWidth2);
    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

    FReal r[3]={dx,dy,dz};

    for(unsigned int i = 0 ; i < 3 ; ++i){
      FReal ri2=r[i]*r[i];
      for(unsigned int j = 0 ; j < 3 ; ++j){
        FReal rj2=r[j]*r[j];
        for(unsigned int k = 0 ; k < 3 ; ++k){

          // potentials
          FReal potentialCoef;
          if(i==j){
            if(j==k) //i=j=k
              potentialCoef = FReal(3.) * ( FReal(-1.) + ri2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
            else //i=j!=k
              potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[k] * inv_distance_pow3;
          }
          else{ //(i!=j)
            if(i==k) //i=k!=j
              potentialCoef = ( FReal(-1.) + FReal(3.) * ri2 * inv_distance_pow2 ) * r[j] * inv_distance_pow3;
            else if(j==k) //i!=k=j
              potentialCoef = ( FReal(-1.) + FReal(3.) * rj2 * inv_distance_pow2 ) * r[i] * inv_distance_pow3;
            else //i!=k!=j
              potentialCoef = FReal(3.) * r[i] * r[j] * r[k] * inv_distance_pow2 * inv_distance_pow3;
          }

          // forces
          FReal coef[3]; 
          for(unsigned int l = 0 ; l < 3 ; ++l){
            coef[l]= -(targetPhysicalValue[j*3+k] * sourcePhysicalValue[j*3+k]);

	    //                    // Grad of RIJK kernel is RIJKL kernel => TODO Implement
	    //                    coef[l] *= ;

          }// l

          targetForceX[i] += coef[0];
          targetForceY[i] += coef[1];
          targetForceZ[i] += coef[2];
          targetPotential[i] += ( potentialCoef * sourcePhysicalValue[j*3+k] );

          sourceForceX[i] -= coef[0];
          sourceForceY[i] -= coef[1];
          sourceForceZ[i] -= coef[2];
          sourcePotential[i] += ( -potentialCoef * targetPhysicalValue[j*3+k] );

        }// k
      }// j
    }// i

  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // PB: Test Tensorial Kernel ID_over_R
  // i.e. [[1 1 1]
  //       [1 1 1] * 1/R
  //       [1 1 1]]
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief FullMutualIOR
   */
  template <class ContainerClass>
  inline void FullMutualIOR(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                            const int limiteNeighbors){

    const int nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        if( inNeighbors[idxNeighbors] ){
          const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
          const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
          const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
          const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

          for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
            FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
            FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
            FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

            FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
            FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
            FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

            FReal r[3]={dx,dy,dz};

            for(unsigned int i = 0 ; i < 3 ; ++i){
              FReal*const targetsPotentials = inTargets->getPotentials(i);
              FReal*const targetsForcesX = inTargets->getForcesX(i);
              FReal*const targetsForcesY = inTargets->getForcesY(i);
              FReal*const targetsForcesZ = inTargets->getForcesZ(i);
              FReal*const sourcesPotentials = inNeighbors[idxNeighbors]->getPotentials(i);
              FReal*const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX(i);
              FReal*const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY(i);
              FReal*const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ(i);

              for(unsigned int j = 0 ; j < 3 ; ++j){
                const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
                const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

                // potentials
                FReal potentialCoef;
                //if(i==j)
                potentialCoef = inv_distance;
                //else
                //  potentialCoef = FReal(0.);

                // forces
                FReal coef[3];
                for(unsigned int k = 0 ; k < 3 ; ++k){
                  //if(i==j){
                  coef[k] = + r[k] * inv_distance_pow3 * (targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);
                  //}
                  //else{
                  // coef[k] = FReal(0.);
                  //}
                }// k

                targetsForcesX[idxTarget] += coef[0];
                targetsForcesY[idxTarget] += coef[1];
                targetsForcesZ[idxTarget] += coef[2];
                targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

                sourcesForcesX[idxSource] -= coef[0];
                sourcesForcesY[idxSource] -= coef[1];
                sourcesForcesZ[idxSource] -= coef[2];
                sourcesPotentials[idxSource] += potentialCoef * targetsPhysicalValues[idxTarget];

                             
              }// j
            }// i
          }
        }
      }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
      for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){
        FReal dx = targetsX[idxSource] - targetsX[idxTarget];
        FReal dy = targetsY[idxSource] - targetsY[idxTarget];
        FReal dz = targetsZ[idxSource] - targetsZ[idxTarget];

        FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
        FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
        FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

        FReal r[3]={dx,dy,dz};

        for(unsigned int i = 0 ; i < 3 ; ++i){
          FReal*const targetsPotentials = inTargets->getPotentials(i);
          FReal*const targetsForcesX = inTargets->getForcesX(i);
          FReal*const targetsForcesY = inTargets->getForcesY(i);
          FReal*const targetsForcesZ = inTargets->getForcesZ(i);

          for(unsigned int j = 0 ; j < 3 ; ++j){
            const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);

            // potentials
            FReal potentialCoef;
            //if(i==j)
            potentialCoef = inv_distance;
            //else
            //  potentialCoef = FReal(0.);

            // forces
            FReal coef[3];
            for(unsigned int k = 0 ; k < 3 ; ++k){
              //if(i==j){
              coef[k] = + r[k] * inv_distance_pow3 * (targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);
              //}
              //else{
              //  coef[k] = FReal(0.);
              //}
            }// k

            targetsForcesX[idxTarget] += coef[0];
            targetsForcesY[idxTarget] += coef[1];
            targetsForcesZ[idxTarget] += coef[2];
            targetsPotentials[idxTarget] += ( potentialCoef * targetsPhysicalValues[idxSource] );

            targetsForcesX[idxSource] -= coef[0];
            targetsForcesY[idxSource] -= coef[1];
            targetsForcesZ[idxSource] -= coef[2];
            targetsPotentials[idxSource] += potentialCoef * targetsPhysicalValues[idxTarget];
                             
          }// j
        }// i

      }
    }
  }

  /**
   * @brief FullRemoteIOR
   */
  template <class ContainerClass>
  inline void FullRemoteIOR(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                            const int limiteNeighbors){

    const int nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
      for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        if( inNeighbors[idxNeighbors] ){
          const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
          const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
          const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
          const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

          for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
            // potential
            FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
            FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
            FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

            FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
            FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
            FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

            FReal r[3]={dx,dy,dz};

            for(unsigned int i = 0 ; i < 3 ; ++i){
              FReal*const targetsPotentials = inTargets->getPotentials(i);
              FReal*const targetsForcesX = inTargets->getForcesX(i);
              FReal*const targetsForcesY = inTargets->getForcesY(i);
              FReal*const targetsForcesZ = inTargets->getForcesZ(i);

              for(unsigned int j = 0 ; j < 3 ; ++j){
                const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
                const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

                // potentials
                FReal potentialCoef;
                //if(i==j)
                potentialCoef = inv_distance ;
                //else
                //  potentialCoef = FReal(0.);

                // forces
                FReal coef[3];
                for(unsigned int k = 0 ; k < 3 ; ++k){
                  //if(i==j){
                  coef[k] = + r[k] * inv_distance_pow3 * (targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);
                  //}
                  //else{
                  //  coef[k] = FReal(0.);
                  //}
                }// k

                targetsForcesX[idxTarget] += coef[0];
                targetsForcesY[idxTarget] += coef[1];
                targetsForcesZ[idxTarget] += coef[2];
                targetsPotentials[idxTarget] += ( potentialCoef * sourcesPhysicalValues[idxSource] );

              }// j
            }// i

          }
        }
      }
    }
  }


  /**
   * @brief MutualParticlesIOR
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param sourceForceX
   * @param sourceForceY
   * @param sourceForceZ
   * @param sourcePotential
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   */
  inline void MutualParticlesIOR(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal* sourcePhysicalValue,
                                 FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                                 const FReal targetX,const FReal targetY,const FReal targetZ, const FReal* targetPhysicalValue,
                                 FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential
                                 ){
    // GradGradR potential
    FReal dx = sourceX - targetX;
    FReal dy = sourceY - targetY;
    FReal dz = sourceZ - targetZ;

    FReal inv_distance_pow2 = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
    FReal inv_distance = FMath::Sqrt(inv_distance_pow2);
    FReal inv_distance_pow3 = inv_distance_pow2 * inv_distance;

    FReal r[3]={dx,dy,dz};

    for(unsigned int i = 0 ; i < 3 ; ++i){
      for(unsigned int j = 0 ; j < 3 ; ++j){

        // potentials
        FReal potentialCoef;
        //if(i==j)
        potentialCoef = inv_distance;
        //else
        //  potentialCoef = FReal(0.);

        // forces
        FReal coef[3];
        for(unsigned int k = 0 ; k < 3 ; ++k){
          //if(i==j){
          coef[k] = + r[k] * inv_distance_pow3 * (targetPhysicalValue[j] * sourcePhysicalValue[j]);
          //}
          //else{
          //  coef[k] = FReal(0.);
          //}
        }// k

        targetForceX[i] += coef[0];
        targetForceY[i] += coef[1];
        targetForceZ[i] += coef[2];
        targetPotential[i] += ( potentialCoef * sourcePhysicalValue[j] );

        sourceForceX[i] -= coef[0];
        sourceForceY[i] -= coef[1];
        sourceForceZ[i] -= coef[2];
        sourcePotential[i] += ( potentialCoef * targetPhysicalValue[j] );

      }// j
    }// i
  }

}

#ifdef ScalFMM_USE_SSE
#include "FP2PSse.hpp"
#elif defined(ScalFMM_USE_AVX)
#include "FP2PAvx.hpp"
#else
#include "FP2PClassic.hpp"
#endif //Includes

#endif // FP2P_HPP
