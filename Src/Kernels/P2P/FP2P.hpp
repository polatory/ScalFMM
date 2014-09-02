#ifndef FP2P_HPP
#define FP2P_HPP

namespace FP2P {

/**
   * @brief MutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   *
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
   * @param MatrixKernel pointer to an interaction kernel evaluator
   */
template <typename MatrixKernelClass>
inline void MutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                            FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                            const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                            FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                            const MatrixKernelClass *const MatrixKernel){

    // Compute kernel of interaction...
    const FPoint sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint targetPoint(targetX,targetY,targetZ);
    FReal Kxy[1];
    FReal dKxy[3];
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint.getX(),sourcePoint.getY(),sourcePoint.getZ(),
                                             targetPoint.getX(),targetPoint.getY(),targetPoint.getZ(),Kxy,dKxy);
    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX) += dKxy[0] * coef;
    (*targetForceY) += dKxy[1] * coef;
    (*targetForceZ) += dKxy[2] * coef;
    (*targetPotential) += ( Kxy[0] * sourcePhysicalValue );

    (*sourceForceX) -= dKxy[0] * coef;
    (*sourceForceY) -= dKxy[1] * coef;
    (*sourceForceZ) -= dKxy[2] * coef;
    (*sourcePotential) += ( Kxy[0] * targetPhysicalValue );
}

/**
   * @brief NonMutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   */
template <typename MatrixKernelClass>
inline void NonMutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                               const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                               FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                               const MatrixKernelClass *const MatrixKernel){

    // Compute kernel of interaction...
    const FPoint sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint targetPoint(targetX,targetY,targetZ);
    FReal Kxy[1];
    FReal dKxy[3];
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint.getX(),sourcePoint.getY(),sourcePoint.getZ(),
                                             targetPoint.getX(),targetPoint.getY(),targetPoint.getZ(),Kxy,dKxy);
    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX) += dKxy[0] * coef;
    (*targetForceY) += dKxy[1] * coef;
    (*targetForceZ) += dKxy[2] * coef;
    (*targetPotential) += ( Kxy[0] * sourcePhysicalValue );
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// R_IJ
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


}

#ifdef ScalFMM_USE_SSE
#include "FP2PSse.hpp"
#elif defined(ScalFMM_USE_AVX)
#include "FP2PAvx.hpp"
#else
#include "FP2PClassic.hpp"
#endif //Includes

#endif // FP2P_HPP
