#ifndef FP2P_HPP
#define FP2P_HPP

#include "../../Utils/FGlobal.hpp"

/**
 * @brief The FP2P class
 */
class FP2P {
public:
    /**
     *
     */
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
                        double dx = sourcesX[idxSource] - targetsX[idxTarget];
                        double dy = sourcesY[idxSource] - targetsY[idxTarget];
                        double dz = sourcesZ[idxSource] - targetsZ[idxTarget];

                        double inv_square_distance = 1.0 / (dx*dx + dy*dy + dz*dz);
                        const double inv_distance = sqrt(inv_square_distance);

                        inv_square_distance *= inv_distance;
                        inv_square_distance *= targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource];

                        dx *= inv_square_distance;
                        dy *= inv_square_distance;
                        dz *= inv_square_distance;

                        targetsForcesX[idxTarget] += dx;
                        targetsForcesY[idxTarget] += dy;
                        targetsForcesZ[idxTarget] += dz;
                        targetsPotentials[idxTarget] += inv_distance * sourcesPhysicalValues[idxSource];

                        sourcesForcesX[idxSource] -= dx;
                        sourcesForcesY[idxSource] -= dy;
                        sourcesForcesZ[idxSource] -= dz;
                        sourcesPotentials[idxSource] += inv_distance * targetsPhysicalValues[idxTarget];
                    }
                }
            }
        }

        for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){
                double dx = targetsX[idxSource] - targetsX[idxTarget];
                double dy = targetsY[idxSource] - targetsY[idxTarget];
                double dz = targetsZ[idxSource] - targetsZ[idxTarget];

                double inv_square_distance = 1.0 / (dx*dx + dy*dy + dz*dz);
                const double inv_distance = sqrt(inv_square_distance);

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

    /**
     *
     */
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

        for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
            for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                if( inNeighbors[idxNeighbors] ){
                    const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
                    const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
                    const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
                    const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
                    const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

                    for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                        double dx = sourcesX[idxSource] - targetsX[idxTarget];
                        double dy = sourcesY[idxSource] - targetsY[idxTarget];
                        double dz = sourcesZ[idxSource] - targetsZ[idxTarget];

                        double inv_square_distance = 1.0 / (dx*dx + dy*dy + dz*dz);
                        const double inv_distance = sqrt(inv_square_distance);

                        inv_square_distance *= inv_distance;
                        inv_square_distance *= targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource];

                        dx *= inv_square_distance;
                        dy *= inv_square_distance;
                        dz *= inv_square_distance;

                        targetsForcesX[idxTarget] += dx;
                        targetsForcesY[idxTarget] += dy;
                        targetsForcesZ[idxTarget] += dz;
                        targetsPotentials[idxTarget] += inv_distance * sourcesPhysicalValues[idxSource];
                    }
                }
            }
        }
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
    static void MutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
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
    static void NonMutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
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
    static void FullMutualLJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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
    static void FullRemoteLJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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
    static void NonMutualParticlesLJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
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
    static void MutualParticlesLJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
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
};

#endif // FP2P_HPP
