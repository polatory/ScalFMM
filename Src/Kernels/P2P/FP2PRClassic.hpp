#ifndef FP2PRCLASSIC_HPP
#define FP2PRCLASSIC_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FMath.hpp"


namespace FP2PR{

template <class ContainerClass>
inline void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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

            for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
                    FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
                    FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

                    FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
                    const FReal inv_distance = FMath::Sqrt(inv_square_distance);

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

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const int nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
            const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
            const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
            const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
            const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

            for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){

                for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    FReal dx = sourcesX[idxSource] - targetsX[idxTarget];
                    FReal dy = sourcesY[idxSource] - targetsY[idxTarget];
                    FReal dz = sourcesZ[idxSource] - targetsZ[idxTarget];

                    FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
                    const FReal inv_distance = FMath::Sqrt(inv_square_distance);

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

}
#endif //FP2PRCLASSIC_HPP
