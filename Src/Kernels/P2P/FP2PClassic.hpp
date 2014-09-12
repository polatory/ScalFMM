// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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
#ifndef FP2PCLASSIC_HPP
#define FP2PCLASSIC_HPP

namespace FP2P {

/*
   * FullMutual (generic version)
   */
template <class ContainerClass, typename MatrixKernelClass>
inline void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                       const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){


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

                    // Compute kernel of interaction and its derivative
                    const FPoint sourcePoint(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
                    const FPoint targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
                    FReal Kxy[1];
                    FReal dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(sourcePoint.getX(),sourcePoint.getY(),sourcePoint.getZ(),
                                                             targetPoint.getX(),targetPoint.getY(),targetPoint.getZ(),Kxy,dKxy);
                    FReal coef = (targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                    targetsForcesX[idxTarget] += dKxy[0] * coef;
                    targetsForcesY[idxTarget] += dKxy[1] * coef;
                    targetsForcesZ[idxTarget] += dKxy[2] * coef;
                    targetsPotentials[idxTarget] += ( Kxy[0] * sourcesPhysicalValues[idxSource] );

                    sourcesForcesX[idxSource] -= dKxy[0] * coef;
                    sourcesForcesY[idxSource] -= dKxy[1] * coef;
                    sourcesForcesZ[idxSource] -= dKxy[2] * coef;
                    sourcesPotentials[idxSource] += ( Kxy[0] * targetsPhysicalValues[idxTarget] );
                }
            }
        }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){

            // Compute kernel of interaction...
            const FPoint sourcePoint(targetsX[idxSource],targetsY[idxSource],targetsZ[idxSource]);
            const FPoint targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
            FReal Kxy[1];
            FReal dKxy[3];
            MatrixKernel->evaluateBlockAndDerivative(sourcePoint.getX(),sourcePoint.getY(),sourcePoint.getZ(),
                                                     targetPoint.getX(),targetPoint.getY(),targetPoint.getZ(),Kxy,dKxy);
            FReal coef = (targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);

            targetsForcesX[idxTarget] += dKxy[0] * coef;
            targetsForcesY[idxTarget] += dKxy[1] * coef;
            targetsForcesZ[idxTarget] += dKxy[2] * coef;
            targetsPotentials[idxTarget] += ( Kxy[0] * targetsPhysicalValues[idxSource] );

            targetsForcesX[idxSource] -= dKxy[0] * coef;
            targetsForcesY[idxSource] -= dKxy[1] * coef;
            targetsForcesZ[idxSource] -= dKxy[2] * coef;
            targetsPotentials[idxSource] += ( Kxy[0] * targetsPhysicalValues[idxTarget] );
        }
    }
}


/**
   * FullRemote (generic version)
   */
template <class ContainerClass, typename MatrixKernelClass>
inline void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                       const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
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

                    // Compute kernel of interaction...
                    const FPoint sourcePoint(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
                    const FPoint targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
                    FReal Kxy[1];
                    FReal dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(sourcePoint.getX(),sourcePoint.getY(),sourcePoint.getZ(),
                                                             targetPoint.getX(),targetPoint.getY(),targetPoint.getZ(),Kxy,dKxy);
                    FReal coef = (targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                    targetsForcesX[idxTarget] += dKxy[0] * coef;
                    targetsForcesY[idxTarget] += dKxy[1] * coef;
                    targetsForcesZ[idxTarget] += dKxy[2] * coef;
                    targetsPotentials[idxTarget] += ( Kxy[0] * sourcesPhysicalValues[idxSource] );
                }
            }
        }
    }
}

}

#endif // FP2PCLASSIC_HPP
