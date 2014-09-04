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

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tensorial Matrix Kernels: K_IJ
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief FullMutualKIJ
 */
  template <class ContainerClass, typename MatrixKernelClass>
  inline void FullMutualKIJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                            const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

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

            // evaluate kernel and its partial derivatives
            const FPoint sourcePoint(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
            const FPoint targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
            FReal Kxy[ncmp];
            FReal dKxy[ncmp][3];
            MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

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

                // update component to be applied
                const int d = applyTab[i*3+j];

                // forces prefactor
                FReal coef = -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                targetsForcesX[idxTarget] += dKxy[d][0] * coef;
                targetsForcesY[idxTarget] += dKxy[d][1] * coef;
                targetsForcesZ[idxTarget] += dKxy[d][2] * coef;
                targetsPotentials[idxTarget] += ( Kxy[d] * sourcesPhysicalValues[idxSource] );

                sourcesForcesX[idxSource] -= dKxy[d][0] * coef;
                sourcesForcesY[idxSource] -= dKxy[d][1] * coef;
                sourcesForcesZ[idxSource] -= dKxy[d][2] * coef;
                sourcesPotentials[idxSource] += Kxy[d] * targetsPhysicalValues[idxTarget];

              }// j
            }// i
          }
        }
      }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
      for(int idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){
 
        // evaluate kernel and its partial derivatives
        const FPoint sourcePoint(targetsX[idxSource],targetsY[idxSource],targetsZ[idxSource]);
        const FPoint targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
        FReal Kxy[ncmp];
        FReal dKxy[ncmp][3];
        MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

        for(unsigned int i = 0 ; i < 3 ; ++i){
          FReal*const targetsPotentials = inTargets->getPotentials(i);
          FReal*const targetsForcesX = inTargets->getForcesX(i);
          FReal*const targetsForcesY = inTargets->getForcesY(i);
          FReal*const targetsForcesZ = inTargets->getForcesZ(i);

          for(unsigned int j = 0 ; j < 3 ; ++j){
            const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);

            // update component to be applied
            const int d = applyTab[i*3+j];

            // forces prefactor
            const FReal coef = -(targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);

            targetsForcesX[idxTarget] += dKxy[d][0] * coef;
            targetsForcesY[idxTarget] += dKxy[d][1] * coef;
            targetsForcesZ[idxTarget] += dKxy[d][2] * coef;
            targetsPotentials[idxTarget] += ( Kxy[d] * targetsPhysicalValues[idxSource] );

            targetsForcesX[idxSource] -= dKxy[d][0] * coef;
            targetsForcesY[idxSource] -= dKxy[d][1] * coef;
            targetsForcesZ[idxSource] -= dKxy[d][2] * coef;
            targetsPotentials[idxSource] += Kxy[d] * targetsPhysicalValues[idxTarget];
          }// j
        }// i

      }
    }
  }

  /**
   * @brief FullRemoteKIJ
   */
  template <class ContainerClass, typename MatrixKernelClass>
  inline void FullRemoteKIJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                            const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

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

            // evaluate kernel and its partial derivatives
            const FPoint sourcePoint(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
            const FPoint targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
            FReal Kxy[ncmp];
            FReal dKxy[ncmp][3];
            MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

            for(unsigned int i = 0 ; i < 3 ; ++i){
              FReal*const targetsPotentials = inTargets->getPotentials(i);
              FReal*const targetsForcesX = inTargets->getForcesX(i);
              FReal*const targetsForcesY = inTargets->getForcesY(i);
              FReal*const targetsForcesZ = inTargets->getForcesZ(i);

              for(unsigned int j = 0 ; j < 3 ; ++j){
                const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
                const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

                // update component to be applied
                const int d = applyTab[i*3+j];

                // forces prefactor
                const FReal coef = -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                targetsForcesX[idxTarget] += dKxy[d][0] * coef;
                targetsForcesY[idxTarget] += dKxy[d][1] * coef;
                targetsForcesZ[idxTarget] += dKxy[d][2] * coef;
                targetsPotentials[idxTarget] += ( Kxy[d] * sourcesPhysicalValues[idxSource] );

              }// j
            }// i

          }
        }
      }
    }
  }

}

#endif // FP2PCLASSIC_HPP
