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
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);
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
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);
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
// Tensorial Matrix Kernels: K_IJ / p_i=\sum_j K_{ij} w_j
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief MutualParticlesKIJ
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
  inline void MutualParticlesKIJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal* sourcePhysicalValue,
                                 FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                                 const FReal targetX,const FReal targetY,const FReal targetZ, const FReal* targetPhysicalValue,
                                 FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                                 const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

    // evaluate kernel and its partial derivatives
    const FPoint sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint targetPoint(targetX,targetY,targetZ);
    FReal Kxy[ncmp];
    FReal dKxy[ncmp][3];
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

    for(unsigned int i = 0 ; i < 3 ; ++i){
      for(unsigned int j = 0 ; j < 3 ; ++j){

        // update component to be applied
        const int d = applyTab[i*3+j];

        // forces prefactor
        const FReal coef = -(targetPhysicalValue[j] * sourcePhysicalValue[j]);

        targetForceX[i] += dKxy[d][0] * coef;
        targetForceY[i] += dKxy[d][1] * coef;
        targetForceZ[i] += dKxy[d][2] * coef;
        targetPotential[i] += ( Kxy[d] * sourcePhysicalValue[j] );

        sourceForceX[i] -= dKxy[d][0] * coef;
        sourceForceY[i] -= dKxy[d][1] * coef;
        sourceForceZ[i] -= dKxy[d][2] * coef;
        sourcePotential[i] += ( Kxy[d] * targetPhysicalValue[j] );

      }// j
    }// i

  }

}

#ifdef ScalFMM_USE_SSE
#include "FP2PSse.hpp"
#elif defined(ScalFMM_USE_AVX)
#include "FP2PAvx.h"
#else
#include "FP2PClassic.hpp"
#endif //Includes

#endif // FP2P_HPP
