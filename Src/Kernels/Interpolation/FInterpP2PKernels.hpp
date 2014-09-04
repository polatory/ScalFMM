#ifndef FINTERPP2PKERNELS_HPP
#define FINTERPP2PKERNELS_HPP


#include "../P2P/FP2P.hpp"
#include "../P2P/FP2PR.hpp"

///////////////////////////////////////////////////////
// P2P Wrappers
///////////////////////////////////////////////////////

template <int NCMP, int NVALS> 
struct DirectInteractionComputer
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2P::FullMutualKIJ(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
    }
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2P::FullRemoteKIJ(inTargets,inNeighbors,inSize,MatrixKernel);
    }
  }
};


/*! Specialization for scalar kernels */
template <int NVALS>
struct DirectInteractionComputer<1, NVALS>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs)
      FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs)
      FP2P::FullRemote(inTargets,inNeighbors,inSize,MatrixKernel);
  }
};

#endif // FINTERPP2PKERNELS_HPP
