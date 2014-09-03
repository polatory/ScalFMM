#ifndef FINTERPP2PKERNELS_HPP
#define FINTERPP2PKERNELS_HPP


#include "../P2P/FP2P.hpp"
#include "../P2P/FP2PR.hpp"


template <KERNEL_FUNCTION_IDENTIFIER Identifier, int NVALS>
struct DirectInteractionComputer
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    FP2P::FullRemote(inTargets,inNeighbors,inSize,MatrixKernel);
  }
};

///////////////////////////////////////////////////////
// P2P Wrappers
///////////////////////////////////////////////////////

/*! Specialization for Laplace potential */
template <>
struct DirectInteractionComputer<ONE_OVER_R, 1>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const /*MatrixKernel*/){
    FP2PR::FullMutual(TargetParticles,NeighborSourceParticles,14);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const /*MatrixKernel*/){
    FP2PR::FullRemote(inTargets,inNeighbors,inSize);
  }
};


/*! Specialization for GradGradR potential */
template <>
struct DirectInteractionComputer<R_IJ, 1>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    FP2P::FullMutualRIJ(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    FP2P::FullRemoteRIJ(inTargets,inNeighbors,inSize,MatrixKernel);
  }
};

/*! Specialization for GradGradGradR potential */
template <>
struct DirectInteractionComputer<R_IJK, 1>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    FP2P::FullMutualRIJK(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    FP2P::FullRemoteRIJK(inTargets,inNeighbors,inSize,MatrixKernel);
  }
};

///////////////////////////////////////////////////////
// In case of multi right hand side
///////////////////////////////////////////////////////



template <int NVALS>
struct DirectInteractionComputer<ONE_OVER_R, NVALS>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P(ContainerClass* const FRestrict TargetParticles,
                  ContainerClass* const NeighborSourceParticles[27],
                  const MatrixKernelClass *const /*MatrixKernel*/){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2PR::FullMutual(TargetParticles,NeighborSourceParticles,14);
    }
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote(ContainerClass* const FRestrict inTargets,
                        ContainerClass* const inNeighbors[27],
                        const int inSize,
                        const MatrixKernelClass *const /*MatrixKernel*/){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2PR::FullRemote(inTargets,inNeighbors,inSize);
    }
  }
};


/*! Specialization for GradGradR potential */
template <int NVALS>
struct DirectInteractionComputer<R_IJ, NVALS>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
    FP2P::FullMutualRIJ(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
    }
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
    FP2P::FullRemoteRIJ(inTargets,inNeighbors,inSize,MatrixKernel);
    }
  }
};

/*! Specialization for GradGradGradR potential */
template <int NVALS>
struct DirectInteractionComputer<R_IJK, NVALS>
{
  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27],
                   const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
    FP2P::FullMutualRIJK(TargetParticles,NeighborSourceParticles,14,MatrixKernel);
    }
  }

  template <typename ContainerClass, typename MatrixKernelClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize,
                         const MatrixKernelClass *const MatrixKernel){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
    FP2P::FullRemoteRIJK(inTargets,inNeighbors,inSize,MatrixKernel);
    }
  }
};

#endif // FINTERPP2PKERNELS_HPP
