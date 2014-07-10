#ifndef FINTERPP2PKERNELS_HPP
#define FINTERPP2PKERNELS_HPP


#include "../P2P/FP2P.hpp"

template <KERNEL_FUNCTION_IDENTIFIER Identifier, int NVALS>
struct DirectInteractionComputer;

///////////////////////////////////////////////////////
// P2P Wrappers
///////////////////////////////////////////////////////

/*! Specialization for Laplace potential */
template <>
struct DirectInteractionComputer<ONE_OVER_R, 1>
{
  template <typename ContainerClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27]){
    FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14);
  }

  template <typename ContainerClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize){
    FP2P::FullRemote(inTargets,inNeighbors,inSize);
  }
};

//! Specialization for Laplace potential on Non uniform domain
template <>
struct DirectInteractionComputer<ONE_OVER_RH, 1>
{
  template <typename ContainerClass>
  static void P2P( ContainerClass* const FRestrict TargetParticles,
                   ContainerClass* const NeighborSourceParticles[27]){
    FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14);
  }
//
  template <typename ContainerClass>
  static void P2PRemote( ContainerClass* const FRestrict inTargets,
                         ContainerClass* const inNeighbors[27],
                         const int inSize){
    FP2P::FullRemote(inTargets,inNeighbors,inSize);
  }
};



/*! Specialization for Lennard-Jones potential */
template <>
struct DirectInteractionComputer<LENNARD_JONES_POTENTIAL, 1>
{
  template <typename ContainerClass>
  static void P2P(ContainerClass* const FRestrict TargetParticles,
                  ContainerClass* const NeighborSourceParticles[27]){
    FP2P::FullMutualLJ(TargetParticles,NeighborSourceParticles,14);
  }

  template <typename ContainerClass>
  static void P2PRemote(ContainerClass* const FRestrict inTargets,
                        ContainerClass* const inNeighbors[27],
                        const int inSize){
    FP2P::FullRemoteLJ(inTargets,inNeighbors,inSize);
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
  template <typename ContainerClass>
  static void P2P(ContainerClass* const FRestrict TargetParticles,
                  ContainerClass* const NeighborSourceParticles[27]){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14);
    }
  }

  template <typename ContainerClass>
  static void P2PRemote(ContainerClass* const FRestrict inTargets,
                        ContainerClass* const inNeighbors[27],
                        const int inSize){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2P::FullRemote(inTargets,inNeighbors,inSize);
    }
  }
};


/*! Specialization for Lennard-Jones potential */
template <int NVALS>
struct DirectInteractionComputer<LENNARD_JONES_POTENTIAL, NVALS>
{
  template <typename ContainerClass>
  static void P2P(ContainerClass* const FRestrict TargetParticles,
                  ContainerClass* const NeighborSourceParticles[27]){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2P::FullMutualLJ(TargetParticles,NeighborSourceParticles,14);
    }
  }

  template <typename ContainerClass>
  static void P2PRemote(ContainerClass* const FRestrict inTargets,
                        ContainerClass* const inNeighbors[27],
                        const int inSize){
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      FP2P::FullRemoteLJ(inTargets,inNeighbors,inSize);
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
