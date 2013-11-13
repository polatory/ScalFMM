#ifndef FCHEBP2PKERNELS_HPP
#define FCHEBP2PKERNELS_HPP


#include "../P2P/FP2P.hpp"

template <KERNEL_FUNCCTION_IDENTIFIER Identifier, int NVALS>
struct DirectInteactionComputer;

///////////////////////////////////////////////////////
// P2P Wrappers
///////////////////////////////////////////////////////

/*! Specialization for Laplace potential */
template <>
struct DirectInteactionComputer<ONE_OVER_R, 1>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
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


/*! Specialization for Leonard-Jones potential */
template <>
struct DirectInteactionComputer<LEONARD_JONES_POTENTIAL, 1>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        FP2P::FullMutualLJ(TargetParticles,NeighborSourceParticles,14);
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        FP2P::FullRemoteLJ(inTargets,inNeighbors,inSize);
    }
};

///////////////////////////////////////////////////////
// In case of multi right hand side
///////////////////////////////////////////////////////



template <int NVALS>
struct DirectInteactionComputer<ONE_OVER_R, NVALS>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14);
        }
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullRemote(inTargets,inNeighbors,inSize);
        }
    }
};


/*! Specialization for Leonard-Jones potential */
template <int NVALS>
struct DirectInteactionComputer<LEONARD_JONES_POTENTIAL, NVALS>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullMutualLJ(TargetParticles,NeighborSourceParticles,14);
        }
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullRemoteLJ(inTargets,inNeighbors,inSize);
        }
    }
};

#endif // FCHEBP2PKERNELS_HPP
