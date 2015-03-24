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
#ifndef FP2PR_HPP
#define FP2PR_HPP

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"


/**
 * @brief The FP2PR namespace
 */
namespace FP2PR{
inline void MutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
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

inline void NonMutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
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


template <class ContainerClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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

    const ComputeClass mOne = FMath::One<ComputeClass>();

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+NbFRealInComputeClass-1)/NbFRealInComputeClass;
            const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const ComputeClass*const sourcesX = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[0];
            const ComputeClass*const sourcesY = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[1];
            const ComputeClass*const sourcesZ = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[2];
            ComputeClass*const sourcesForcesX = (ComputeClass*)inNeighbors[idxNeighbors]->getForcesX();
            ComputeClass*const sourcesForcesY = (ComputeClass*)inNeighbors[idxNeighbors]->getForcesY();
            ComputeClass*const sourcesForcesZ = (ComputeClass*)inNeighbors[idxNeighbors]->getForcesZ();
            ComputeClass*const sourcesPotentials = (ComputeClass*)inNeighbors[idxNeighbors]->getPotentials();

            for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
                const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
                const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
                const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
                ComputeClass  tfx = FMath::Zero<ComputeClass>();
                ComputeClass  tfy = FMath::Zero<ComputeClass>();
                ComputeClass  tfz = FMath::Zero<ComputeClass>();
                ComputeClass  tpo = FMath::Zero<ComputeClass>();

                for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    ComputeClass dx = sourcesX[idxSource] - tx;
                    ComputeClass dy = sourcesY[idxSource] - ty;
                    ComputeClass dz = sourcesZ[idxSource] - tz;

                    ComputeClass inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
                    const ComputeClass inv_distance = FMath::Sqrt(inv_square_distance);

                    inv_square_distance *= inv_distance;
                    inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

                    dx *= inv_square_distance;
                    dy *= inv_square_distance;
                    dz *= inv_square_distance;

                    tfx += dx;
                    tfy += dy;
                    tfz += dz;
                    tpo += inv_distance * sourcesPhysicalValues[idxSource];

                    sourcesForcesX[idxSource] -= dx;
                    sourcesForcesY[idxSource] -= dy;
                    sourcesForcesZ[idxSource] -= dz;
                    sourcesPotentials[idxSource] += inv_distance * tv;
                }

                targetsForcesX[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfx);
                targetsForcesY[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfy);
                targetsForcesZ[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfz);
                targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
            }
        }
    }

    {//In this part, we compute (vectorially) the interaction
        //within the target leaf.

        const int nbParticlesSources = (nbParticlesTargets+NbFRealInComputeClass-1)/NbFRealInComputeClass;
        const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)targetsPhysicalValues;
        const ComputeClass*const sourcesX = (const ComputeClass*)targetsX;
        const ComputeClass*const sourcesY = (const ComputeClass*)targetsY;
        const ComputeClass*const sourcesZ = (const ComputeClass*)targetsZ;
        ComputeClass*const sourcesForcesX = (ComputeClass*)targetsForcesX;
        ComputeClass*const sourcesForcesY = (ComputeClass*)targetsForcesY;
        ComputeClass*const sourcesForcesZ = (ComputeClass*)targetsForcesZ;
        ComputeClass*const sourcesPotentials = (ComputeClass*)targetsPotentials;

        for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
            const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
            const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
            const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
            ComputeClass  tfx = FMath::Zero<ComputeClass>();
            ComputeClass  tfy = FMath::Zero<ComputeClass>();
            ComputeClass  tfz = FMath::Zero<ComputeClass>();
            ComputeClass  tpo = FMath::Zero<ComputeClass>();

            for(int idxSource = (idxTarget+NbFRealInComputeClass)/NbFRealInComputeClass ; idxSource < nbParticlesSources ; ++idxSource){

                ComputeClass dx = sourcesX[idxSource] - tx;
                ComputeClass dy = sourcesY[idxSource] - ty;
                ComputeClass dz = sourcesZ[idxSource] - tz;
                ComputeClass inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
                const ComputeClass inv_distance = FMath::Sqrt(inv_square_distance);

                inv_square_distance *= inv_distance;
                inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

                dx *= inv_square_distance;
                dy *= inv_square_distance;
                dz *= inv_square_distance;

                tfx += dx;
                tfy += dy;
                tfz += dz;
                tpo += inv_distance * sourcesPhysicalValues[idxSource];

                sourcesForcesX[idxSource] -= dx;
                sourcesForcesY[idxSource] -= dy;
                sourcesForcesZ[idxSource] -= dz;
                sourcesPotentials[idxSource] += inv_distance * tv;
            }

            targetsForcesX[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfx);
            targetsForcesY[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfy);
            targetsForcesZ[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfz);
            targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
        }
    }

    for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        const int limitForTarget = NbFRealInComputeClass-(idxTarget%NbFRealInComputeClass);
        for(int idxS = 1 ; idxS < limitForTarget ; ++idxS){
            const int idxSource = idxTarget + idxS;
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

template <class ContainerClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
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

    const ComputeClass mOne = FMath::One<ComputeClass>();

    for(int idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const int nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+NbFRealInComputeClass-1)/NbFRealInComputeClass;
            const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const ComputeClass*const sourcesX = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[0];
            const ComputeClass*const sourcesY = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[1];
            const ComputeClass*const sourcesZ = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[2];

            for(int idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
                const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
                const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
                const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
                ComputeClass  tfx = FMath::Zero<ComputeClass>();
                ComputeClass  tfy = FMath::Zero<ComputeClass>();
                ComputeClass  tfz = FMath::Zero<ComputeClass>();
                ComputeClass  tpo = FMath::Zero<ComputeClass>();

                for(int idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    ComputeClass dx = sourcesX[idxSource] - tx;
                    ComputeClass dy = sourcesY[idxSource] - ty;
                    ComputeClass dz = sourcesZ[idxSource] - tz;

                    ComputeClass inv_square_distance = mOne / (dx*dx + dy*dy + dz*dz);
                    const ComputeClass inv_distance = FMath::Sqrt(inv_square_distance);

                    inv_square_distance *= inv_distance;
                    inv_square_distance *= tv * sourcesPhysicalValues[idxSource];

                    dx *= inv_square_distance;
                    dy *= inv_square_distance;
                    dz *= inv_square_distance;

                    tfx += dx;
                    tfy += dy;
                    tfz += dz;
                    tpo += inv_distance * sourcesPhysicalValues[idxSource];
                }

                targetsForcesX[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfx);
                targetsForcesY[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfy);
                targetsForcesZ[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfz);
                targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
            }
        }
    }
}

} // End namespace

template <class ValueClass>
struct FP2PRT{
};

#if defined(ScalFMM_USE_AVX)

template <>
struct FP2PRT<double>{
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors){
        FP2PR::GenericFullMutual<ContainerClass, __m256d, 4>(inTargets, inNeighbors, limiteNeighbors);
    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
               const int limiteNeighbors){
        FP2PR::GenericFullRemote<ContainerClass, __m256d, 4>(inTargets, inNeighbors, limiteNeighbors);
    }
};

template <>
struct FP2PRT<float>{
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors){
        FP2PR::GenericFullMutual<ContainerClass, __m256, 8>(inTargets, inNeighbors, limiteNeighbors);
    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
               const int limiteNeighbors){
        FP2PR::GenericFullRemote<ContainerClass, __m256, 8>(inTargets, inNeighbors, limiteNeighbors);
    }
};

#elif defined(ScalFMM_USE_SSE)
template <>
struct FP2PRT<double>{
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors){
        FP2PR::GenericFullMutual<ContainerClass, __m128d, 2>(inTargets, inNeighbors, limiteNeighbors);
    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
               const int limiteNeighbors){
        FP2PR::GenericFullRemote<ContainerClass, __m128d, 2>(inTargets, inNeighbors, limiteNeighbors);
    }
};

template <>
struct FP2PRT<float>{
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors){
        FP2PR::GenericFullMutual<ContainerClass, __m128, 4>(inTargets, inNeighbors, limiteNeighbors);
    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
               const int limiteNeighbors){
        FP2PR::GenericFullRemote<ContainerClass, __m128, 4>(inTargets, inNeighbors, limiteNeighbors);
    }
};

#else
template <>
struct FP2PRT<double>{
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors){
        FP2PR::GenericFullMutual<ContainerClass, FReal, 1>(inTargets, inNeighbors, limiteNeighbors);
    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
               const int limiteNeighbors){
        FP2PR::GenericFullRemote<ContainerClass, FReal, 1>(inTargets, inNeighbors, limiteNeighbors);
    }
};

template <>
struct FP2PRT<float>{
    template <class ContainerClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors){
        FP2PR::GenericFullMutual<ContainerClass, FReal, 1>(inTargets, inNeighbors, limiteNeighbors);
    }

    template <class ContainerClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
               const int limiteNeighbors){
        FP2PR::GenericFullRemote<ContainerClass, FReal, 1>(inTargets, inNeighbors, limiteNeighbors);
    }
};
#endif




#endif // FP2PR_HPP
