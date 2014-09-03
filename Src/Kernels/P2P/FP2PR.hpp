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

}

#ifdef ScalFMM_USE_SSE
#include "FP2PRSse.hpp"
#elif defined(ScalFMM_USE_AVX)
#include "FP2PRAvx.hpp"
#else
#include "FP2PRClassic.hpp"
#endif //Includes

#endif // FP2PR_HPP
