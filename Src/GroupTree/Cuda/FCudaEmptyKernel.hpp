/// @SCALFMM_PRIVATE
#ifndef FCUDAEMPTYKERNEL_HPP
#define FCUDAEMPTYKERNEL_HPP

#include "FCudaGlobal.hpp"
#include "FCudaGroupAttachedLeaf.hpp"

/**
 * This class defines what should be a Cuda kernel.
 */
template <class ContainerClass = FCudaGroupAttachedLeaf<0, int>>
class FCudaEmptyKernel {
protected:
public:
    __device__ void P2M(unsigned char* const /*pole*/, const ContainerClass* const /*particles*/) {
    }

    __device__ void M2M(unsigned char* const  /*pole*/, const unsigned char *const  *const  /*child*/, const int /*level*/) {
    }

    __device__ void M2L(unsigned char* const  /*pole*/, const unsigned char* /*distantNeighbors*/[343],
        const int /*size*/, const int /*level*/) {
    }

    __device__ void L2L(const unsigned char*const  /*local*/, unsigned char*  *const  /*child*/, const int /*level*/) {
    }

    __device__ void L2P(const unsigned char* const  /*local*/, ContainerClass*const /*particles*/){
    }

    __device__ void P2P(const int3& ,
                 ContainerClass* const  /*targets*/, const ContainerClass* const  /*sources*/,
                 ContainerClass* const /*directNeighborsParticles*/[27], const int ){
    }

    __device__ void P2PRemote(const int3& ,
                 ContainerClass* const  /*targets*/, const ContainerClass* const  /*sources*/,
                 ContainerClass* const /*directNeighborsParticles*/[27], const int ){
    }

    __device__ MortonIndex getMortonIndex(const unsigned char* /*cell*/) const{
        return 0;
    }

    __device__ int3 getCoordinate(const unsigned char* /*cell*/) const{
        int3 coord;
        coord.x = coord.y = coord.z = 0;
        return coord;
    }

    __host__ static FCudaEmptyKernel* InitKernelKernel(void*){
        return nullptr;
    }

    __host__ static void ReleaseKernel(FCudaEmptyKernel* /*todealloc*/){
        // nothing to do
    }
};

#endif // FCUDAEMPTYKERNEL_HPP

