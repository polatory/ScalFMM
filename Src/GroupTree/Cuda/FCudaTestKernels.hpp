// @SCALFMM_PRIVATE
#ifndef FCUDATESTKERNEL_HPP
#define FCUDATESTKERNEL_HPP

#include "FCudaGlobal.hpp"

// We need to describe this cell
#include "../../Components/FTestCell.hpp"


template< class CellClass, class ContainerClass >
class FTestCudaKernels {
public:
    /** Before upward */
    __device__ void P2M(CellClass* pole, const ContainerClass* const particles) {
        // the pole represents all particles under
        if(threadIdx.x == 0){
            pole->dataUp += particles->getNbParticles();
        }
    }

    /** During upward */
    __device__ void M2M(CellClass*  pole, const CellClass*  child[8], const int /*level*/) {
        if(threadIdx.x == 0) {
            // A parent represents the sum of the child
            for(int idx = 0 ; idx < 8 ; ++idx){
                if(child[idx]){
                    pole->dataUp += child[idx]->dataUp;
                }
            }
        }
    }

    /** Before Downward */
    __device__ void M2L(CellClass*  local, const CellClass* distantNeighbors[343], const int /*size*/, const int /*level*/) {
        if(threadIdx.x == 0) {
            // The pole is impacted by what represent other poles
            for(int idx = 0 ; idx < 343 ; ++idx){
                if(distantNeighbors[idx]){
                    local->dataDown += distantNeighbors[idx]->dataUp;
                }
            }
        }
    }

    /** During Downward */
    __device__ void L2L(const CellClass* local, CellClass*  child[8], const int /*level*/) {
        if(threadIdx.x == 0) {
            // Each child is impacted by the father
            for(int idx = 0 ; idx < 8 ; ++idx){
                if(child[idx]){
                    child[idx]->dataDown += local->dataDown;
                }
            }
        }
    }

    /** After Downward */
    __device__ void L2P(const CellClass* local, ContainerClass*const particles){
        if(threadIdx.x == 0) {
            // The particles is impacted by the parent cell
            long long int*const particlesAttributes = particles->template getAttribute<0>();
            for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
                particlesAttributes[idxPart] += local->dataDown;
            }
        }
    }


    /** After Downward */
    __device__ void P2P(const int3& ,
                 ContainerClass* const  targets, const ContainerClass* const  sources,
                 ContainerClass* const directNeighborsParticles[27], const int ){
        if(threadIdx.x == 0) {
            // Each particles targeted is impacted by the particles sources
            long long int inc = sources->getNbParticles();
            if(targets == sources){
                inc -= 1;
            }
            for(int idx = 0 ; idx < 27 ; ++idx){
                if( directNeighborsParticles[idx] ){
                    inc += directNeighborsParticles[idx]->getNbParticles();
                }
            }

            long long int*const particlesAttributes = targets->template getAttribute<0>();
            for(int idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
                particlesAttributes[idxPart] += inc;
            }
        }
    }

    /** After Downward */
    __device__ void P2PRemote(const int3& ,
                 ContainerClass* const  targets, const ContainerClass* const  sources,
                 ContainerClass* const directNeighborsParticles[27], const int ){
        if(threadIdx.x == 0) {
            // Each particles targeted is impacted by the particles sources
            long long int inc = 0;
            for(int idx = 0 ; idx < 27 ; ++idx){
                if( directNeighborsParticles[idx] ){
                    inc += directNeighborsParticles[idx]->getNbParticles();
                }
            }

            long long int*const particlesAttributes = targets->template getAttribute<0>();
            for(int idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
                particlesAttributes[idxPart] += inc;
            }
        }
    }

    __host__ static FTestCudaKernels* InitKernelKernel(void*){
        return nullptr;
    }

    __host__ static void ReleaseKernel(FTestCudaKernels* /*todealloc*/){
        // nothing to do
    }
};


#endif // FCUDATESTKERNEL_HPP

