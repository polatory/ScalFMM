// @SCALFMM_PRIVATE
#ifndef FCUDATESTKERNEL_HPP
#define FCUDATESTKERNEL_HPP

#include "FCudaGlobal.hpp"

// We need to describe this cell
#include "../../Components/FTestCell.hpp"
#include "../../Utils/FOffetOf.hpp"

struct FTestCellDescriptor{
    static const int offset_mortonIndex =  FOffsetOf(FTestCell, mortonIndex);
    static const int offset_FTreeCoordinate_data =  FOffsetOf(FTestCell, coordinate);

    static const int offset_dataUp = FOffsetOf(FTestCell, dataUp);
    static const int offset_dataDown = FOffsetOf(FTestCell, dataDown);

    static const int size = alignof(FTestCell)*((offset_dataDown + sizeof(FTestCell::dataDown) + alignof(FTestCell) - 1)/alignof(FTestCell));
    static_assert(size == sizeof(FTestCell), "Error in attribute shift.");
};

template< class ContainerClass >
class FTestCudaKernels {
public:
    /** Before upward */
    __device__ void P2M(unsigned char* const pole, const ContainerClass* const particles) {
        // the pole represents all particles under
        if(threadIdx.x == 0){
            long long int* dataUp = reinterpret_cast<long long int*>(&pole[FTestCellDescriptor::offset_dataUp]);
            (*dataUp) += particles->getNbParticles();
        }
    }

    /** During upward */
    __device__ void M2M(unsigned char* const  pole, const unsigned char*const*const  child, const int /*level*/) {
        if(threadIdx.x == 0) {
            long long int* dataUp = reinterpret_cast<long long int*>(&pole[FTestCellDescriptor::offset_dataUp]);
            // A parent represents the sum of the child
            for(int idx = 0 ; idx < 8 ; ++idx){
                if(child[idx]){
                    const long long int* child_dataUp = reinterpret_cast<const long long int*>(&child[idx][FTestCellDescriptor::offset_dataUp]);
                    (*dataUp) += (*child_dataUp);
                }
            }
        }
    }

    /** Before Downward */
    __device__ void M2L(unsigned char* const  local, const unsigned char* distantNeighbors[343], const int /*size*/, const int /*level*/) {
        if(threadIdx.x == 0) {
            long long int* dataDown = reinterpret_cast<long long int*>(&local[FTestCellDescriptor::offset_dataDown]);
            // The pole is impacted by what represent other poles
            for(int idx = 0 ; idx < 343 ; ++idx){
                if(distantNeighbors[idx]){
                    const long long int* dataUp = reinterpret_cast<const long long int*>(&distantNeighbors[idx][FTestCellDescriptor::offset_dataUp]);
                    (*dataDown) += (*dataUp);
                }
            }
        }
    }

    /** During Downward */
    __device__ void L2L(const unsigned char*const  local, unsigned char**const  child, const int /*level*/) {
        if(threadIdx.x == 0) {
            const long long int* dataDown = reinterpret_cast<const long long int*>(&local[FTestCellDescriptor::offset_dataDown]);
            // Each child is impacted by the father
            for(int idx = 0 ; idx < 8 ; ++idx){
                if(child[idx]){
                    long long int* child_dataDown = reinterpret_cast<long long int*>(&child[idx][FTestCellDescriptor::offset_dataDown]);
                    (*child_dataDown) = (*dataDown);
                }
            }
        }
    }

    /** After Downward */
    __device__ void L2P(const unsigned char* const  local, ContainerClass*const particles){
        if(threadIdx.x == 0) {
            const long long int* dataDown = reinterpret_cast<const long long int*>(&local[FTestCellDescriptor::offset_dataDown]);
            // The particles is impacted by the parent cell
            long long int*const particlesAttributes = particles->template getAttribute<0>();
            for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
                particlesAttributes[idxPart] += (*dataDown);
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

    __device__ MortonIndex getMortonIndex(const unsigned char* cell) const{
        return 0;// TODO
    }

    __device__ int3 getCoordinate(const unsigned char* cell) const{
        int3 coord;
        coord.x = coord.y = coord.z = 0;
        return coord;// TODO
    }

    __host__ static FTestCudaKernels* InitKernelKernel(void*){
        return nullptr;
    }

    __host__ static void ReleaseKernel(FTestCudaKernels* /*todealloc*/){
        // nothing to do
    }
};


#endif // FCUDATESTKERNEL_HPP

