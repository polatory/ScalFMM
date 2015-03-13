// @SCALFMM_PRIVATE
#ifndef FCUDATESTKERNEL_HPP
#define FCUDATESTKERNEL_HPP

#include "FCudaGlobal.hpp"

// We need to describe this cell
#include "../../Components/FTestCell.hpp"

class FTestCellCudaDescriptor {
    FTestCell* ptr;
public:
    __device__ FTestCellCudaDescriptor(unsigned char* inPtr)
        : ptr(reinterpret_cast<FTestCell*>(inPtr)){
    }

    __device__ long long int& dataUp(){
        return ptr->dataUp;
    }

    __device__ long long int& dataDown(){
        return ptr->dataDown;
    }

    __device__ MortonIndex getMortonIndex() const{
        return ptr->mortonIndex;
    }

    __device__ int3 getCoordinate() const{
        const int* coordinate = (const int*)&ptr->coordinate;
        int3 coord;
        coord.x = coordinate[0];
        coord.y = coordinate[1];
        coord.z = coordinate[2];
        return coord;
    }
};

class FTestCellCudaConstDescriptor {
    const FTestCell* ptr;
public:
    __device__ FTestCellCudaConstDescriptor(const unsigned char* inPtr)
        : ptr(reinterpret_cast<const FTestCell*>(inPtr)){
    }

    __device__ const long long int& dataUp()const {
        return ptr->dataUp;
    }

    __device__ const long long int& dataDown()const {
        return ptr->dataDown;
    }

    __device__ MortonIndex getMortonIndex() const{
        return ptr->mortonIndex;
    }

    __device__ int3 getCoordinate() const{
        const int* coordinate = (const int*)&ptr->coordinate;
        int3 coord;
        coord.x = coordinate[0];
        coord.y = coordinate[1];
        coord.z = coordinate[2];
        return coord;
    }
};

template< class ContainerClass >
class FTestCudaKernels {
public:
    /** Before upward */
    __device__ void P2M(unsigned char* const pole, const ContainerClass* const particles) {
        // the pole represents all particles under
        if(threadIdx.x == 0){
            FTestCellCudaDescriptor cell(pole);
            cell.dataUp() += particles->getNbParticles();
        }
    }

    /** During upward */
    __device__ void M2M(unsigned char* const  pole, const unsigned char*const*const  child, const int /*level*/) {
        if(threadIdx.x == 0) {
            FTestCellCudaDescriptor cell(pole);
            // A parent represents the sum of the child
            for(int idx = 0 ; idx < 8 ; ++idx){
                if(child[idx]){
                    FTestCellCudaConstDescriptor childCell(child[idx]);
                    cell.dataUp() += childCell.dataUp();
                }
            }
        }
    }

    /** Before Downward */
    __device__ void M2L(unsigned char* const  local, const unsigned char* distantNeighbors[343], const int /*size*/, const int /*level*/) {
        if(threadIdx.x == 0) {
            FTestCellCudaDescriptor cell(local);
            // The pole is impacted by what represent other poles
            for(int idx = 0 ; idx < 343 ; ++idx){
                if(distantNeighbors[idx]){
                    FTestCellCudaConstDescriptor interCell(distantNeighbors[idx]);
                    cell.dataDown() += interCell.dataUp();
                }
            }
        }
    }

    /** During Downward */
    __device__ void L2L(const unsigned char*const  local, unsigned char**const  child, const int /*level*/) {
        if(threadIdx.x == 0) {
            FTestCellCudaConstDescriptor cell(local);
            // Each child is impacted by the father
            for(int idx = 0 ; idx < 8 ; ++idx){
                if(child[idx]){
                    FTestCellCudaDescriptor cellChild(child[idx]);
                    cellChild.dataDown() += cell.dataDown();
                }
            }
        }
    }

    /** After Downward */
    __device__ void L2P(const unsigned char* const  local, ContainerClass*const particles){
        if(threadIdx.x == 0) {
            FTestCellCudaConstDescriptor cell(local);
            // The particles is impacted by the parent cell
            long long int*const particlesAttributes = particles->template getAttribute<0>();
            for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
                particlesAttributes[idxPart] += cell.dataDown();
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
        FTestCellCudaConstDescriptor cellAccess(cell);
        return cellAccess.getMortonIndex();
    }

    __device__ int3 getCoordinate(const unsigned char* cell) const{
        FTestCellCudaConstDescriptor cellAccess(cell);
        return cellAccess.getCoordinate();
    }

    __host__ static FTestCudaKernels* InitKernelKernel(void*){
        return nullptr;
    }

    __host__ static void ReleaseKernel(FTestCudaKernels* /*todealloc*/){
        // nothing to do
    }
};


#endif // FCUDATESTKERNEL_HPP

