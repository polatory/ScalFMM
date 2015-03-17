// @SCALFMM_PRIVATE

#include "FCudaDeviceWrapper.hpp"
#include "FCudaTreeCoordinate.hpp"
#include "FCudaStructParams.hpp"

static const int nbCudaThreads = 32;
static const int nbCudaBlocks = 1;

static void FCudaCheckCore(cudaError_t code, const char *file, int line) {
   if (code != cudaSuccess) {
      fprintf(stderr,"Cuda Error %d : %s %s %d\n", code, cudaGetErrorString(code), file, line);
      exit(code);
   }
}
#define FCudaCheck( test ) { FCudaCheckCore((test), __FILE__, __LINE__); }
#define FCudaCheckAfterCall() { FCudaCheckCore((cudaGetLastError()), __FILE__, __LINE__); }
#define FCudaAssertLF(ARGS) if(!(ARGS)){\
                                printf("Error line %d\n", __LINE__);\
                            }

#define FMGetOppositeNeighIndex(index) (27-(index)-1)
#define FMGetOppositeInterIndex(index) (343-(index)-1)

#define FCudaMax(x,y) ((x)<(y) ? (y) : (x))
#define FCudaMin(x,y) ((x)>(y) ? (y) : (x))

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__bottomPassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                         unsigned char* containersPtr, std::size_t containersSize,
                                         CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr, leafCellsSize);
    ParticleContainerGroupClass containers(containersPtr, containersSize);

    const MortonIndex blockStartIdx = leafCells.getStartingIndex();
    const MortonIndex blockEndIdx = leafCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        CellClass* cell = leafCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(cell->mortonIndex == mindex);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
            FCudaAssertLF(particles.isAttachedToSomething());
            kernel->P2M(cell, &particles);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__bottomPassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                 unsigned char* containersPtr, std::size_t containersSize,
                                 CudaKernelClass* kernel, cudaStream_t 	currentStream){
    FCuda__bottomPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
                            (leafCellsPtr, leafCellsSize,
                             containersPtr, containersSize,
                             kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Upward Pass
/////////////////////////////////////////////////////////////////////////////////////

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__upwardPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                  FCudaParams<unsigned char*, 9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                  CudaKernelClass* kernel, int nbSubCellGroups, int idxLevel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx].reset(subCellGroupsPtr.values[idx], subCellGroupsSize.values[idx]);
    }

    FCudaAssertLF(nbSubCellGroups != 0);
    const MortonIndex blockStartIdx = FCudaMax(currentCells.getStartingIndex(),
                                          subCellGroups[0].getStartingIndex()>>3);
    const MortonIndex blockEndIdx   = FCudaMin(currentCells.getEndingIndex(),
                                          ((subCellGroups[nbSubCellGroups-1].getEndingIndex()-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        CellClass* cell = currentCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(cell->mortonIndex == mindex);
            const CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( subCellGroups[idxSubCellGroup].getEndingIndex() <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                child[idxChild] = subCellGroups[idxSubCellGroup].getCell((mindex<<3)+idxChild);
                FCudaAssertLF(child[idxChild] == nullptr || child[idxChild]->mortonIndex == ((mindex<<3)+idxChild));
            }

            kernel->M2M(cell, child, idxLevel);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__upwardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                 FCudaParams<unsigned char*, 9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                 int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream){

    FCuda__upwardPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
                            (currentCellsPtr, currentCellsSize,
                             subCellGroupsPtr, subCellGroupsSize,
                             kernel, nbSubCellGroups, idxLevel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}



/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass Mpi
/////////////////////////////////////////////////////////////////////////////////////

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInoutPassPerformMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                             unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                             CudaKernelClass* kernel, int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                             int nbOutsideInteractions){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        CellClass* interCell = cellsOther.getCell(outsideInteractions[outInterIdx].outIndex);
        if(interCell){
            FCudaAssertLF(interCell->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            CellClass* cell = currentCells.getCell(outsideInteractions[outInterIdx].insideIndex);
            FCudaAssertLF(cell);
            FCudaAssertLF(cell->mortonIndex == outsideInteractions[outInterIdx].insideIndex);

            const CellClass* interactions[343];
            memset(interactions, 0, 343*sizeof(CellClass*));
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            kernel->M2L( cell , interactions, counter, idxLevel);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInoutPassCallbackMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                           unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                           int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                           int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__transferInoutPassPerformMpi
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
                                       externalCellsPtr, externalCellsSize,
                                       kernel, idxLevel, cuOutsideInteractions,
                                       nbOutsideInteractions);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass
/////////////////////////////////////////////////////////////////////////////////////


template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                       CudaKernelClass* kernel, int idxLevel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);

    const MortonIndex blockStartIdx = currentCells.getStartingIndex();
    const MortonIndex blockEndIdx = currentCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        CellClass* cell = currentCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(cell->mortonIndex == mindex);
            MortonIndex interactionsIndexes[189];
            int interactionsPosition[189];
            const int3 coord = (FCudaTreeCoordinate::ConvertCoordinate(cell->coordinates));
            int counter = FCudaTreeCoordinate::GetInteractionNeighbors(coord, idxLevel,interactionsIndexes,interactionsPosition);

            const CellClass* interactions[343];
            memset(interactions, 0, 343*sizeof(CellClass*));
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    CellClass* interCell = currentCells.getCell(interactionsIndexes[idxInter]);
                    if(interCell){
                        FCudaAssertLF(interCell->mortonIndex == interactionsIndexes[idxInter]);
                        FCudaAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                        interactions[interactionsPosition[idxInter]] = interCell;
                        counterExistingCell += 1;
                    }
                }
            }

            kernel->M2L( cell , interactions, counterExistingCell, idxLevel);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                    int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream){

    FCuda__transferInPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
                                 kernel, idxLevel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}


template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__transferInoutPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                         unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                         CudaKernelClass* kernel, int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                         int nbOutsideInteractions){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        CellClass* interCell = cellsOther.getCell(outsideInteractions[outInterIdx].outIndex);
        if(interCell){
            FCudaAssertLF(interCell->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            CellClass* cell = currentCells.getCell(outsideInteractions[outInterIdx].insideIndex);
            FCudaAssertLF(cell);
            FCudaAssertLF(cell->mortonIndex == outsideInteractions[outInterIdx].insideIndex);

            const CellClass* interactions[343];
            memset(interactions, 0, 343*sizeof(CellClass*));
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            kernel->M2L( cell , interactions, counter, idxLevel);

            interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
            interactions[FMGetOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
            kernel->M2L( interCell , interactions, counter, idxLevel);
        }
    }
}


template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInoutPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                        unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                        int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                        int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );


    FCuda__transferInoutPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
                                    externalCellsPtr, externalCellsSize,
                                    kernel, idxLevel, cuOutsideInteractions,
                                    nbOutsideInteractions);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Downard Pass
/////////////////////////////////////////////////////////////////////////////////////

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__downardPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                   FCudaParams<unsigned char*, 9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                   CudaKernelClass* kernel, int nbSubCellGroups, int idxLevel){
    FCudaAssertLF(nbSubCellGroups != 0);
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx].reset(subCellGroupsPtr.values[idx], subCellGroupsSize.values[idx]);
    }

    const MortonIndex blockStartIdx = FCudaMax(currentCells.getStartingIndex(),
                                          subCellGroups[0].getStartingIndex()>>3);
    const MortonIndex blockEndIdx   = FCudaMin(currentCells.getEndingIndex(),
                                          ((subCellGroups[nbSubCellGroups-1].getEndingIndex()-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        CellClass* cell = currentCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(cell->mortonIndex == mindex);
            CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( subCellGroups[idxSubCellGroup].getEndingIndex() <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                child[idxChild] = subCellGroups[idxSubCellGroup].getCell((mindex<<3)+idxChild);
                FCudaAssertLF(child[idxChild] == nullptr || child[idxChild]->mortonIndex == ((mindex<<3)+idxChild));
            }

            kernel->L2L(cell, child, idxLevel);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__downardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                  FCudaParams<unsigned char*, 9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                  int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream){

    FCuda__downardPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
            (currentCellsPtr, currentCellsSize, subCellGroupsPtr, subCellGroupsSize, kernel, nbSubCellGroups, idxLevel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInoutPassPerformMpi(unsigned char* containersPtr, std::size_t containersSize,
                                          unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                          CudaKernelClass* kernel, const OutOfBlockInteraction* outsideInteractions,
                                          int nbOutsideInteractions, const int treeHeight){
    ParticleContainerGroupClass containers(containersPtr, containersSize);
    ParticleContainerGroupClass containersOther(externalContainersPtr, externalContainersSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        ParticleGroupClass interParticles = containersOther.template getLeaf<ParticleGroupClass>(outsideInteractions[outInterIdx].outIndex);
        if(interParticles.isAttachedToSomething()){
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(outsideInteractions[outInterIdx].insideIndex);
            FCudaAssertLF(particles.isAttachedToSomething());
            ParticleGroupClass* interactions[27];
            memset(interactions, 0, 27*sizeof(ParticleGroupClass*));
            interactions[outsideInteractions[outInterIdx].outPosition] = &interParticles;
            const int counter = 1;
            kernel->P2PRemote( FCudaTreeCoordinate::GetPositionFromMorton(outsideInteractions[outInterIdx].insideIndex, treeHeight-1), &particles, &particles , interactions, counter);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallbackMpi(unsigned char* containersPtr, std::size_t containersSize,
                                         unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                         const OutOfBlockInteraction* outsideInteractions,
                                         int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerformMpi
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize,
                                  externalContainersPtr, externalContainersSize,
                                  kernel, cuOutsideInteractions,
                                  nbOutsideInteractions, treeHeight);

    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}

/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass
/////////////////////////////////////////////////////////////////////////////////////


template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInPassPerform(unsigned char* containersPtr, std::size_t containersSize,
                                    CudaKernelClass* kernel, const int treeHeight){
    ParticleContainerGroupClass containers(containersPtr, containersSize);

    const MortonIndex blockStartIdx = containers.getStartingIndex();
    const MortonIndex blockEndIdx = containers.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
        if(particles.isAttachedToSomething()){
            MortonIndex interactionsIndexes[26];
            int interactionsPosition[26];
            const int3 coord = FCudaTreeCoordinate::GetPositionFromMorton(mindex, treeHeight-1);
            int counter = FCudaTreeCoordinate::GetNeighborsIndexes(coord, treeHeight,interactionsIndexes,interactionsPosition);

            ParticleGroupClass interactionsObjects[27];
            ParticleGroupClass* interactions[27];
            memset(interactions, 0, 27*sizeof(ParticleGroupClass*));
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    interactionsObjects[counterExistingCell] = containers.template getLeaf<ParticleGroupClass>(interactionsIndexes[idxInter]);
                    if(interactionsObjects[counterExistingCell].isAttachedToSomething()){
                        FCudaAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                        interactions[interactionsPosition[idxInter]] = &interactionsObjects[counterExistingCell];
                        counterExistingCell += 1;
                    }
                }
            }

            kernel->P2P( coord, &particles, &particles , interactions, counterExistingCell);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInPassCallback(unsigned char* containersPtr, std::size_t containersSize,
                                   const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream){
    FCuda__directInPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize,
                               kernel, treeHeight);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInoutPassPerform(unsigned char* containersPtr, std::size_t containersSize,
                                       unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                       CudaKernelClass* kernel, const OutOfBlockInteraction* outsideInteractions,
                                       int nbOutsideInteractions, const int treeHeight){
    ParticleContainerGroupClass containers(containersPtr, containersSize);
    ParticleContainerGroupClass containersOther(externalContainersPtr, externalContainersSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        ParticleGroupClass interParticles = containersOther.template getLeaf<ParticleGroupClass>(outsideInteractions[outInterIdx].outIndex);
        if(interParticles.isAttachedToSomething()){
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(outsideInteractions[outInterIdx].insideIndex);
            FCudaAssertLF(particles.isAttachedToSomething());
            ParticleGroupClass* interactions[27];
            memset(interactions, 0, 27*sizeof(ParticleGroupClass*));
            interactions[outsideInteractions[outInterIdx].outPosition] = &interParticles;
            const int counter = 1;
            kernel->P2PRemote( FCudaTreeCoordinate::GetPositionFromMorton(outsideInteractions[outInterIdx].insideIndex, treeHeight-1), &particles, &particles , interactions, counter);

            interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
            interactions[FMGetOppositeNeighIndex(outsideInteractions[outInterIdx].outPosition)] = &particles;
            kernel->P2PRemote( FCudaTreeCoordinate::GetPositionFromMorton(outsideInteractions[outInterIdx].outIndex, treeHeight-1), &interParticles, &interParticles , interactions, counter);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallback(unsigned char* containersPtr, std::size_t containersSize,
                                      unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                      const OutOfBlockInteraction* outsideInteractions,
                                      int nbOutsideInteractions, const int treeHeight,
                                      CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize,
                                  externalContainersPtr, externalContainersSize,
                                  kernel, cuOutsideInteractions,
                                  nbOutsideInteractions, treeHeight);

    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Merge Pass
/////////////////////////////////////////////////////////////////////////////////////


template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__mergePassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                 unsigned char* containersPtr, std::size_t containersSize,
                                 CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr,leafCellsSize);
    ParticleContainerGroupClass containers(containersPtr,containersSize);

    const MortonIndex blockStartIdx = leafCells.getStartingIndex();
    const MortonIndex blockEndIdx = leafCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        CellClass* cell = leafCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(cell->mortonIndex == mindex);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
            FCudaAssertLF(particles.isAttachedToSomething());
            kernel->L2P(cell, &particles);
        }
    }
}

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__mergePassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                unsigned char* containersPtr, std::size_t containersSize,
                                CudaKernelClass* kernel, cudaStream_t currentStream){
    FCuda__mergePassPerform
            < CellClass, CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(leafCellsPtr, leafCellsSize,
                            containersPtr, containersSize,
                            kernel);

    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}


template <class CudaKernelClass>
CudaKernelClass* FCuda__BuildCudaKernel(void* kernel){
    return CudaKernelClass::InitKernelKernel(kernel);
}

template <class CudaKernelClass>
void FCuda__ReleaseCudaKernel(CudaKernelClass* cukernel){
    CudaKernelClass::ReleaseKernel(cukernel);
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "FCudaGroupOfCells.hpp"
#include "FCudaGroupAttachedLeaf.hpp"
#include "FCudaGroupOfParticles.hpp"

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "FCudaEmptyKernel.hpp"
#include "FCudaEmptyCell.hpp"

template void FCuda__bottomPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__upwardPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
     int nbSubCellGroups, int idxLevel, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallbackMpi<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    int idxLevel, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__downardPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    int nbSubCellGroups, int idxLevel, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallbackMpi<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    const int treeHeight, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
     const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__mergePassCallback<FCudaEmptyCell, FCudaGroupOfCells<FCudaEmptyCell>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaEmptyCell,FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);


template FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int> >* FCuda__BuildCudaKernel< FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int> > >(void* kernel);
template void FCuda__ReleaseCudaKernel< FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int> > >(FCudaEmptyKernel< FCudaEmptyCell, FCudaGroupAttachedLeaf<0, int> >* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "FCudaTestKernels.hpp"
#include "../FTestCellPOD.hpp"

template void FCuda__bottomPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__upwardPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
     int nbSubCellGroups, int idxLevel, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallbackMpi<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
     int idxLevel, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__downardPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
     int nbSubCellGroups, int idxLevel, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallbackMpi<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    const int treeHeight, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__mergePassCallback<FTestCellPODCore, FCudaGroupOfCells<FTestCellPODCore>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template FTestCudaKernels<FTestCellPODCore,FCudaGroupAttachedLeaf<2, long long int> >* FCuda__BuildCudaKernel<FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int> > >(void* kernel);
template void FCuda__ReleaseCudaKernel< FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int> > >(FTestCudaKernels<FTestCellPODCore, FCudaGroupAttachedLeaf<2, long long int> >* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
