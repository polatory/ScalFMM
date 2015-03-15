// @SCALFMM_PRIVATE

#include "FCudaDeviceWrapper.hpp"
#include "FCudaTreeCoordinate.hpp"


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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__bottomPassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                         unsigned char* containersPtr, std::size_t containersSize,
                                         CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr, leafCellsSize);
    ParticleContainerGroupClass containers(containersPtr, containersSize);

    const MortonIndex blockStartIdx = leafCells.getStartingIndex();
    const MortonIndex blockEndIdx = leafCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        unsigned char* cell = leafCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(kernel->getMortonIndex(cell) == mindex);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
            FCudaAssertLF(particles.isAttachedToSomething());
            kernel->P2M(cell, &particles);
        }
    }
}

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__bottomPassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                 unsigned char* containersPtr, std::size_t containersSize,
                                 CudaKernelClass* kernel, cudaStream_t 	currentStream){
    FCuda__bottomPassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__upwardPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                  unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
                                  CudaKernelClass* kernel, int nbSubCellGroups, int idxLevel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx].reset(subCellGroupsPtr[idx], subCellGroupsSize[idx]);
    }

    FCudaAssertLF(nbSubCellGroups != 0);
    const MortonIndex blockStartIdx = FCudaMax(currentCells.getStartingIndex(),
                                          subCellGroups[0].getStartingIndex()>>3);
    const MortonIndex blockEndIdx   = FCudaMin(currentCells.getEndingIndex(),
                                          ((subCellGroups[nbSubCellGroups-1].getEndingIndex()-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        unsigned char* cell = currentCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(kernel->getMortonIndex(cell) == mindex);
            unsigned char* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( subCellGroups[idxSubCellGroup].getEndingIndex() <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                child[idxChild] = subCellGroups[idxSubCellGroup].getCell((mindex<<3)+idxChild);
                FCudaAssertLF(child[idxChild] == nullptr || kernel->getMortonIndex(child[idxChild]) == ((mindex<<3)+idxChild));
            }

            kernel->M2M(cell, child, idxLevel);
        }
    }
}

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__upwardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                 unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
                                 int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream){

    FCuda__upwardPassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInoutPassPerformMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                             unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                             CudaKernelClass* kernel, int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                             int nbOutsideInteractions){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        unsigned char* interCell = cellsOther.getCell(outsideInteractions[outInterIdx].outIndex);
        if(interCell){
            FCudaAssertLF(kernel->getMortonIndex(interCell) == outsideInteractions[outInterIdx].outIndex);
            unsigned char* cell = currentCells.getCell(outsideInteractions[outInterIdx].insideIndex);
            FCudaAssertLF(cell);
            FCudaAssertLF(kernel->getMortonIndex(cell) == outsideInteractions[outInterIdx].insideIndex);

            const unsigned char* interactions[343];
            memset(interactions, 0, 343*sizeof(unsigned char*));
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            kernel->M2L( cell , interactions, counter, idxLevel);
        }
    }
}

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInoutPassCallbackMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                           unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                           int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                           int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__transferInoutPassPerformMpi
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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


template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                       CudaKernelClass* kernel, int idxLevel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);

    const MortonIndex blockStartIdx = currentCells.getStartingIndex();
    const MortonIndex blockEndIdx = currentCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        unsigned char* cell = currentCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(kernel->getMortonIndex(cell) == mindex);
            MortonIndex interactionsIndexes[189];
            int interactionsPosition[189];
            const int3 coord = (kernel->getCoordinate(cell));
            int counter = FCudaTreeCoordinate::GetInteractionNeighbors(coord, idxLevel,interactionsIndexes,interactionsPosition);

            const unsigned char* interactions[343];
            memset(interactions, 0, 343*sizeof(unsigned char*));
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    unsigned char* interCell = currentCells.getCell(interactionsIndexes[idxInter]);
                    if(interCell){
                        FCudaAssertLF(kernel->getMortonIndex(interCell) == interactionsIndexes[idxInter]);
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                    int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream){

    FCuda__transferInPassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
                                 kernel, idxLevel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}


template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__transferInoutPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                         unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                         CudaKernelClass* kernel, int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                         int nbOutsideInteractions){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        unsigned char* interCell = cellsOther.getCell(outsideInteractions[outInterIdx].outIndex);
        if(interCell){
            FCudaAssertLF(kernel->getMortonIndex(interCell) == outsideInteractions[outInterIdx].outIndex);
            unsigned char* cell = currentCells.getCell(outsideInteractions[outInterIdx].insideIndex);
            FCudaAssertLF(cell);
            FCudaAssertLF(kernel->getMortonIndex(cell) == outsideInteractions[outInterIdx].insideIndex);

            const unsigned char* interactions[343];
            memset(interactions, 0, 343*sizeof(unsigned char*));
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            kernel->M2L( cell , interactions, counter, idxLevel);

            interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
            interactions[FMGetOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
            kernel->M2L( interCell , interactions, counter, idxLevel);
        }
    }
}


template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInoutPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                        unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                        int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                        int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );


    FCuda__transferInoutPassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__downardPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                   unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
                                   CudaKernelClass* kernel, int nbSubCellGroups, int idxLevel){
    FCudaAssertLF(nbSubCellGroups != 0);
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize);
    CellContainerClass subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx].reset(subCellGroupsPtr[idx], subCellGroupsSize[idx]);
    }

    const MortonIndex blockStartIdx = FCudaMax(currentCells.getStartingIndex(),
                                          subCellGroups[0].getStartingIndex()>>3);
    const MortonIndex blockEndIdx   = FCudaMin(currentCells.getEndingIndex(),
                                          ((subCellGroups[nbSubCellGroups-1].getEndingIndex()-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        unsigned char* cell = currentCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(kernel->getMortonIndex(cell) == mindex);
            unsigned char* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( subCellGroups[idxSubCellGroup].getEndingIndex() <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                child[idxChild] = subCellGroups[idxSubCellGroup].getCell((mindex<<3)+idxChild);
                FCudaAssertLF(child[idxChild] == nullptr || kernel->getMortonIndex(child[idxChild]) == ((mindex<<3)+idxChild));
            }

            kernel->L2L(cell, child, idxLevel);
        }
    }
}

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__downardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                  unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
                                  int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream){

    FCuda__downardPassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
            (currentCellsPtr, currentCellsSize, subCellGroupsPtr, subCellGroupsSize, kernel, nbSubCellGroups, idxLevel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallbackMpi(unsigned char* containersPtr, std::size_t containersSize,
                                         unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                         const OutOfBlockInteraction* outsideInteractions,
                                         int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerformMpi
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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


template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInPassCallback(unsigned char* containersPtr, std::size_t containersSize,
                                   const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream){
    FCuda__directInPassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize,
                               kernel, treeHeight);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
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

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
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
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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


template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__mergePassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                 unsigned char* containersPtr, std::size_t containersSize,
                                 CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr,leafCellsSize);
    ParticleContainerGroupClass containers(containersPtr,containersSize);

    const MortonIndex blockStartIdx = leafCells.getStartingIndex();
    const MortonIndex blockEndIdx = leafCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        unsigned char* cell = leafCells.getCell(mindex);
        if(cell){
            FCudaAssertLF(kernel->getMortonIndex(cell) == mindex);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
            FCudaAssertLF(particles.isAttachedToSomething());
            kernel->L2P(cell, &particles);
        }
    }
}

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__mergePassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
                                unsigned char* containersPtr, std::size_t containersSize,
                                CudaKernelClass* kernel, cudaStream_t currentStream){
    FCuda__mergePassPerform
            < CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
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

template void FCuda__bottomPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__upwardPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
     int nbSubCellGroups, int idxLevel, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallbackMpi<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    int idxLevel, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__downardPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
    int nbSubCellGroups, int idxLevel, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallbackMpi<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    const int treeHeight, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
     const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);

template void FCuda__mergePassCallback<FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int>>* kernel, cudaStream_t currentStream);


template FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int> >* FCuda__BuildCudaKernel< FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int> > >(void* kernel);
template void FCuda__ReleaseCudaKernel< FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int> > >(FCudaEmptyKernel< FCudaGroupAttachedLeaf<0, int> >* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "FCudaTestKernels.hpp"

template void FCuda__bottomPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__upwardPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
     int nbSubCellGroups, int idxLevel, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallbackMpi<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
     int idxLevel, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
     int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__downardPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
     int nbSubCellGroups, int idxLevel, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallbackMpi<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    const int treeHeight, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template void FCuda__mergePassCallback<FCudaGroupOfCells<sizeof(FTestCell)>, FCudaGroupOfParticles<2, long long int>, FCudaGroupAttachedLeaf<2, long long int>, FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int>>* kernel, cudaStream_t currentStream);

template FTestCudaKernels<FCudaGroupAttachedLeaf<2, long long int> >* FCuda__BuildCudaKernel<FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int> > >(void* kernel);
template void FCuda__ReleaseCudaKernel< FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int> > >(FTestCudaKernels< FCudaGroupAttachedLeaf<2, long long int> >* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
