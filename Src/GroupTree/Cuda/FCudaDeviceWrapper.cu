// @SCALFMM_PRIVATE

#include "FCudaDeviceWrapper.hpp"
#include "FCudaTreeCoordinate.hpp"
#include "FCudaStructParams.hpp"


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

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__bottomPassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
                                         unsigned char* containersPtr, std::size_t containersSize,
                                         CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr, leafCellsSize, leafCellsUpPtr, nullptr);
    ParticleContainerGroupClass containers(containersPtr, containersSize, nullptr);

    for(int leafIdx = 0 ; leafIdx < leafCells.getNumberOfCellsInBlock() ; ++leafIdx){
        typename CellContainerClass::CompleteCellClass cell = leafCells.getUpCell(leafIdx);
        ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(leafIdx);
        FCudaAssertLF(leafCells.getCellMortonIndex(leafIdx) == containers.getLeafMortonIndex(leafIdx));
        kernel->P2M(cell, &particles);
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__bottomPassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
                                        unsigned char* containersPtr, std::size_t containersSize,
                                        CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    FCuda__bottomPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>
                            (leafCellsPtr, leafCellsSize,leafCellsUpPtr,
                             containersPtr, containersSize,
                             kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Upward Pass
/////////////////////////////////////////////////////////////////////////////////////

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__upwardPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
                                         FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                         FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
                                         int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize,currentCellsUpPtr,nullptr);
    CellContainerClass subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx].reset(subCellGroupsPtr.values[idx], subCellGroupsSize.values[idx], subCellGroupsUpPtr.values[idx], nullptr);
    }

    FCudaAssertLF(nbSubCellGroups != 0);
    int idxSubCellGroup = 0;
    int idxChildCell = subCellGroups[0].getFistChildIdx(currentCells.getCellMortonIndex(0));
    FCudaAssertLF(idxChildCell != -1);

    for(int cellIdx = 0 ; cellIdx < currentCells.getNumberOfCellsInBlock() ; ++cellIdx){
        typename CellContainerClass::CompleteCellClass cell = currentCells.getUpCell(cellIdx);
        FCudaAssertLF(cell.symb->mortonIndex == currentCells.getCellMortonIndex(cellIdx));
        typename CellContainerClass::CompleteCellClass child[8];

        FCudaAssertLF(idxSubCellGroup != nbSubCellGroups);

        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            child[idxChild].symb = nullptr;
        }

        while(idxSubCellGroup != nbSubCellGroups
              && (subCellGroups[idxSubCellGroup].getCellMortonIndex(idxChildCell)>>3) == cell.symb->mortonIndex){
            const int idxChild = ((subCellGroups[idxSubCellGroup].getCellMortonIndex(idxChildCell)) & 7);
            FCudaAssertLF(child[idxChild].symb == nullptr);
            child[idxChild] = subCellGroups[idxSubCellGroup].getUpCell(idxChildCell);

            idxChildCell += 1;
            if(idxChildCell == subCellGroups[idxSubCellGroup].getNumberOfCellsInBlock()){
                idxChildCell = 0;
                idxSubCellGroup += 1;
            }
        }

        kernel->M2M(cell, child, idxLevel);
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__upwardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
                                        FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                        FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
                                        int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){

    FCuda__upwardPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>
                            (currentCellsPtr, currentCellsSize,currentCellsUpPtr,
                             subCellGroupsPtr, subCellGroupsSize,subCellGroupsUpPtr,
                             nbSubCellGroups, idxLevel, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}



/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass Mpi
/////////////////////////////////////////////////////////////////////////////////////
#ifdef SCALFMM_USE_MPI
template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInoutPassPerformMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
                                                  unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
                                                  int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                                  int nbOutsideInteractions, CudaKernelClass* kernel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize, nullptr, currentCellsDownPtr);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize, externalCellsUpPtr, nullptr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        const int cellPos = cellsOther.getCellIndex(outsideInteractions[outInterIdx].outIndex);
        if(cellPos != -1){
            typename CellContainerClass::CompleteCellClass interCell = cellsOther.getUpCell(cellPos);
            FCudaAssertLF(interCell.symb->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            typename CellContainerClass::CompleteCellClass cell = currentCells.getDownCell(outsideInteractions[outInterIdx].insideIdxInBlock);
            FCudaAssertLF(cell.symb->mortonIndex == outsideInteractions[outInterIdx].insideIndex);

            typename CellContainerClass::CompleteCellClass interactions[343];
            memset(interactions, 0, 343*sizeof(interactions[0]));
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            kernel->M2L( cell , interactions, counter, idxLevel);
        }
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInoutPassCallbackMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
                                                  unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
                                                  int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                                  int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__transferInoutPassPerformMpi
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(currentCellsPtr, currentCellsSize, currentCellsDownPtr,
                                       externalCellsPtr, externalCellsSize, externalCellsUpPtr,
                                       idxLevel, cuOutsideInteractions, nbOutsideInteractions, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}
#endif

/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass
/////////////////////////////////////////////////////////////////////////////////////


template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                              unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
                                              int idxLevel, CudaKernelClass* kernel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize, currentCellsUpPtr, currentCellsDownPtr);

    const MortonIndex blockStartIdx = currentCells.getStartingIndex();
    const MortonIndex blockEndIdx = currentCells.getEndingIndex();

    for(int cellIdx = 0 ; cellIdx < currentCells.getNumberOfCellsInBlock() ; ++cellIdx){
        typename CellContainerClass::CompleteCellClass cell = currentCells.getDownCell(cellIdx);
        FCudaAssertLF(cell.symb->mortonIndex == currentCells.getCellMortonIndex(cellIdx));
        MortonIndex interactionsIndexes[189];
        int interactionsPosition[189];
        const int3 coord = (FCudaTreeCoordinate::ConvertCoordinate(cell.symb->coordinates));
        int counter = FCudaTreeCoordinate::GetInteractionNeighbors(coord, idxLevel,interactionsIndexes,interactionsPosition);

        typename CellContainerClass::CompleteCellClass interactions[343];
        memset(interactions, 0, 343*sizeof(interactions[0]));
        int counterExistingCell = 0;

        for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
            if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                const int cellPos = currentCells.getCellIndex(interactionsIndexes[idxInter]);
                if(cellPos != -1){
                    typename CellContainerClass::CompleteCellClass interCell = currentCells.getUpCell(cellPos);
                    FCudaAssertLF(interCell.symb->mortonIndex == interactionsIndexes[idxInter]);
                    FCudaAssertLF(interactions[interactionsPosition[idxInter]].symb == nullptr);
                    interactions[interactionsPosition[idxInter]] = interCell;
                    counterExistingCell += 1;
                }
            }
        }

        kernel->M2L( cell , interactions, counterExistingCell, idxLevel);
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                            unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
                                            int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){

    FCuda__transferInPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
                                                                currentCellsUpPtr, currentCellsDownPtr,
                                                                idxLevel, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}


template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__transferInoutPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                                unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
                                                unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                                unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
                                                int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                                int nbOutsideInteractions, CudaKernelClass* kernel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize, currentCellsUpPtr, currentCellsDownPtr);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize, externalCellsUpPtr, externalCellsDownPtr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        const int cellPos = cellsOther.getCellIndex(outsideInteractions[outInterIdx].outIndex);
        if(cellPos != -1){
            typename CellContainerClass::CompleteCellClass interCell = cellsOther.getCompleteCell(outsideInteractions[outInterIdx].outIndex);
            FCudaAssertLF(interCell.symb->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            typename CellContainerClass::CompleteCellClass cell = currentCells.getCompleteCell(outsideInteractions[outInterIdx].insideIdxInBlock);
            FCudaAssertLF(cell.symb);
            FCudaAssertLF(cell.symb->mortonIndex == outsideInteractions[outInterIdx].insideIndex);

            typename CellContainerClass::CompleteCellClass interactions[343];
            memset(interactions, 0, 343*sizeof(interactions[0]));
            interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
            const int counter = 1;
            kernel->M2L( cell , interactions, counter, idxLevel);

            interactions[outsideInteractions[outInterIdx].outPosition].symb = nullptr;
            interactions[FMGetOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
            kernel->M2L( interCell , interactions, counter, idxLevel);
        }
    }
}


template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInoutPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                               unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
                                               unsigned char* externalCellsPtr, std::size_t externalCellsSize,
                                               unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
                                               int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                               int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );


    FCuda__transferInoutPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
                                                                currentCellsUpPtr, currentCellsDownPtr,
                                                                externalCellsPtr, externalCellsSize,
                                                                externalCellsUpPtr, externalCellsDownPtr,
                                                                idxLevel, cuOutsideInteractions,
                                                                nbOutsideInteractions, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Downard Pass
/////////////////////////////////////////////////////////////////////////////////////

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__downardPassPerform(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
                                          FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
                                          FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
                                          int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel){
    FCudaAssertLF(nbSubCellGroups != 0);
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize, nullptr, currentCellsDownPtr);
    CellContainerClass subCellGroups[9];
    for(int idx = 0 ; idx < nbSubCellGroups ; ++idx){
        subCellGroups[idx].reset(subCellGroupsPtr.values[idx], subCellGroupsSize.values[idx], nullptr, subCellGroupsDownPtr.values[idx]);
    }

    int idxSubCellGroup = 0;
    int idxChildCell = subCellGroups[0].getFistChildIdx(currentCells.getCellMortonIndex(0));
    FCudaAssertLF(idxChildCell != -1);

    for(int cellIdx = 0 ; cellIdx < currentCells.getNumberOfCellsInBlock() ; ++cellIdx){
        typename CellContainerClass::CompleteCellClass cell = currentCells.getDownCell(cellIdx);
        FCudaAssertLF(cell.symb->mortonIndex == currentCells.getCellMortonIndex(cellIdx));
        typename CellContainerClass::CompleteCellClass child[8];

        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            child[idxChild].symb = nullptr;
        }

        while(idxSubCellGroup != nbSubCellGroups
              && (subCellGroups[idxSubCellGroup].getCellMortonIndex(idxChildCell)>>3) == cell.symb->mortonIndex){
            const int idxChild = ((subCellGroups[idxSubCellGroup].getCellMortonIndex(idxChildCell)) & 7);
            FCudaAssertLF(child[idxChild].symb == nullptr);
            child[idxChild] = subCellGroups[idxSubCellGroup].getDownCell(idxChildCell);

            idxChildCell += 1;
            if(idxChildCell == subCellGroups[idxSubCellGroup].getNumberOfCellsInBlock()){
                idxChildCell = 0;
                idxSubCellGroup += 1;
            }
        }

        kernel->L2L(cell, child, idxLevel);
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__downardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
                                         FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
                                         FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
                                         int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){

    FCuda__downardPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>
            (currentCellsPtr, currentCellsSize, currentCellsDownPtr, subCellGroupsPtr, subCellGroupsSize, subCellGroupsDownPtr,
             nbSubCellGroups, idxLevel, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////
#ifdef SCALFMM_USE_MPI
template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInoutPassPerformMpi(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                                 unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                                 const OutOfBlockInteraction* outsideInteractions,
                                                 int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel){
    ParticleContainerGroupClass containers(containersPtr, containersSize, containersDownPtr);
    ParticleContainerGroupClass containersOther(externalContainersPtr, externalContainersSize, nullptr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        const int leafPos = containersOther.getLeafIndex(outsideInteractions[outInterIdx].outIndex);
        if(leafPos != -1){
            ParticleGroupClass interParticles = containersOther.template getLeaf<ParticleGroupClass>(leafPos);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(outsideInteractions[outInterIdx].insideIdxInBlock);
            ParticleGroupClass* interactions[27];
            memset(interactions, 0, 27*sizeof(ParticleGroupClass*));
            interactions[outsideInteractions[outInterIdx].outPosition] = &interParticles;
            const int counter = 1;
            kernel->P2PRemote( FCudaTreeCoordinate::GetPositionFromMorton(outsideInteractions[outInterIdx].insideIndex, treeHeight-1), &particles, &particles , interactions, counter);
        }
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallbackMpi(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                                unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                                const OutOfBlockInteraction* outsideInteractions,
                                                int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerformMpi
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(containersPtr, containersSize, containersDownPtr,
                                  externalContainersPtr, externalContainersSize,
                                  cuOutsideInteractions, nbOutsideInteractions, treeHeight, kernel);

    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}
#endif
/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass
/////////////////////////////////////////////////////////////////////////////////////


template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInPassPerform(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                           const int treeHeight, CudaKernelClass* kernel){
    ParticleContainerGroupClass containers(containersPtr, containersSize, containersDownPtr);

    const MortonIndex blockStartIdx = containers.getStartingIndex();
    const MortonIndex blockEndIdx = containers.getEndingIndex();

    for(int leafIdx = 0 ; leafIdx < containers.getNumberOfLeavesInBlock() ; ++leafIdx){
        ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(leafIdx);
        const MortonIndex mindex = containers.getLeafMortonIndex(leafIdx);
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
                const int leafPos = containers.getLeafIndex(interactionsIndexes[idxInter]);
                if(leafPos != -1){
                    interactionsObjects[counterExistingCell] = containers.template getLeaf<ParticleGroupClass>(leafPos);
                    FCudaAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                    interactions[interactionsPosition[idxInter]] = &interactionsObjects[counterExistingCell];
                    counterExistingCell += 1;
                }
            }
        }

        kernel->P2P( coord, &particles, &particles , interactions, counterExistingCell);
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInPassCallback(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                          const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    FCuda__directInPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(containersPtr, containersSize, containersDownPtr,
                               treeHeight, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInoutPassPerform(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                              unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
                                              const OutOfBlockInteraction* outsideInteractions,
                                              int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel){
    ParticleContainerGroupClass containers(containersPtr, containersSize, containersDownPtr);
    ParticleContainerGroupClass containersOther(externalContainersPtr, externalContainersSize, externalContainersPtr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){
        const int leafPos = containersOther.getLeafIndex(outsideInteractions[outInterIdx].outIndex);
        if(leafPos != -1){
            ParticleGroupClass interParticles = containersOther.template getLeaf<ParticleGroupClass>(leafPos);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(outsideInteractions[outInterIdx].insideIdxInBlock);

            FCudaAssertLF(containersOther.getLeafMortonIndex(leafPos) == outsideInteractions[outInterIdx].outIndex);
            FCudaAssertLF(containers.getLeafMortonIndex(outsideInteractions[outInterIdx].insideIdxInBlock) == outsideInteractions[outInterIdx].insideIndex);

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

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallback(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                             unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
                                             const OutOfBlockInteraction* outsideInteractions,
                                             int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(containersPtr, containersSize,containersDownPtr,
                                  externalContainersPtr, externalContainersSize,externalContainersDownPtr,
                                  cuOutsideInteractions, nbOutsideInteractions, treeHeight, kernel);

    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}


/////////////////////////////////////////////////////////////////////////////////////
/// Merge Pass
/////////////////////////////////////////////////////////////////////////////////////


template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__mergePassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
                                        unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                        CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr,leafCellsSize, nullptr, leafCellsDownPtr);
    ParticleContainerGroupClass containers(containersPtr,containersSize, containersDownPtr);

    for(int cellIdx = 0 ; cellIdx < leafCells.getNumberOfCellsInBlock() ; ++cellIdx){
        typename CellContainerClass::CompleteCellClass cell = leafCells.getDownCell(cellIdx);
        FCudaAssertLF(cell.symb->mortonIndex == leafCells.getCellMortonIndex(cellIdx));
        ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(cellIdx);
        FCudaAssertLF(leafCells.getCellMortonIndex(cellIdx) == containers.getLeafMortonIndex(cellIdx));
        kernel->L2P(cell, &particles);
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__mergePassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
                                       unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                       CudaKernelClass* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize){
    FCuda__mergePassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<inGridSize, inBlocksSize, 0, currentStream>>>(leafCellsPtr, leafCellsSize,leafCellsDownPtr,
                            containersPtr, containersSize,containersDownPtr,
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

template <class CudaKernelClass>
dim3 FCuda__GetGridSize(CudaKernelClass* /*kernel*/, int intervalSize){
    return CudaKernelClass::GetGridSize(intervalSize);
}

template <class CudaKernelClass>
dim3 FCuda__GetBlockSize(CudaKernelClass* /*kernel*/){
    return CudaKernelClass::GetBlocksSize();
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "FCudaGroupOfCells.hpp"
#include "FCudaGroupAttachedLeaf.hpp"
#include "FCudaGroupOfParticles.hpp"

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "FCudaEmptyKernel.hpp"
#include "FCudaEmptyCellSymb.hpp"

template void FCuda__bottomPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
    unsigned char* containersPtr, std::size_t containersSize,
    FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__upwardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
int nbSubCellGroups, int idxLevel, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__transferInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__transferInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__transferInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__downardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__directInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__directInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__directInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__mergePassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                       FCudaGroupOfParticles<int,0,0,int>, FCudaGroupAttachedLeaf<int,0,0,int>, FCudaEmptyKernel<int> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FCudaEmptyKernel<int>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template FCudaEmptyKernel<int>* FCuda__BuildCudaKernel< FCudaEmptyKernel<int> >(void* kernel);
template void FCuda__ReleaseCudaKernel< FCudaEmptyKernel<int> >(FCudaEmptyKernel<int>* cukernel);
template dim3 FCuda__GetGridSize< FCudaEmptyKernel<int> >(FCudaEmptyKernel<int>* kernel, int intervalSize);
template dim3 FCuda__GetBlockSize< FCudaEmptyKernel<int> >(FCudaEmptyKernel<int>* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "../TestKernel/FCudaTestKernels.hpp"
#include "../TestKernel/FTestCellPOD.hpp"

template void FCuda__bottomPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__upwardPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
    int nbSubCellGroups, int idxLevel, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__transferInoutPassCallbackMpi<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__transferInPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__transferInoutPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__downardPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__directInoutPassCallbackMpi<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__directInPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__directInoutPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__mergePassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<float,0, 1, long long int>, FCudaGroupAttachedLeaf<float,0, 1, long long int>, FTestCudaKernels<float> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FTestCudaKernels<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template FTestCudaKernels<float>* FCuda__BuildCudaKernel<FTestCudaKernels<float>>(void* kernel);
template void FCuda__ReleaseCudaKernel<FTestCudaKernels<float>>(FTestCudaKernels<float>* cukernel);

template dim3 FCuda__GetGridSize< FTestCudaKernels<float> >(FTestCudaKernels<float>* kernel, int intervalSize);
template dim3 FCuda__GetBlockSize< FTestCudaKernels<float> >(FTestCudaKernels<float>* cukernel);




template void FCuda__bottomPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__upwardPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
    int nbSubCellGroups, int idxLevel, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__transferInoutPassCallbackMpi<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__transferInPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__transferInoutPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__downardPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__directInoutPassCallbackMpi<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__directInPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__directInoutPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__mergePassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<double,0, 1, long long int>, FCudaGroupAttachedLeaf<double,0, 1, long long int>, FTestCudaKernels<double> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FTestCudaKernels<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template FTestCudaKernels<double>* FCuda__BuildCudaKernel<FTestCudaKernels<double>>(void* kernel);
template void FCuda__ReleaseCudaKernel<FTestCudaKernels<double>>(FTestCudaKernels<double>* cukernel);

template dim3 FCuda__GetGridSize< FTestCudaKernels<double> >(FTestCudaKernels<double>* kernel, int intervalSize);
template dim3 FCuda__GetBlockSize< FTestCudaKernels<double> >(FTestCudaKernels<double>* cukernel);


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "../P2P/FCudaP2P.hpp"

template void FCuda__bottomPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
unsigned char* containersPtr, std::size_t containersSize,
    FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__upwardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
    int nbSubCellGroups, int idxLevel, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__transferInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__transferInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__transferInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__downardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__directInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__directInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__directInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__mergePassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<float,4, 4, float>, FCudaGroupAttachedLeaf<float,4, 4, float>, FCudaP2P<float> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FCudaP2P<float>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template FCudaP2P<float>* FCuda__BuildCudaKernel<FCudaP2P<float>>(void* kernel);
template void FCuda__ReleaseCudaKernel<FCudaP2P<float>>(FCudaP2P<float>* cukernel);

template dim3 FCuda__GetGridSize< FCudaP2P<float> >(FCudaP2P<float>* kernel, int intervalSize);
template dim3 FCuda__GetBlockSize< FCudaP2P<float> >(FCudaP2P<float>* cukernel);




template void FCuda__bottomPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
unsigned char* containersPtr, std::size_t containersSize,
    FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__upwardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
    int nbSubCellGroups, int idxLevel, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__transferInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__transferInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__transferInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__downardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#ifdef SCALFMM_USE_MPI
template void FCuda__directInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);
#endif
template void FCuda__directInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__directInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template void FCuda__mergePassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<double,4, 4, double>, FCudaGroupAttachedLeaf<double,4, 4, double>, FCudaP2P<double> >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FCudaP2P<double>* kernel, cudaStream_t currentStream,
                                        const dim3 inGridSize, const dim3 inBlocksSize);

template FCudaP2P<double>* FCuda__BuildCudaKernel<FCudaP2P<double>>(void* kernel);
template void FCuda__ReleaseCudaKernel<FCudaP2P<double>>(FCudaP2P<double>* cukernel);

template dim3 FCuda__GetGridSize< FCudaP2P<double> >(FCudaP2P<double>* kernel, int intervalSize);
template dim3 FCuda__GetBlockSize< FCudaP2P<double> >(FCudaP2P<double>* cukernel);
