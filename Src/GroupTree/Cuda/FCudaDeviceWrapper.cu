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

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__bottomPassPerform(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
                                         unsigned char* containersPtr, std::size_t containersSize,
                                         CudaKernelClass* kernel){
    CellContainerClass leafCells(leafCellsPtr, leafCellsSize, leafCellsUpPtr, nullptr);
    ParticleContainerGroupClass containers(containersPtr, containersSize, nullptr);

    const MortonIndex blockStartIdx = leafCells.getStartingIndex();
    const MortonIndex blockEndIdx = leafCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        if(leafCells.exists(mindex)){
            typename CellContainerClass::CompleteCellClass cell = leafCells.getUpCell(mindex);
            FCudaAssertLF(cell.symb->mortonIndex == mindex);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
            FCudaAssertLF(particles.isAttachedToSomething());
            kernel->P2M(cell, &particles);
        }
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__bottomPassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
                                        unsigned char* containersPtr, std::size_t containersSize,
                                        CudaKernelClass* kernel, cudaStream_t 	currentStream){
    FCuda__bottomPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
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
    const MortonIndex blockStartIdx = FCudaMax(currentCells.getStartingIndex(),
                                          subCellGroups[0].getStartingIndex()>>3);
    const MortonIndex blockEndIdx   = FCudaMin(currentCells.getEndingIndex(),
                                          ((subCellGroups[nbSubCellGroups-1].getEndingIndex()-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        if(currentCells.exists(mindex)){
            typename CellContainerClass::CompleteCellClass cell = currentCells.getUpCell(mindex);
            FCudaAssertLF(cell.symb->mortonIndex == mindex);
            typename CellContainerClass::CompleteCellClass child[8];

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( subCellGroups[idxSubCellGroup].getEndingIndex() <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                if(subCellGroups[idxSubCellGroup].exists((mindex<<3)+idxChild)){
                    child[idxChild] = subCellGroups[idxSubCellGroup].getUpCell((mindex<<3)+idxChild);
                }
                else{
                    child[idxChild].symb = nullptr;
                }
            }

            kernel->M2M(cell, child, idxLevel);
        }
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__upwardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
                                        FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
                                        FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
                                        int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t 	currentStream){

    FCuda__upwardPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
                            (currentCellsPtr, currentCellsSize,currentCellsUpPtr,
                             subCellGroupsPtr, subCellGroupsSize,subCellGroupsUpPtr,
                             nbSubCellGroups, idxLevel, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}



/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass Mpi
/////////////////////////////////////////////////////////////////////////////////////

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__  void FCuda__transferInoutPassPerformMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
                                                  unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
                                                  int idxLevel, const OutOfBlockInteraction* outsideInteractions,
                                                  int nbOutsideInteractions, CudaKernelClass* kernel){
    CellContainerClass currentCells(currentCellsPtr, currentCellsSize, nullptr, currentCellsDownPtr);
    CellContainerClass cellsOther(externalCellsPtr, externalCellsSize, externalCellsUpPtr, nullptr);

    for(int outInterIdx = 0 ; outInterIdx < nbOutsideInteractions ; ++outInterIdx){        
        if(cellsOther.exists(outsideInteractions[outInterIdx].outIndex)){
            typename CellContainerClass::CompleteCellClass interCell = cellsOther.getUpCell(outsideInteractions[outInterIdx].outIndex);
            FCudaAssertLF(interCell.symb->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            typename CellContainerClass::CompleteCellClass cell = currentCells.getDownCell(outsideInteractions[outInterIdx].insideIndex);
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
                                                  int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t 	currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__transferInoutPassPerformMpi
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize, currentCellsDownPtr,
                                       externalCellsPtr, externalCellsSize, externalCellsUpPtr,
                                       idxLevel, cuOutsideInteractions, nbOutsideInteractions, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}


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

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        if(currentCells.exists(mindex)){
            typename CellContainerClass::CompleteCellClass cell = currentCells.getDownCell(mindex);
            FCudaAssertLF(cell.symb->mortonIndex == mindex);
            MortonIndex interactionsIndexes[189];
            int interactionsPosition[189];
            const int3 coord = (FCudaTreeCoordinate::ConvertCoordinate(cell.symb->coordinates));
            int counter = FCudaTreeCoordinate::GetInteractionNeighbors(coord, idxLevel,interactionsIndexes,interactionsPosition);

            typename CellContainerClass::CompleteCellClass interactions[343];
            memset(interactions, 0, 343*sizeof(interactions[0]));
            int counterExistingCell = 0;

            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                    if(currentCells.exists(interactionsIndexes[idxInter])){
                        typename CellContainerClass::CompleteCellClass interCell = currentCells.getUpCell(interactionsIndexes[idxInter]);
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
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__transferInPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
                                            unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
                                            int idxLevel, CudaKernelClass* kernel, cudaStream_t 	currentStream){

    FCuda__transferInPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
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
        if(cellsOther.exists(outsideInteractions[outInterIdx].outIndex)){
            typename CellContainerClass::CompleteCellClass interCell = cellsOther.getCompleteCell(outsideInteractions[outInterIdx].outIndex);
            FCudaAssertLF(interCell.symb->mortonIndex == outsideInteractions[outInterIdx].outIndex);
            typename CellContainerClass::CompleteCellClass cell = currentCells.getCompleteCell(outsideInteractions[outInterIdx].insideIndex);
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
                                               int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t 	currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );


    FCuda__transferInoutPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(currentCellsPtr, currentCellsSize,
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

    const MortonIndex blockStartIdx = FCudaMax(currentCells.getStartingIndex(),
                                          subCellGroups[0].getStartingIndex()>>3);
    const MortonIndex blockEndIdx   = FCudaMin(currentCells.getEndingIndex(),
                                          ((subCellGroups[nbSubCellGroups-1].getEndingIndex()-1)>>3)+1);

    int idxSubCellGroup = 0;

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
        if(currentCells.exists(mindex)){
            typename CellContainerClass::CompleteCellClass cell = currentCells.getDownCell(mindex);
            FCudaAssertLF(cell.symb->mortonIndex == mindex);
            typename CellContainerClass::CompleteCellClass child[8];

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if( subCellGroups[idxSubCellGroup].getEndingIndex() <= ((mindex<<3)+idxChild) ){
                    idxSubCellGroup += 1;
                }
                if( idxSubCellGroup == nbSubCellGroups ){
                    break;
                }
                if(subCellGroups[idxSubCellGroup].exists((mindex<<3)+idxChild)){
                    child[idxChild] = subCellGroups[idxSubCellGroup].getDownCell((mindex<<3)+idxChild);
                }
                else{
                    child[idxChild].symb = nullptr;
                }
            }

            kernel->L2L(cell, child, idxLevel);
        }
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__downardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
                                         FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
                                         FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
                                         int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t 	currentStream){

    FCuda__downardPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>
            (currentCellsPtr, currentCellsSize, currentCellsDownPtr, subCellGroupsPtr, subCellGroupsSize, subCellGroupsDownPtr,
             nbSubCellGroups, idxLevel, kernel);
    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));
}

/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__global__ void FCuda__directInoutPassPerformMpi(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                                 unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                                 const OutOfBlockInteraction* outsideInteractions,
                                                 int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel){
    ParticleContainerGroupClass containers(containersPtr, containersSize, containersDownPtr);
    ParticleContainerGroupClass containersOther(externalContainersPtr, externalContainersSize, nullptr);

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

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallbackMpi(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                                unsigned char* externalContainersPtr, std::size_t externalContainersSize,
                                                const OutOfBlockInteraction* outsideInteractions,
                                                int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t 	currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerformMpi
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize, containersDownPtr,
                                  externalContainersPtr, externalContainersSize,
                                  cuOutsideInteractions, nbOutsideInteractions, treeHeight, kernel);

    FCudaCheckAfterCall();
    FCudaCheck(cudaStreamSynchronize(currentStream));

    FCudaCheck(cudaFree(cuOutsideInteractions));
}

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

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInPassCallback(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                          const int treeHeight, CudaKernelClass* kernel, cudaStream_t 	currentStream){
    FCuda__directInPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize, containersDownPtr,
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

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__directInoutPassCallback(unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                             unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
                                             const OutOfBlockInteraction* outsideInteractions,
                                             int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t 	currentStream){
    OutOfBlockInteraction* cuOutsideInteractions;
    FCudaCheck( cudaMalloc(&cuOutsideInteractions,nbOutsideInteractions*sizeof(OutOfBlockInteraction)) );
    FCudaCheck( cudaMemcpy( cuOutsideInteractions, outsideInteractions, nbOutsideInteractions*sizeof(OutOfBlockInteraction),
                cudaMemcpyHostToDevice ) );

    FCuda__directInoutPassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(containersPtr, containersSize,containersDownPtr,
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

    const MortonIndex blockStartIdx = leafCells.getStartingIndex();
    const MortonIndex blockEndIdx = leafCells.getEndingIndex();

    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
        if(leafCells.exists(mindex)){
            typename CellContainerClass::CompleteCellClass cell = leafCells.getDownCell(mindex);
            FCudaAssertLF(cell.symb->mortonIndex == mindex);
            ParticleGroupClass particles = containers.template getLeaf<ParticleGroupClass>(mindex);
            FCudaAssertLF(particles.isAttachedToSomething());
            kernel->L2P(cell, &particles);
        }
    }
}

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
__host__ void FCuda__mergePassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
                                       unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
                                       CudaKernelClass* kernel, cudaStream_t 	currentStream){
    FCuda__mergePassPerform
            <SymboleCellClass, PoleCellClass, LocalCellClass,
            CellContainerClass, ParticleContainerGroupClass, ParticleGroupClass, CudaKernelClass>
            <<<nbCudaThreads, nbCudaBlocks, 0, currentStream>>>(leafCellsPtr, leafCellsSize,leafCellsDownPtr,
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
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
    unsigned char* containersPtr, std::size_t containersSize,
    FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__upwardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
int nbSubCellGroups, int idxLevel, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__transferInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__downardPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallbackMpi<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__directInPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FCudaEmptyKernel* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                        FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FCudaEmptyKernel* kernel, cudaStream_t 	currentStream);

template void FCuda__mergePassCallback<FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                       FCudaGroupOfParticles<0,0,int>, FCudaGroupAttachedLeaf<0,0,int>, FCudaEmptyKernel >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FCudaEmptyKernel* kernel, cudaStream_t currentStream);


template FCudaEmptyKernel* FCuda__BuildCudaKernel< FCudaEmptyKernel >(void* kernel);
template void FCuda__ReleaseCudaKernel< FCudaEmptyKernel >(FCudaEmptyKernel* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#include "../TestKernel/FCudaTestKernels.hpp"
#include "../TestKernel/FTestCellPOD.hpp"

template void FCuda__bottomPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsUpPtr,
unsigned char* containersPtr, std::size_t containersSize,
    FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__upwardPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsUpPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsUpPtr,
    int nbSubCellGroups, int idxLevel, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallbackMpi<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize, unsigned char* externalCellsUpPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__transferInPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    int idxLevel, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__transferInoutPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* currentCellsUpPtr, unsigned char* currentCellsDownPtr,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    unsigned char* externalCellsUpPtr, unsigned char* externalCellsDownPtr,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__downardPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* currentCellsPtr, std::size_t currentCellsSize, unsigned char* currentCellsDownPtr,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    FCudaParams<unsigned char*,9> subCellGroupsDownPtr,
    int nbSubCellGroups, int idxLevel, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallbackMpi<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__directInPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    const int treeHeight, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__directInoutPassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize, unsigned char* externalContainersDownPtr,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, FTestCudaKernels* kernel, cudaStream_t currentStream);

template void FCuda__mergePassCallback<FTestCellPODCore, FTestCellPODData, FTestCellPODData, FCudaGroupOfCells<FTestCellPODCore, FTestCellPODData, FTestCellPODData>,
                                        FCudaGroupOfParticles<0, 1, long long int>, FCudaGroupAttachedLeaf<0, 1, long long int>, FTestCudaKernels >
    (unsigned char* leafCellsPtr, std::size_t leafCellsSize, unsigned char* leafCellsDownPtr,
    unsigned char* containersPtr, std::size_t containersSize, unsigned char* containersDownPtr,
    FTestCudaKernels* kernel, cudaStream_t currentStream);

template FTestCudaKernels* FCuda__BuildCudaKernel<FTestCudaKernels>(void* kernel);
template void FCuda__ReleaseCudaKernel<FTestCudaKernels>(FTestCudaKernels* cukernel);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
