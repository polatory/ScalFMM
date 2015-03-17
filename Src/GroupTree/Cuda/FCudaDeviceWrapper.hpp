// @SCALFMM_PRIVATE
#ifndef FCUDADEVICEWRAPPER_HPP
#define FCUDADEVICEWRAPPER_HPP


#include "../../Utils/FGlobal.hpp"
#include "../Core/FOutOfBlockInteraction.hpp"
#include "FCudaStructParams.hpp"

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__bottomPassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__upwardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t, 9> subCellGroupsSize,
    int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__transferInoutPassCallbackMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__transferInPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    int idxLevel, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__transferInoutPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__downardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    FCudaParams<unsigned char*,9> subCellGroupsPtr, FCudaParams<std::size_t,9> subCellGroupsSize,
    int nbSubCellGroups, int idxLevel, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__directInoutPassCallbackMpi(unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__directInPassCallback(unsigned char* containersPtr, std::size_t containersSize,
    const int treeHeight, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__directInoutPassCallback(unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight, CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CellClass, class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__mergePassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    CudaKernelClass* kernel, cudaStream_t 	currentStream);

template <class CudaKernelClass>
CudaKernelClass* FCuda__BuildCudaKernel(void*);

template <class CudaKernelClass>
void FCuda__ReleaseCudaKernel(CudaKernelClass*);

#endif
