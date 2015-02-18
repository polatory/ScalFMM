// @SCALFMM_PRIVATE
#ifndef FCUDADEVICEWRAPPER_HPP
#define FCUDADEVICEWRAPPER_HPP


#include "../../Utils/FGlobal.hpp"
#include "../FOutOfBlockInteraction.hpp"

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__bottomPassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    CudaKernelClass* kernel);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__upwardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
    CudaKernelClass* kernel, int nbSubCellGroups, int idxLevel);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__transferInoutPassCallbackMpi(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    CudaKernelClass* kernel, int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__transferInPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    CudaKernelClass* kernel, int idxLevel);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__transferInoutPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* externalCellsPtr, std::size_t externalCellsSize,
    CudaKernelClass* kernel, int idxLevel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__downardPassCallback(unsigned char* currentCellsPtr, std::size_t currentCellsSize,
    unsigned char* subCellGroupsPtr[9], std::size_t subCellGroupsSize[9],
    CudaKernelClass* kernel, int nbSubCellGroups, int idxLevel);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__directInoutPassCallbackMpi(unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    CudaKernelClass* kernel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__directInPassCallback(unsigned char* containersPtr, std::size_t containersSize,
    CudaKernelClass* kernel, const int treeHeight);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__directInoutPassCallback(unsigned char* containersPtr, std::size_t containersSize,
    unsigned char* externalContainersPtr, std::size_t externalContainersSize,
    CudaKernelClass* kernel, const OutOfBlockInteraction* outsideInteractions,
    int nbOutsideInteractions, const int treeHeight);

template <class CellContainerClass, class ParticleContainerGroupClass, class ParticleGroupClass, class CudaKernelClass>
void FCuda__mergePassCallback(unsigned char* leafCellsPtr, std::size_t leafCellsSize,
    unsigned char* containersPtr, std::size_t containersSize,
    CudaKernelClass* kernel);

template <class CudaKernelClass>
CudaKernelClass* FCuda__BuildCudaKernel(void*);

template <class CudaKernelClass>
void FCuda__ReleaseCudaKernel(CudaKernelClass*);

#endif
