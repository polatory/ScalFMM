// @SCALFMM_PRIVATE
#ifndef FEMPTYOPENCLCODE_HPP
#define FEMPTYOPENCLCODE_HPP

// Return the same thing as FEmptyKernel.cl

#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

class FEmptyOpenCLCode{

public:
    FEmptyOpenCLCode(){
    }

    const char* getKernelCode(const int /*inDevId*/){
        const char* kernelcode =
"typedef long long int MortonIndex; \
#define DefaultStructAlign " FStarPUDefaultAlignStr "\
typedef struct OutOfBlockInteraction{ \
    MortonIndex outIndex; \
    MortonIndex insideIndex; \
    int outPosition; \
} __attribute__ ((aligned (DefaultStructAlign))); \
struct Uptr9{ \
    __global unsigned char* ptrs[9]; \
} __attribute__ ((aligned (DefaultStructAlign))); \
struct size_t9{ \
    size_t v[9]; \
}__attribute__ ((aligned (DefaultStructAlign))); \
__kernel void FOpenCL__bottomPassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize, \
                                         __global unsigned char* containersPtr, size_t containersSize, \
                                         __global void* userkernel ){ \
} \
__kernel void FOpenCL__upwardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize, \
                                  struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize, \
                                  int nbSubCellGroups, int idxLevel, __global void* userkernel){ \
} \
__kernel  void FOpenCL__transferInoutPassPerformMpi(__global unsigned char* currentCellsPtr, size_t currentCellsSize, \
                                             __global unsigned char* externalCellsPtr, size_t externalCellsSize, \
                                             int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions, \
                                             size_t nbOutsideInteractions, __global void* userkernel){ \
} \
__kernel  void FOpenCL__transferInPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize, \
                                       int idxLevel, __global void* userkernel){ \
} \
__kernel void FOpenCL__transferInoutPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize, \
                                         __global unsigned char* externalCellsPtr, size_t externalCellsSize, \
                                         int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions, \
                                         size_t nbOutsideInteractions, __global void* userkernel){ \
} \
__kernel void FOpenCL__downardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize, \
                                   struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize, \
                                   int nbSubCellGroups, int idxLevel, __global void* userkernel){ \
} \
__kernel void FOpenCL__directInoutPassPerformMpi(__global unsigned char* containersPtr, size_t containersSize, \
                                          __global unsigned char* externalContainersPtr, size_t externalContainersSize, \
                                          const __global struct OutOfBlockInteraction* outsideInteractions, \
                                          size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){ \
} \
__kernel void FOpenCL__directInPassPerform(__global unsigned char* containersPtr, size_t containersSize, \
                                    const int treeHeight, __global void* userkernel){ \
} \
__kernel void FOpenCL__directInoutPassPerform(__global unsigned char* containersPtr, size_t containersSize, \
                                       __global unsigned char* externalContainersPtr, size_t externalContainersSize, \
                                       const __global struct OutOfBlockInteraction* outsideInteractions, \
                                       size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){ \
} \
__kernel void FOpenCL__mergePassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize, \
                                 __global unsigned char* containersPtr, size_t containersSize, \
                                 __global void* userkernel){ \
}";
        return kernelcode;
    }

    void releaseKernelCode(){
    }

    size_t getNbDims() const {
        return 0;
    }

    const size_t* getNbGroups() const {
        return nullptr;
    }

    const size_t* getGroupSize() const {
        return nullptr;
    }
};

#endif // FEMPTYOPENCLCODE_HPP

