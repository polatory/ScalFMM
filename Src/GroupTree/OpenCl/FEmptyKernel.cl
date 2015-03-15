/** This file contains the prototype for a kernel in opencl */
// @SCALFMM_PRIVATE

typedef long long int MortonIndex;

#define DefaultStructAlign ___DefaultStructAlign___

typedef struct OutOfBlockInteraction{
    MortonIndex outIndex;
    MortonIndex insideIndex;
    int outPosition;
} __attribute__ ((aligned (DefaultStructAlign)));
struct Uptr9{
    __global unsigned char* ptrs[9];
} __attribute__ ((aligned (DefaultStructAlign)));
struct size_t9{
    size_t v[9];
}__attribute__ ((aligned (DefaultStructAlign)));
__kernel void FOpenCL__bottomPassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,
                                         __global unsigned char* containersPtr, size_t containersSize,
                                         __global void* userkernel ){
}
__kernel void FOpenCL__upwardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                  struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize,
                                  int nbSubCellGroups, int idxLevel, __global void* userkernel){
}
__kernel  void FOpenCL__transferInoutPassPerformMpi(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                             __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                             int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions,
                                             size_t nbOutsideInteractions, __global void* userkernel){
}
__kernel  void FOpenCL__transferInPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                       int idxLevel, __global void* userkernel){
}
__kernel void FOpenCL__transferInoutPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                         __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                         int idxLevel, const __global struct OutOfBlockInteraction* outsideInteractions,
                                         size_t nbOutsideInteractions, __global void* userkernel){
}
__kernel void FOpenCL__downardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                   struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize,
                                   int nbSubCellGroups, int idxLevel, __global void* userkernel){
}
__kernel void FOpenCL__directInoutPassPerformMpi(__global unsigned char* containersPtr, size_t containersSize,
                                          __global unsigned char* externalContainersPtr, size_t externalContainersSize,
                                          const __global struct OutOfBlockInteraction* outsideInteractions,
                                          size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){
}
__kernel void FOpenCL__directInPassPerform(__global unsigned char* containersPtr, size_t containersSize,
                                    const int treeHeight, __global void* userkernel){
}
__kernel void FOpenCL__directInoutPassPerform(__global unsigned char* containersPtr, size_t containersSize,
                                       __global unsigned char* externalContainersPtr, size_t externalContainersSize,
                                       const __global struct OutOfBlockInteraction* outsideInteractions,
                                       size_t nbOutsideInteractions, const int treeHeight, __global void* userkernel){
}
__kernel void FOpenCL__mergePassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,
                                 __global unsigned char* containersPtr, size_t containersSize,
                                 __global void* userkernel){
}

