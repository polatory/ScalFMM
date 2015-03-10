/** This file contains the prototype for a kernel in opencl */

typedef long long int MortonIndex;


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

typedef struct OutOfBlockInteraction_t{
    MortonIndex outIndex;
    MortonIndex insideIndex;
    int outPosition;
} OutOfBlockInteraction;


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

struct Uptr8{
    __global unsigned char* ptrs[8];
};

struct Uptr9{
    __global unsigned char* ptrs[9];
};

struct size_t9{
    size_t v[9];
};


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/


__kernel void FOpenCL__bottomPassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,
                                         __global unsigned char* containersPtr, size_t containersSize,
                                         __global void* userkernel ){
}


/////////////////////////////////////////////////////////////////////////////////////
/// Upward Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__upwardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                  struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize,
                                  int nbSubCellGroups, int idxLevel, __global void* userkernel){
}


/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass Mpi
/////////////////////////////////////////////////////////////////////////////////////

__kernel  void FOpenCL__transferInoutPassPerformMpi(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                             __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                             int idxLevel, const __global OutOfBlockInteraction* outsideInteractions,
                                             int nbOutsideInteractions, __global void* userkernel){
}


/////////////////////////////////////////////////////////////////////////////////////
/// Transfer Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel  void FOpenCL__transferInPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                       int idxLevel, __global void* userkernel){
}

__kernel void FOpenCL__transferInoutPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                         __global unsigned char* externalCellsPtr, size_t externalCellsSize,
                                         int idxLevel, const __global OutOfBlockInteraction* outsideInteractions,
                                         int nbOutsideInteractions, __global void* userkernel){
}

/////////////////////////////////////////////////////////////////////////////////////
/// Downard Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__downardPassPerform(__global unsigned char* currentCellsPtr, size_t currentCellsSize,
                                   struct Uptr9 subCellGroupsPtr, struct size_t9 subCellGroupsSize,
                                   int nbSubCellGroups, int idxLevel, __global void* userkernel){
}



/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass MPI
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__directInoutPassPerformMpi(__global unsigned char* containersPtr, size_t containersSize,
                                          __global unsigned char* externalContainersPtr, size_t externalContainersSize,
                                          const __global OutOfBlockInteraction* outsideInteractions,
                                          int nbOutsideInteractions, const int treeHeight, __global void* userkernel){
}


/////////////////////////////////////////////////////////////////////////////////////
/// Direct Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__directInPassPerform(__global unsigned char* containersPtr, size_t containersSize,
                                    const int treeHeight, __global void* userkernel){
}

__kernel void FOpenCL__directInoutPassPerform(__global unsigned char* containersPtr, size_t containersSize,
                                       __global unsigned char* externalContainersPtr, size_t externalContainersSize,
                                       const __global OutOfBlockInteraction* outsideInteractions,
                                       int nbOutsideInteractions, const int treeHeight, __global void* userkernel){
}

/////////////////////////////////////////////////////////////////////////////////////
/// Merge Pass
/////////////////////////////////////////////////////////////////////////////////////

__kernel void FOpenCL__mergePassPerform(__global unsigned char* leafCellsPtr, size_t leafCellsSize,
                                 __global unsigned char* containersPtr, size_t containersSize,
                                 __global void* userkernel){
}

