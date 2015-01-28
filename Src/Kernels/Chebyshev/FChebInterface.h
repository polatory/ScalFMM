// ===================================================================================
// Copyright ScalFmm 2011 INRIA,
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================


/**
 * This file provide an interface to the Chebyshev kernel, in order to
 * call it from C code (and thus use it through API's user defined
 * kernel feature).
 */

#ifndef FCHEBINTERFACE_H
#define FCHEBINTERFACE_H



//To access a cell
struct FChebCell_struct;
typedef struct FChebCell_struct ChebCellStruct;
ChebCellStruct * ChebCellStruct_create(long long int index,int * tree_position);
void ChebCellStruct_free(ChebCellStruct * cell);

/* //To access a Leaf */
/* struct FSimpleLeaf_struct; */
/* typedef struct FSimpleLeaf_struct ChebLeafStruct; */
/* ChebLeafStruct * ChebLeafStruct_create(); */
/* void ChebLeafStruct_free(ChebLeafStruct* leaf); */

//To access the kernel
struct FChebKernel_struct;
typedef struct FChebKernel_struct ChebKernelStruct;
ChebKernelStruct * ChebKernelStruct_create(int inTreeHeight,
                                           double inBoxWidth,
                                           double* inBoxCenter);

void ChebKernelStruct_free(void * kernel);
//To access kernel member function

void ChebKernel_P2M(void * leafCell, int nbParticles,const int* particleIndexes, void* kernel);
void ChebKernel_M2M(int level, void* parentCell, int childPosition, void* childCell, void* kernel);
void ChebKernel_M2L(int level, void* targetCell, void* sourceCell[343], void* kernel);
void ChebKernel_L2L(int level, void* parentCell, int childPosition, void* childCell, void* kernel);
void ChebKernel_L2P(void* leafCell, int nbParticles, const int* particleIndexes, void* kernel);
void ChebKernel_P2P(int nbParticles, const int* particleIndexes,
                    const int * sourceParticleIndexes[27], int sourceNbPart[27],void* userData);

typedef struct myUserDatas{
    ChebKernelStruct * kernelStruct;
    double * insertedPositions;
    double * myPhyValues;
    double ** forcesComputed;
}UserData;


#endif //FCHEBINTERFACE_H
