// ===================================================================================
// Copyright ScalFmm 2014 I
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
//
#ifndef CKERNELAPI_H
#define CKERNELAPI_H

/**
 * This file defines the API for the USER.
 * We briefly comment all the functions.
 * The objective of the C Kernel API is to give a quick and easy way
 * to anyone (who can program in C) to implement a kernel.
 * Using C++ is advised but this is a simple alternative.
 */


///////////////////////////////////////////////////////////////////////////
/// Init part
///////////////////////////////////////////////////////////////////////////

//< For the user an handle is a void*
typedef void* Scalfmm_Handle;

//< Function to init the cells (should be given by the user when calling Scalfmm_init_cell)
//< it gives the level of the cell, its morton index, it position in term of box at that level
//< and the spatial position of its center
typedef void* (*Callback_init_cell)(int level, long long morton_index, int* tree_position, double* spatial_position);
//< Function to destroy what have bee initialized by the user (should be give in Scalfmm_dealloc_handle)
typedef void (*Callback_free_cell)(void*);

//< This function init an handle (and an tree based on the given properties)
Scalfmm_Handle Scalfmm_init_handle(int treeHeight, double boxWidth, double* boxCenter);
//< This function should be used to dealloc our handle
void Scalfmm_dealloc_handle(Scalfmm_Handle handle, Callback_free_cell cellDestroyer);

//< This function should be used to insert an array of particle in the tree
//< The indexes are the one used on the particles operator
//< The position of the particles should be composed of one triple per particle:
//< xyzxyzxyz...
void Scalfmm_insert_array_of_particles(Scalfmm_Handle handle, int nbParticles, int* particleIndexes, double* particleXYZ);
//< To insert one particle only
void Scalfmm_one_particle(Scalfmm_Handle handle, int particleIndexe, double x, double y, double z);

//< This function should be called to init the cells
//< It must be called after all the particles have been inserted!
void Scalfmm_init_cell(Scalfmm_Handle handle, Callback_init_cell cellInitializer);


///////////////////////////////////////////////////////////////////////////
/// Kernel part
///////////////////////////////////////////////////////////////////////////

//< These function are the callbacks of the FMM operators
typedef void (*Callback_P2M)(void* leafCell, int nbParticles, const int* particleIndexes, void* userData);
typedef void (*Callback_M2M)(int level, void* parentCell, int childPosition, void* childCell, void* userData);
typedef void (*Callback_M2L)(int level, void* targetCell, int sourceCellPosition, void* sourceCell, void* userData);
typedef void (*Callback_L2L)(int level, void* parentCell, int childPosition, void* childCell, void* userData);
typedef void (*Callback_L2P)(void* leafCell, int nbParticles, int* particleIndexes, void* userData);
typedef void (*Callback_P2P)(int nbParticles, const int* particleIndexes, int nbSourceParticles, const int* sourceParticleIndexes, void* userData);
typedef void (*Callback_P2PInner)(int nbParticles, int* particleIndexes, void* userData);

//< This structure should be filled (or filled with null) to call the FMM
struct Scalfmm_Kernel_Descriptor {
    Callback_P2M p2m;
    Callback_M2M m2m;
    Callback_M2L m2l;
    Callback_L2L l2l;
    Callback_L2P l2p;
    Callback_P2P p2p;
    Callback_P2PInner p2pinner;
};

//< Execute one FMM using the given kernel, the userData is the one given in all the operators as last
//< parameter.
void Scalfmm_execute_kernel(Scalfmm_Handle handle, struct Scalfmm_Kernel_Descriptor userKernel, void* userData);


///////////////////////////////////////////////////////////////////////////
/// Util functions
///////////////////////////////////////////////////////////////////////////


//< This function fill the childFullPosition[3] with [0;1] to know the position of a child relatively to
//< its position from its parent
inline void Scalfmm_utils_parentChildPosition(int childPosition, int* childFullPosition){
    childFullPosition[0] = childPosition%2;
    childFullPosition[1] = (childPosition/2)%2;
    childFullPosition[2] = (childPosition/4)%2;
}

//< This function fill the childFullPosition[3] with [-3;3] to know the position of a interaction
//< cell relatively to its position from the target
inline void Scalfmm_utils_interactionPosition(int interactionPosition, int* srcPosition){
    srcPosition[0] = interactionPosition%7 - 3;
    srcPosition[1] = (interactionPosition/7)%7 - 3;
    srcPosition[2] = (interactionPosition/49)%7 - 3;
}

#endif // CKERNELAPI_H
