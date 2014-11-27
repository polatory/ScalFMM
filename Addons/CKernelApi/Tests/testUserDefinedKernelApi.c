// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
  * In this file we implement a example kernel which is simply printing information
  * about the cells and particles.
  * It is recommanded to compile and execute this code in order to understand the behavior
  * of the application.
  * We mark some part with "JUST-PUT-HERE:" to give instruction to user to create its own kernel.
  **/

// Include the FMM API (should be in a different folder for you)
#include "../Src/CScalfmmApi.h"

// Uncomment the next line to avoid the print in the kernel and verbosity
// #define NOT_TOO_MUCH_VERBOSE
#ifdef NOT_TOO_MUCH_VERBOSE
#define VerbosePrint(X...)
#else
#define VerbosePrint(X...) printf(X)
#endif

///////////////////////////////////////////////////////////////////////////
/// Cell struct and functions
///////////////////////////////////////////////////////////////////////////

// This represent a cell
struct MyCellDescriptor{
    int level;
    long long mortonIndex;
    int coord[3];
    double position[3];
    // JUST-PUT-HERE:
    // You local and multipole arrays
};

// This is our function that init a cell (struct MyCellDescriptor)
void* my_Callback_init_cell(int level, long long mortonIndex, int* coord, double* position){
    VerbosePrint("\tAllocating cell for level %d, morton index %lld, coord %d/%d/%d\n", level, mortonIndex, coord[0], coord[1], coord[2]);
    struct MyCellDescriptor* cellData = (struct MyCellDescriptor*)malloc(sizeof(struct MyCellDescriptor));
    memset(cellData, 0, sizeof(struct MyCellDescriptor));
    cellData->level       = level;
    cellData->mortonIndex = mortonIndex;
    memcpy(cellData->coord, coord, sizeof(int)*3);
    memcpy(cellData->position, position, sizeof(double)*3);
    // JUST-PUT-HERE:
    // Fill your structure
    return cellData;
}

// To dealloc a cell
void my_Callback_free_cell(void* cellData){
    free(cellData);
}


///////////////////////////////////////////////////////////////////////////
/// Kernel struct and functions
///////////////////////////////////////////////////////////////////////////

// This is the data passed to our kernel
struct MyData {
    // We simply strore the number of call the and particle indexes
    double* insertedPositions;
    int countP2M;
    int countM2M;
    int countM2L;
    int countL2L;
    int countL2P;
    int countP2PInner;
    int countP2P;
    // JUST-PUT-HERE:
    // everything your kernel will need for its computation
    // pre-computed matrices, etc....
};


// Our P2M
void my_Callback_P2M(void* cellData, int nbParticlesInLeaf, const int* particleIndexes, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countP2M += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    VerbosePrint("Cell morton %lld is doing P2M with %d particles :\n", my_cell->mortonIndex, nbParticlesInLeaf);
    int idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        VerbosePrint("\t particle idx %d position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
        // JUST-PUT-HERE:
        // Your real P2M computation
    }
}

void my_Callback_M2M(int level, void* cellData, int childPosition, void* childData, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countM2M += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    struct MyCellDescriptor* my_child = (struct MyCellDescriptor*) childData;

    int childFullPosition[3];
    Scalfmm_utils_parentChildPosition(childPosition, childFullPosition);

    VerbosePrint("Doing a M2M at level %d, between cells %lld and %lld (position %d/%d/%d)\n",
           level, my_cell->mortonIndex, my_child->mortonIndex,
           childFullPosition[0], childFullPosition[1], childFullPosition[2]);
    // JUST-PUT-HERE: your real M2M computation
}

void my_Callback_M2L(int level, void* cellData, int srcPosition, void* srcData, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countM2L += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    struct MyCellDescriptor* my_src_cell = (struct MyCellDescriptor*) srcData;

    int interactionFullPosition[3];
    Scalfmm_utils_interactionPosition(srcPosition, interactionFullPosition);

    VerbosePrint("Doing a M2L at level %d, between cells %lld and %lld (position %d/%d/%d)\n",
           level, my_cell->mortonIndex, my_src_cell->mortonIndex,
           interactionFullPosition[0], interactionFullPosition[1], interactionFullPosition[2]);
    // JUST-PUT-HERE: Your M2L
}

void my_Callback_L2L(int level, void* cellData, int childPosition, void* childData, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countL2L += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    struct MyCellDescriptor* my_child = (struct MyCellDescriptor*) childData;

    int childFullPosition[3];
    Scalfmm_utils_parentChildPosition(childPosition, childFullPosition);

    VerbosePrint("Doing a L2L at level %d, between cells %lld and %lld (position %d/%d/%d)\n",
           level, my_cell->mortonIndex, my_child->mortonIndex,
           childFullPosition[0], childFullPosition[1], childFullPosition[2]);
    // JUST-PUT-HERE: Your L2L
}

void my_Callback_L2P(void* cellData, int nbParticlesInLeaf, const int* particleIndexes, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countL2P += 1;

    struct MyCellDescriptor* my_cell = (struct MyCellDescriptor*) cellData;
    VerbosePrint("Cell morton %lld is doing L2P with %d particles :\n", my_cell->mortonIndex, nbParticlesInLeaf);
    int idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        VerbosePrint("\t particle idx %d position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
        // JUST-PUT-HERE: Your L2P
    }
}

void my_Callback_P2P(int nbParticlesInLeaf, const int* particleIndexes, int nbParticlesInSrc, const int* particleIndexesSrc, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countP2P += 1;

    VerbosePrint("Doing P2P between two leaves of %d particles and %d particles :\n", nbParticlesInLeaf, nbParticlesInSrc);
    int idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        VerbosePrint("\t Target >> particle idx %d position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
    }
    for(idxPart = 0 ; idxPart < nbParticlesInSrc ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexesSrc[idxPart]*3];
        VerbosePrint("\t Target >> Src idx %d position %e/%e%e\n", particleIndexesSrc[idxPart],
               position[0], position[1], position[2]);
    }

    // JUST-PUT-HERE:
    // Put one loop into the other to make all particles from source
    // interacting with the target particles
}

void my_Callback_P2PInner(int nbParticlesInLeaf, const int* particleIndexes, void* userData){
    struct MyData* my_data = (struct MyData*)userData;
    my_data->countP2PInner += 1;

    VerbosePrint("Doing P2P inside a leaf of %d particles :\n", nbParticlesInLeaf);
    int idxPart;
    for(idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
        double* position = &my_data->insertedPositions[particleIndexes[idxPart]*3];
        VerbosePrint("\t particle idx %d position %e/%e%e\n", particleIndexes[idxPart],
               position[0], position[1], position[2]);
        // JUST-PUT-HERE: another loop to have all particles interacting with
        // the other from the same leaf
    }
}

///////////////////////////////////////////////////////////////////////////
/// Main
///////////////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    // The properties of our tree
    int treeHeight = 4;
    double boxWidth = 1.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;

    // Create random particles
    int nbParticles = 10;
    int particleIndexes[nbParticles];
    double particleXYZ[nbParticles*3];
    {
        printf("Creating Particles:\n");
        int idxPart;
        for(idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            particleIndexes[idxPart] = idxPart;
            particleXYZ[idxPart*3]   = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
            particleXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
            particleXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
            VerbosePrint("\t %d] %e/%e/%e\n", idxPart, particleXYZ[idxPart*3], particleXYZ[idxPart*3+1], particleXYZ[idxPart*3+2]);
        }
    }

    // Init the handle
    scalfmm_handle handle = scalfmm_init(treeHeight, boxWidth, boxCenter,user_defined_kernel);
    // Insert particles
    printf("Inserting particles...\n");
    Scalfmm_insert_array_of_particles(handle, nbParticles, particleIndexes, particleXYZ);
    // Init cell
    printf("Initizalizing the cells:\n");
    Scalfmm_init_cell(handle, my_Callback_init_cell);

    // Init our callback struct
    struct User_Scalfmm_Kernel_Descriptor kernel;
    kernel.p2m = my_Callback_P2M;
    kernel.m2m = my_Callback_M2M;
    kernel.m2l = my_Callback_M2L;
    kernel.l2l = my_Callback_L2L;
    kernel.l2p = my_Callback_L2P;
    kernel.p2pinner = my_Callback_P2PInner;
    kernel.p2p = my_Callback_P2P;

    // Init the data to pass to all our callbacks
    struct MyData my_data;
    memset(&my_data, 0, sizeof(struct MyData));
    my_data.insertedPositions = particleXYZ;

    // Execute the FMM
    Scalfmm_execute_kernel(handle, kernel, &my_data);

    // Print the results store in our callback
    printf("There was %d P2M\n", my_data.countP2M);
    printf("There was %d M2M\n", my_data.countM2M);
    printf("There was %d M2L\n", my_data.countM2L);
    printf("There was %d L2L\n", my_data.countL2L);
    printf("There was %d L2P\n", my_data.countL2P);
    printf("There was %d P2PInner\n", my_data.countP2PInner);
    printf("There was %d P2P\n", my_data.countP2P);

    // Dealloc the handle
    Scalfmm_dealloc_handle(handle,my_Callback_free_cell);

    return 0;
}
