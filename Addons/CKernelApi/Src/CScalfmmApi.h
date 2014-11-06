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
 * @file
 * This file defines the API for the USER.  The objective is to
 * provide an unified way to compute the Fast Multipole Method with
 * ScalFMM's algorithm. The mathematical kernel used can be a user's
 * one or ScalFMM implementation of Chebyshev Interpolation or
 * Lagrange Interpolation.
 */



/////////////////////////////////////////////////////////////////////
//////////////////     Init  Part                      //////////////
/////////////////////////////////////////////////////////////////////

/**
 * Enum over the different kernel usable
 */
typedef enum kernel_type {
    user_defined_kernel = 0, /** Case if user provides a complete
			      * kernel, ie all the needed function to
			      * compute FMM */
    chebyshev = 1, /** Case if user uses ScalFMM's implementation of
		    *  Chebyshev Interpolation */
    lagrange = 2,  /** Case if user uses ScalFMM's implementation of
		    *  Lagrange Interpolation*/
} scalfmm_kernel_type;


/**
 * Handle for the user
 */
typedef void* scalfmm_handle;



/**
 * @brief This function initialize scalfmm datas.
 * @param TreeHeight Height of the octree.
 * @param BoxWidth Width of the entire simulation box.
 * @param BoxCenter Coordinate of the center of the box (ie array)
 * @param KernelType kernel to be used
 * @return scalfmm_handle (ie void *). This handle will be given to
 * every other scalfmm functions.
 *
 * Every data will be stored in order to crete later through a builder
 * what is needed for the simulation
 */
scalfmm_handle scalfmm_init(int TreeHeight,double BoxWidth,double* BoxCenter, scalfmm_kernel_type KernelType);


/////////////////////////////////////////////////////////////////////
//////////////////       Tree  Part                    //////////////
/////////////////////////////////////////////////////////////////////

/**
 * @brief This function insert an array of position into the octree
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param NbPositions Number of position to be inserted
 * @param arrayX Array containing the X coordinate for all the parts, size : NbPositions
 * @param arrayY Array containing the Y coordinate for all the parts, size : NbPositions
 * @param arrayZ Array containing the Z coordinate for all the parts, size : NbPositions
 *
 * The parts will be inserted with their indices inside the
 * array. Each index will be unique.
 */
void scalfmm_tree_insert_arrays(scalfmm_handle Handle, int NbPositions, double * arrayX, double * arrayY, double * arrayZ);


/**
 * @brief This function insert an array of position into the octree
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param NbPositions Number of position to be inserted
 * @param XYZ Array containing the X,Y,Z coordinates for all the
 * parts, sequentially stored. Size : NbPositions*3
 *
 * The parts will be inserted with their indices inside the
 * array. Each index will be unique.
 */
void scalfmm_tree_insert_positions(scalfmm_handle Handle, int NbPositions, double * XYZ);



/**
 * @brief This function insert an array of physical values into the octree
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbPhysicalValues Number of physical values to be
 * inserted. Must be equal to the number of positions previously
 * inserted.
 * @param physicalValues Array containing the physical values to be
 * associated to each parts.
 *
 * The physical values will be stored according to their indices in
 * the array
 *
 */
void scalfmm_set_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues);


/**
 * @brief This function give back the resulting forces
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbParts total number of particles to retrieve (must be equal
 * to the number of parts inserted)
 * @param forcesToFill array of size nbParts*3, that will contains the
 * forces. WARNING : User must allocate the array before call.
 *
 * Forces will be stored sequentially, according to the indices in the
 * array. (ie fx1,fy1,fz1,fx2,fy2,fz2,fx3 ....)
 */
void scalfmm_get_forces_xyz(scalfmm_handle Handle, int nbParts, double * forcesToFill);


/**
 * @brief This function give back the resulting forces
 * @param Handle scalfmm_handle provided by scalfmm_init.
 * @param nbParts total number of particles to retrieve (must be equal
 * to the number of parts inserted)
 * @param forcesX array of size nbParts, that will contains the
 * forces . WARNING : User must allocate the array before call.
 * @param forcesY array of size nbParts, that will contains the
 * forces . WARNING : User must allocate the array before call.
 * @param forcesZ array of size nbParts, that will contains the
 * forces . WARNING : User must allocate the array before call.
 *
 */
void scalfmm_get_forces(scalfmm_handle Handle, int nbParts, double * forcesX, double * forcesY, double * forcesZ);



/**
 * @brief This function give back the resulting potentials
 * @param Handle scalfmm_handle provided by scalfmm_init
 * @param nbParts total number of particles to retrieve (must be equal
 * to the number of parts inserted)
 * @param potentialsToFill array of potentials to be filled. WARNING :
 * User must allocate the array before call.
 *
 * Potentials will be stored sequentially, according to the indices in the
 * array.
 */
void scalfmm_get_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToFill);


/**
 * @brief This function update the positions inside the tree, in case
 * of multiple runs of the FMM.
 * @param Handle scalfmm_handle provided by scalfmm_init
 * @param NbPositions total number of particles (must be equal to the
 * number of parts inserted)
 * @param updatedXYZ array of displacement (ie
 * dx1,dy1,dz1,dx2,dy2,dz2,dx3 ...)
 */
void scalfmm_update_positions(scalfmm_handle Handle, int NbPositions, double * updatedXYZ)


/**
 * @brief This function set again the positions inside the tree, in case
 * of multiple runs of the FMM.
 * @param Handle scalfmm_handle provided by scalfmm_init
 * @param NbPositions total number of particles (must be equal to the
 * number of parts inserted)
 * @param newXYZ array of new positions (ie
 * dx1,dy1,dz1,dx2,dy2,dz2,dx3 ...)
 *
 * @return Error code, a parts may move out of the simulation
 * box. ScalFMM cannot deals with that specific case. Error code :
 * 0. Success code 1. Could be an arg in order to be Fortran
 * compliant.
 *
 */
int scalfmm_change_positions(scalfmm_handle Handle, int NbPositions, double * newXYZ)



/////////////////////////////////////////////////////////////////////
//////////////////       Kernel  Part                  //////////////
/////////////////////////////////////////////////////////////////////


///////////////// User kernel part :


/**
 * @brief Function to be filled by user's P2M
 * @param nbParticles number of particle in current leaf
 * @param leafCell current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2M)(void* leafCell, int nbParticles, const int* particleIndexes, void* userData);

/**
 * @brief Function to be filled by user's M2M
 * @param level current level in the tree
 * @prama parentCell cell to be filled
 * @param childPosition number of child (child position in the tree can inferred from its number (refer to doc))
 * @param userData datas specific to the user's kernel
 * @param childCell array of cells to be read
 */
typedef void (*Callback_M2M)(int level, void* parentCell, int childPosition, void* childCell, void* userData);

/**
 * @brief Function to be filled by user's M2L
 * @param level current level in the tree
 * @param targetCell pointer to cell to be filled
 * @param sourceCellPosition number of source cell (source cell
 * position in the tree can inferred from its number (refer to doc))
 * @param sourceCell array of cell to be read
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_M2L)(int level, void* targetCell, int sourceCellPosition, void* sourceCell, void* userData);

/**
 * @brief Function to be filled by user's L2L
 * @param level current level in the tree
 * @param parentCell cell to be read
 * @param childPosition number of child (child position in the tree can inferred from its number (refer to doc))
 * @param childCell cell to be filled
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_L2L)(int level, void* parentCell, int childPosition, void* childCell, void* userData);

/**
 * @brief Function to be filled by user's L2P
 * @param leafCell leaf to be filled
 * @param nbParticles number of particles in the current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_L2P)(void* leafCell, int nbParticles, int* particleIndexes, void* userData);

/**
 * @brief Function to be filled by user's P2P
 * @param nbParticles number of particle in current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param nbSourceParticles number of particles in source leaf
 * @param sourceParticleIndexes indexes of cource particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2P)(int nbParticles, const int* particleIndexes, int nbSourceParticles, const int* sourceParticleIndexes, void* userData);

/**
 * @brief Function to be filled by user's P2P inside the leaf
 * @param nbParticles number of particle in current leaf
 * @param particleIndexes indexes of particles currently computed
 * @param userData datas specific to the user's kernel
 */
typedef void (*Callback_P2PInner)(int nbParticles, int* particleIndexes, void* userData);



/**
 * @brief Structure containing callbacks to fill in order to define
 * user kernel.
 *
 */
struct Scalfmm_Kernel_Descriptor {
    Callback_P2M p2m;
    Callback_M2M m2m;
    Callback_M2L m2l;
    Callback_L2L l2l;
    Callback_L2P l2p;
    Callback_P2P p2p;
    Callback_P2PInner p2pinner;
};

//To fill gestion des cellules utilisateurs




///////////////// Common kernel part : /////////////////


/**
 * @brief This function launch the fmm on the parameters given
 * @param Handle scalfmm_handle provided by scalfmm_init

 * @param user_kernel Scalfmm_Kernel_Descriptor containing callbacks
 * to user Fmm function. Meaningless if using one of our kernel
 * (Chebyshev or Interpolation).

 * @param userDatas Data that will be passed to each FMM
 * function. Can be anything, but allocation/deallocation is user
 * side.
 *
 */
void scalfmm_execute_kernel(scalfmm_handle Handle, Scalfmm_Kernel_Descriptor userKernel, void * userDatas)
