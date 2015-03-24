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
#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FBasicKernels.hpp"
#include "../../Src/Components/FBasicParticleContainer.hpp"

/** It should be compiled with C export */
extern "C" {
#include "CKernelApi.h"
}

/**
 * This file contains the code to have a C Kernel API
 * working. It is simply a wrapper with light cell and
 * kernel in order to delegate the calls to the API.
 *
 * The user should allocate the cell data (and deallocate them).
 */

/**
 * The core cell is a basic cell + data from the user.
 * The cell is not responsible to allocate or deallocate the user data.
 */
class CoreCell : public FBasicCell {
    // Mutable in order to work with the API
    mutable void* userData;

public:
    CoreCell() : userData(nullptr) {
    }

    ~CoreCell(){
    }

    void setContainer(void* inContainer) const {
        userData = inContainer;
    }

    void* getContainer() const {
        return userData;
    }
};


/**
 * This class simply call the function pointers from Scalfmm_Kernel_Descriptor.
 * If not pointer is set the calls are skipped.
 * The userData is given at any calls.
 */
template< class CellClass, class ContainerClass>
class CoreKernel : public FAbstractKernels<CellClass,ContainerClass> {
    struct Scalfmm_Kernel_Descriptor kernel;
    void* userData;

public:
    CoreKernel(struct Scalfmm_Kernel_Descriptor inKernel, void* inUserData) : kernel(inKernel) , userData(inUserData){
    }

    /** Default destructor */
    virtual ~CoreKernel(){
    }

    /** Do nothing */
    virtual void P2M(CellClass* const cell, const ContainerClass* const container) {
        if(kernel.p2m) kernel.p2m(cell->getContainer(), container->getNbParticles(), container->getAttribute(0), userData);
    }

    /** Do nothing */
    virtual void M2M(CellClass* const FRestrict cell, const CellClass*const FRestrict *const FRestrict children, const int level) {
        if(kernel.m2m){
            for(int idx = 0 ; idx < 8 ; ++idx){
                if( children[idx] ){
                    kernel.m2m(level, cell->getContainer(), idx, children[idx]->getContainer(), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void M2L(CellClass* const FRestrict cell, const CellClass* interactions[], const int , const int level) {
        if(kernel.m2l){
            for(int idx = 0 ; idx < 343 ; ++idx){
                if( interactions[idx] ){
                    kernel.m2l(level, cell->getContainer(), idx, interactions[idx]->getContainer(), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void L2L(const CellClass* const FRestrict cell, CellClass* FRestrict *const FRestrict children, const int level) {
        if(kernel.l2l){
            for(int idx = 0 ; idx < 8 ; ++idx){
                if( children[idx] ){
                    kernel.l2l(level, cell->getContainer(), idx, children[idx]->getContainer(), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void L2P(const CellClass* const cell, ContainerClass* const container){
        if(kernel.l2p) kernel.l2p(cell->getContainer(), container->getNbParticles(), container->getAttribute(0), userData);
    }


    /** Do nothing */
    virtual void P2P(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
                     ContainerClass* const neighbors[27], const int ){
        if(kernel.p2pinner) kernel.p2pinner(targets->getNbParticles(), targets->getAttribute(0), userData);

        if(kernel.p2p){
            for(int idx = 0 ; idx < 27 ; ++idx){
                if( neighbors[idx] ){
                    kernel.p2p(targets->getNbParticles(), targets->getAttribute(0),
                                    neighbors[idx]->getNbParticles(), neighbors[idx]->getAttribute(0), userData);
                }
            }
        }
    }

    /** Do nothing */
    virtual void P2PRemote(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [27], const int ){
    }

};

// Here are all the type definition
typedef FBasicParticleContainer<1, int>   CoreContainerClass;
typedef FSimpleLeaf<FReal, CoreContainerClass >                        LeafClass;
typedef FOctree<CoreCell, CoreContainerClass , LeafClass >     OctreeClass;
typedef CoreKernel<CoreCell, CoreContainerClass>         CoreKernelClass;
typedef FFmmAlgorithmThread<OctreeClass, CoreCell, CoreContainerClass, CoreKernelClass, LeafClass >     FmmClassThread;


// Our scalfmm handle
struct ScalFmmCoreHandle {
    struct ScalFmmCoreConfig {
        // Read/Write parameter
        int treeHeight;     //  Number of level in the octree
        FReal boxWidth;    // Simulation box size (root level)
        FPoint<FReal> boxCenter; // Center position of the box simulation(FReal[3])
    };

    ScalFmmCoreConfig config;
    OctreeClass* octree;
};


// Allocate and init the handle
extern "C" Scalfmm_Handle Scalfmm_init_handle(int treeHeight, double boxWidth, double* boxCenter){
    ScalFmmCoreHandle* corehandle = new ScalFmmCoreHandle;
    memset(corehandle, 0, sizeof(ScalFmmCoreHandle));

    corehandle->config.treeHeight = treeHeight;
    corehandle->config.boxWidth   = boxWidth;
    corehandle->config.boxCenter  = FPoint<FReal>(boxCenter[0],boxCenter[1],boxCenter[2]);

    corehandle->octree = new OctreeClass(corehandle->config.treeHeight, FMath::Min(3,corehandle->config.treeHeight-1),
                                         corehandle->config.boxWidth, corehandle->config.boxCenter);

    return corehandle;
}

// Init the cells (once the particles have been pushed)
extern "C" void Scalfmm_init_cell(Scalfmm_Handle handle, Callback_init_cell cellInitializer){
    if(cellInitializer){
        ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)handle;

        double boxWidth = corehandle->config.boxWidth;
        double boxCorner[3];
        boxCorner[0] = corehandle->config.boxCenter.getX() - boxWidth/2.0;
        boxCorner[1] = corehandle->config.boxCenter.getY() - boxWidth/2.0;
        boxCorner[2] = corehandle->config.boxCenter.getZ() - boxWidth/2.0;

        corehandle->octree->forEachCellWithLevel([&](CoreCell* cell, int level){
            int coord[3] = {cell->getCoordinate().getX(),cell->getCoordinate().getY(),cell->getCoordinate().getZ() };

            double position[3];
            position[0] = boxCorner[0] + coord[0]*boxWidth/double(1<<level);
            position[1] = boxCorner[1] + coord[1]*boxWidth/double(1<<level);
            position[2] = boxCorner[2] + coord[2]*boxWidth/double(1<<level);

            cell->setContainer(cellInitializer(level, cell->getMortonIndex(), coord, position));
        });
    }
}

// Dealloc the handle and call the destroy function on the cells if given
extern "C" void Scalfmm_dealloc_handle(Scalfmm_Handle handle, Callback_free_cell cellDestroyer){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)handle;

    if(cellDestroyer){
        corehandle->octree->forEachCell([&](CoreCell* cell){
            cellDestroyer(cell->getContainer());
        });
    }

    delete corehandle->octree;
    delete corehandle;
}

// Insert a full array of particles. Posisions are xyz.xyz.xyz...
extern "C" void Scalfmm_insert_array_of_particles(Scalfmm_Handle handle, int nbParticles, int* particleIndexes, double* particleXYZ){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)handle;
    for(int idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
        corehandle->octree->insert(FPoint<FReal>(particleXYZ[idxPart*3],particleXYZ[idxPart*3+1],particleXYZ[idxPart*3+2]),
                                particleIndexes[idxPart]);
    }
}

// Push one particle in the tree
extern "C" void Scalfmm_one_particle(Scalfmm_Handle handle, int particleIndex, double x, double y, double z){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)handle;
    corehandle->octree->insert(FPoint<FReal>(x,y,z), particleIndex);
}

// Execute a kernel (by using the user functions callback).
extern "C" void Scalfmm_execute_kernel(Scalfmm_Handle handle, struct Scalfmm_Kernel_Descriptor userKernel, void *userData){
    ScalFmmCoreHandle* corehandle = (ScalFmmCoreHandle*)handle;

    CoreKernelClass kernel(userKernel, userData);
    FmmClassThread algorithm(corehandle->octree,&kernel);

    algorithm.execute();
}

//< This function fill the childFullPosition[3] with [-1;1] to know the position of a child relatively to
//< its position from its parent
extern "C" void Scalfmm_utils_parentChildPosition(int childPosition, int* childFullPosition){
    childFullPosition[2] = (childPosition%2 ? 1 : -1);
    childFullPosition[1] = ((childPosition/2)%2 ? 1 : -1);
    childFullPosition[0] = ((childPosition/4)%2 ? 1 : -1);
}

//< This function fill the childFullPosition[3] with [-3;3] to know the position of a interaction
//< cell relatively to its position from the target
extern "C" void Scalfmm_utils_interactionPosition(int interactionPosition, int* srcPosition){
    srcPosition[2] = interactionPosition%7 - 3;
    srcPosition[1] = (interactionPosition/7)%7 - 3;
    srcPosition[0] = (interactionPosition/49)%7 - 3;
}


