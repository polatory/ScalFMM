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


/**
 * @file This file contains a class that inherits from FScalFMMEngine,
 * and will implement the API functions for Interpolations kernels.
 */

#ifndef FINTERENGINE_HPP
#define FINTERENGINE_HPP

#include "FScalFMMEngine.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
//#include "Kernels/P2P/FP2PLeafInterface.hpp"
#include "Arranger/FOctreeArranger.hpp"
#include "Arranger/FArrangerPeriodic.hpp"
#include "Arranger/FBasicParticleContainerIndexedMover.hpp"

#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmPeriodic.hpp"


/**
 * @class FInterEngine implements API for Interpolations kernels, its
 * templates can be ChebCell/ChebKernel or UnifCell/UnifKernel
 */
template<class InterCell,class InterKernel,
         class ContainerClass = FP2PParticleContainerIndexed<FReal>,
         class LeafClass = FSimpleLeaf<FReal, FP2PParticleContainerIndexed<FReal> >,
         class MatrixKernelClass = FInterpMatrixKernelR>
class FInterEngine : public FScalFMMEngine{
private:
    //Typedef on the octree class, in order to clarify following code
    typedef FOctree<InterCell,ContainerClass,LeafClass> OctreeClass;
    //typedef FP2PLeafInterface<OctreeClass>            LeafInterface;


    //Typedef on Octree Arranger, in order to clarify following code
    typedef FBasicParticleContainerIndexedMover<OctreeClass, ContainerClass> MoverClass;
    typedef FOctreeArranger<FReal,OctreeClass, ContainerClass, MoverClass> ArrangerClass;
    typedef FArrangerPeriodic<OctreeClass, ContainerClass, MoverClass> ArrangerClassPeriodic;

    //Pointer to the kernel to be executed
    InterKernel * kernel;
    MatrixKernelClass * matrix;
    //Link to the tree
    OctreeClass * octree;
    ArrangerClass * arranger;

public:
    /**
     * @brief Constructor : build the tree and the interpolation
     * kernel
     * @param TreeHeight Height of the tree
     * @param BoxWidth box Width
     * @param BoxCenter double[3] coordinate of the center of the
     * simulation box
     */
    FInterEngine(scalfmm_kernel_type KernelType) :
        kernel(nullptr), matrix(nullptr), octree(nullptr),arranger(nullptr){
        kernelType = KernelType;
    }

    void build_tree(int TreeHeight, double BoxWidth , double * BoxCenter,User_Scalfmm_Cell_Descriptor notUsedHere){
        octree = new OctreeClass(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint<FReal>(BoxCenter));
        this->matrix = new MatrixKernelClass();
        this->kernel = new InterKernel(TreeHeight,BoxWidth,FPoint<FReal>(BoxCenter),matrix);
    }


    //TODO free kernel too
    ~FInterEngine(){
        delete matrix;
        delete octree;
        delete kernel;
        if(arranger){
            delete arranger;
        }
    }

    //Inserting array of position
    void tree_insert_particles_xyz(int NbPositions, double * XYZ){
        for(int idPart = 0; idPart<NbPositions ; ++idPart){
            octree->insert(FPoint<FReal>(&XYZ[3*idPart]),idPart);
        }
        nbPart += NbPositions;
    }

    //Inserting arrayS of position
    void tree_insert_particles(int NbPositions, double * X, double * Y, double * Z){
        for(int idPart = 0; idPart<NbPositions ; ++idPart){
            octree->insert(FPoint<FReal>(X[idPart],Y[idPart],Z[idPart]),idPart);
        }
        nbPart += NbPositions;
    }

    //Set the physical values
    void set_physical_values(int nbPhysicalValues,double * physicalValues){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getPhysicalValues()[idxPart] = physicalValues[indexes[idxPart]];
                }
            });
    }

    //Set only a subpart of physical values
    //Algorithm : loop over each leaf, and then search in user array
    //if any index matches
    void set_physical_values_npart( int nbPhysicalValues, int* idxOfParticles, double * physicalValues){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbPhysicalValues && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            sources->getPhysicalValues()[idxPart] = physicalValues[indexes[idxPart]];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    //get back the physical values
    void get_physical_values( int nbPhysicalValues, double * physicalValues){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    physicalValues[indexes[idxPart]] = sources->getPhysicalValues()[idxPart];
                }
            });
    }

    //Same algorithm as in set_physical_values_npart
    void get_physical_values_npart( int nbPhysicalValues, int* idxOfParticles, double * physicalValues){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbPhysicalValues && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            physicalValues[indexes[idxPart]] = sources->getPhysicalValues()[idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    void get_forces_xyz( int nbParts, double * forcesToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    forcesToFill[indexes[idxPart]*3+0] = sources->getForcesX()[idxPart];
                    forcesToFill[indexes[idxPart]*3+1] = sources->getForcesY()[idxPart];
                    forcesToFill[indexes[idxPart]*3+2] = sources->getForcesZ()[idxPart];
                }
            });
    }

    void get_forces_xyz_npart(int nbParts, int* idxOfParticles , double * forcesToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbParts && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            forcesToFill[indexes[idxPart]*3+0] = sources->getForcesX()[idxPart];
                            forcesToFill[indexes[idxPart]*3+1] = sources->getForcesY()[idxPart];
                            forcesToFill[indexes[idxPart]*3+2] = sources->getForcesZ()[idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    void get_forces( int nbParts, double * fX, double* fY, double* fZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    fX[indexes[idxPart]] = sources->getForcesX()[idxPart];
                    fY[indexes[idxPart]] = sources->getForcesY()[idxPart];
                    fZ[indexes[idxPart]] = sources->getForcesZ()[idxPart];
                }
            });
    }

    void get_forces_npart(int nbParts, int* idxOfParticles ,double * fX, double* fY, double* fZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbParts && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            fX[indexes[idxPart]] = sources->getForcesX()[idxPart];
                            fY[indexes[idxPart]] = sources->getForcesY()[idxPart];
                            fZ[indexes[idxPart]] = sources->getForcesZ()[idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    //To set initial condition
    void set_forces_xyz( int nbParts, double * forcesToRead){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getForcesX()[idxPart] = forcesToRead[indexes[idxPart]*3+0];
                    sources->getForcesY()[idxPart] = forcesToRead[indexes[idxPart]*3+1];
                    sources->getForcesZ()[idxPart] = forcesToRead[indexes[idxPart]*3+2];
                }
            });
    }
    void set_forces_xyz_npart( int nbParts, int* idxOfParticles, double * forcesToRead){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbParts && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            sources->getForcesX()[idxPart] = forcesToRead[indexes[idxPart]*3+0];
                            sources->getForcesY()[idxPart] = forcesToRead[indexes[idxPart]*3+1];
                            sources->getForcesZ()[idxPart] = forcesToRead[indexes[idxPart]*3+2];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }
    void set_forces( int nbParts, double * fX, double* fY, double* fZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getForcesX()[idxPart] = fX[indexes[idxPart]];
                    sources->getForcesY()[idxPart] = fY[indexes[idxPart]];
                    sources->getForcesZ()[idxPart] = fZ[indexes[idxPart]];
                }
            });
    }
    void set_forces_npart( int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbParts && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            sources->getForcesX()[idxPart] = fX[indexes[idxPart]];
                            sources->getForcesY()[idxPart] = fY[indexes[idxPart]];
                            sources->getForcesZ()[idxPart] = fZ[indexes[idxPart]];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    //Set the potentials
    void set_potentials(int nbPotentials,double * potentialsToRead){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getPotentials()[idxPart] = potentialsToRead[indexes[idxPart]];
                }
            });
    }

    //Set only a subpart of potentials
    //Algorithm : loop over each leaf, and then search in user array
    //if any index matches
    void set_potentials_npart( int nbPotentials, int* idxOfParticles, double * potentialsToRead){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbPotentials && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            sources->getPotentials()[idxPart] = potentialsToRead[indexes[idxPart]];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    //get back the potentials
    void get_potentials( int nbPotentials, double * potentialsToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    potentialsToFill[indexes[idxPart]] = sources->getPotentials()[idxPart];
                }
            });
    }

    //Same algorithm as in set_potentials_npart
    void get_potentials_npart( int nbPotentials, int* idxOfParticles, double * potentialsToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < nbPotentials && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            potentialsToFill[indexes[idxPart]] = sources->getPotentials()[idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    void get_positions_xyz(int NbPositions, double * positionsToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    positionsToFill[indexes[idxPart]*3+0] = sources->getPositions()[0][idxPart];
                    positionsToFill[indexes[idxPart]*3+1] = sources->getPositions()[1][idxPart];
                    positionsToFill[indexes[idxPart]*3+2] = sources->getPositions()[2][idxPart];
                }
            });
    }

    void get_positions_xyz_npart(int NbPositions, int * idxOfParticles, double * positionsToFill){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < NbPositions && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            positionsToFill[indexes[idxPart]*3+0] =  sources->getPositions()[0][idxPart];
                            positionsToFill[indexes[idxPart]*3+1] =  sources->getPositions()[1][idxPart];
                            positionsToFill[indexes[idxPart]*3+2] =  sources->getPositions()[2][idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }

    void get_positions( int NbPositions, double * X, double * Y , double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    X[indexes[idxPart]] = sources->getPositions()[0][idxPart];
                    Y[indexes[idxPart]] = sources->getPositions()[1][idxPart];
                    Z[indexes[idxPart]] = sources->getPositions()[2][idxPart];
                }
            });
    }

    void get_positions_npart(int NbPositions, int * idxOfParticles,double * X, double * Y , double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    int iterPart = 0;
                    bool notFoundYet = true;
                    while(iterPart < NbPositions && notFoundYet){
                        if(indexes[idxPart] == idxOfParticles[iterPart]){
                            X[indexes[idxPart]] =  sources->getPositions()[0][idxPart];
                            Y[indexes[idxPart]] =  sources->getPositions()[1][idxPart];
                            Z[indexes[idxPart]] =  sources->getPositions()[2][idxPart];
                            notFoundYet = false;
                        }
                        else{
                            ++iterPart;
                        }
                    }
                }
            });
    }


    //Arranger parts : following function provide a way to move parts
    //inside the tree
    void add_to_positions_xyz(int NbPositions,double * updatedXYZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] += updatedXYZ[indexes[idxPart]*3+0];
                    sources->getWPositions()[1][idxPart] += updatedXYZ[indexes[idxPart]*3+1];
                    sources->getWPositions()[2][idxPart] += updatedXYZ[indexes[idxPart]*3+2];
                }
            });
    }

    void add_to_positions(int NbPositions,double * X, double * Y , double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] += X[indexes[idxPart]];
                    sources->getWPositions()[1][idxPart] += Y[indexes[idxPart]];
                    sources->getWPositions()[2][idxPart] += Z[indexes[idxPart]];
                }
            });
    }


    void set_positions_xyz(int NbPositions, double * updatedXYZ){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] = updatedXYZ[indexes[idxPart]*3+0];
                    sources->getWPositions()[1][idxPart] = updatedXYZ[indexes[idxPart]*3+1];
                    sources->getWPositions()[2][idxPart] = updatedXYZ[indexes[idxPart]*3+2];
                }
            });
    }

    void set_positions(int NbPositions, double * X, double * Y, double * Z){
        octree->forEachLeaf([&](LeafClass* leaf){
                ContainerClass * sources = leaf->getSrc();
                const FVector<int>& indexes = sources->getIndexes();
                int nbPartThere = sources->getNbParticles();
                for(int idxPart = 0 ; idxPart<nbPartThere ; ++idxPart){
                    sources->getWPositions()[0][idxPart] = X[indexes[idxPart]];
                    sources->getWPositions()[1][idxPart] = Y[indexes[idxPart]];
                    sources->getWPositions()[2][idxPart] = Z[indexes[idxPart]];
                }
            });
    }

    void update_tree(){
        if(arranger){
            arranger->rearrange();
        }
        else{
            if(Algorithm == 2){ //case in wich the periodic algorithm is used
                arranger = new ArrangerClassPeriodic(octree);
                arranger->rearrange();
            }
            else{
                arranger = new ArrangerClass(octree);
                arranger->rearrange();
            }
        }
    }


    void execute_fmm(){
        switch(Algorithm){
        case 0:
            {
                typedef FFmmAlgorithm<OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassSeq;
                AlgoClassSeq algoSeq(octree,kernel);
                algoSeq.execute();
                break;
            }
        case 1:
            {
                typedef FFmmAlgorithmThread<OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassThread;
                AlgoClassThread algoThread(octree,kernel);
                algoThread.execute();
                break;
            }
        case 2:
            {
                typedef FFmmAlgorithmPeriodic<FReal,OctreeClass,InterCell,ContainerClass,InterKernel,LeafClass> AlgoClassPeriodic;
                AlgoClassPeriodic algoPeriod(octree,2);
                algoPeriod.setKernel(kernel);
                algoPeriod.execute();
                break;
            }
        default :
            std::cout<< "No algorithm found (probably for strange reasons) : "<< Algorithm <<" exiting" << std::endl;
        }
    }

    void intern_dealloc_handle(Callback_free_cell unUsed){
        //this->~FInterEngine();
    }
};


#endif
