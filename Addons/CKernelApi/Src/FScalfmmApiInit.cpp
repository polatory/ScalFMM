/** It should be compiled with C export */
extern "C" {
#include "CScalfmmApi.h"
}

#include "FInterEngine.hpp"
#include "FUserKernelEngine.hpp"


extern "C" scalfmm_handle scalfmm_init(/*int TreeHeight,double BoxWidth,double* BoxCenter, */scalfmm_kernel_type KernelType,
                                       scalfmm_algorithm algo){
    ScalFmmCoreHandle<double> * handle = new ScalFmmCoreHandle<double>();
    typedef double FReal;

    if(algo == source_target){

        switch(KernelType){
        case 0:
            typedef FP2PParticleContainerIndexed<FReal>           ContainerClass;
            typedef FTypedLeaf<FReal,ContainerClass>                                         LeafClass;

            handle->engine = new FUserKernelEngine<FReal,LeafClass>(/*TreeHeight, BoxWidth, BoxCenter, */KernelType);
            break;

        case 1:
            //TODO typedefs
            typedef FP2PParticleContainerIndexed<FReal>                                 ContainerClass;
            typedef FTypedChebCell<FReal,7>                                                   ChebCell;
            typedef FTypedLeaf<FReal,ContainerClass>                                         LeafClass;

            typedef FInterpMatrixKernelR<FReal>                                        MatrixKernelClass;
            typedef FChebSymKernel<FReal,ChebCell,ContainerClass,MatrixKernelClass,7>        ChebKernel;

            handle->engine = new FInterEngine<FReal,ChebCell,ChebKernel,LeafClass>(/*TreeHeight,BoxWidth,BoxCenter, */KernelType,algo);
            break;
            // case 2:
            //     //TODO typedefs
            //     typedef FP2PParticleContainerIndexed<FReal>                                 ContainerClass;
            //     typedef FUnifCell<7>                                                         UnifCell;

            //     typedef FInterpMatrixKernelR<FReal>                                        MatrixKernelClass;
            //     typedef FUnifKernel<UnifCell,ContainerClass,MatrixKernelClass,7>           UnifKernel;

            //     handle->engine = new FInterEngine<UnifCell,UnifKernel>(/*TreeHeight,BoxWidth,BoxCenter, */KernelType);
            //     break;

        default:
            std::cout<< "Kernel type unsupported" << std::endl;
            exit(0);
            break;
        }
    }
    else{
        // if(algo == adaptiv){
        //     //Temporary
        //     handle->engine = new FAdaptEngine<FReal,4>(KernelType,algo);
        //}else{
            //No Source/Targets distinction
            switch(KernelType){
            case 0:
                typedef FP2PParticleContainerIndexed<FReal>                            ContainerClass;
                typedef FSimpleLeaf<FReal,ContainerClass>                                   LeafClass;

                handle->engine = new FUserKernelEngine<FReal,LeafClass>(/*TreeHeight, BoxWidth, BoxCenter, */KernelType);
                break;

            case 1:
                //TODO typedefs
                typedef FP2PParticleContainerIndexed<FReal>                                 ContainerClass;
                //typedef FChebCell<FReal,7>                                                   ChebCell;
                typedef FTypedChebCell<FReal,7>                                                   ChebCell;
                typedef FSimpleLeaf<FReal,ContainerClass>                                         LeafClass;

                typedef FInterpMatrixKernelR<FReal>                                        MatrixKernelClass;
                typedef FChebSymKernel<FReal,ChebCell,ContainerClass,MatrixKernelClass,7>        ChebKernel;

                handle->engine = new FInterEngine<FReal,ChebCell,ChebKernel,LeafClass>(/*TreeHeight,BoxWidth,BoxCenter, */KernelType,algo);
                break;
                // case 2:
                //     //TODO typedefs
                //     typedef FP2PParticleContainerIndexed<FReal>                                 ContainerClass;
                //     typedef FUnifCell<7>                                                         UnifCell;

                //     typedef FInterpMatrixKernelR<FReal>                                        MatrixKernelClass;
                //     typedef FUnifKernel<UnifCell,ContainerClass,MatrixKernelClass,7>           UnifKernel;

                //     handle->engine = new FInterEngine<UnifCell,UnifKernel>(/*TreeHeight,BoxWidth,BoxCenter, */KernelType);
                //     break;
            default:
                std::cout<< "Kernel type unsupported" << std::endl;
                exit(0);
                break;
            }
            // }
    }
   return handle;
}

extern "C" void scalfmm_dealloc_handle(scalfmm_handle handle, Callback_free_cell userDeallocator){
    ((ScalFmmCoreHandle<double> *) handle)->engine->intern_dealloc_handle(userDeallocator);
    delete ((ScalFmmCoreHandle<double> *) handle)->engine ;
    delete (ScalFmmCoreHandle<double> *) handle;
}

/**
 * This parts implements all the function defined in FChebInterface.h
 * using the Chebyshev classes
 */
#ifndef CHEBINTERFACE_HPP
#define CHEBINTERFACE_HPP

#warning "Compiling Cheb Interface"


extern "C" {
#include "Kernels/Chebyshev/FChebInterface.h"
}


#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Components/FSimpleLeaf.hpp"


typedef struct FChebCell_struct{
    //Store what's needed
    FChebCell<double,7> * cell;
}ChebCellStruct;


//How to create/destroy cells
extern "C" ChebCellStruct * ChebCellStruct_create(long long int inIndex,int * position){
    ChebCellStruct * newCell = new ChebCellStruct();
    newCell->cell = new FChebCell<double,7>();
    newCell->cell->setMortonIndex(inIndex);
    newCell->cell->setCoordinate(position[0],position[1],position[2]);
    return newCell;
}

extern "C" void ChebCellStruct_free(ChebCellStruct * inCell){
    if(inCell->cell) {
        delete inCell->cell;
    }
    delete inCell;
}



typedef struct FChebKernel_struct{
    //Be ready full duplication go there !!!
    FChebSymKernel<double,FChebCell<double,7>,FP2PParticleContainerIndexed<double>,FInterpMatrixKernelR<double>,7> ** kernel;
    FInterpMatrixKernelR<double> * matrix;
} ChebKernelStruct;

//Kernel functions
extern "C" ChebKernelStruct * ChebKernelStruct_create(int inTreeHeight,
                                                      double inBoxWidth,
                                                      double* inBoxCenter){
    //typedef to lighten the writing
    typedef FChebSymKernel<double,FChebCell<double,7>, FP2PParticleContainerIndexed<double>, FInterpMatrixKernelR<double>,7> KernelClass;
    ChebKernelStruct * newKernel = new ChebKernelStruct();
    newKernel->matrix= new FInterpMatrixKernelR<double>();
    int nb_threads = omp_get_max_threads();
    newKernel->kernel = new KernelClass*[nb_threads];
    newKernel->kernel[0]=new KernelClass(inTreeHeight,
                                         inBoxWidth,
                                         FPoint<double>(inBoxCenter[0], inBoxCenter[1], inBoxCenter[2]),
                                         newKernel->matrix);

    for(int idThreads=1 ; idThreads<nb_threads ; ++idThreads){
        newKernel->kernel[idThreads] = new KernelClass(*newKernel->kernel[0]);
    }
    return newKernel;
}

extern "C" void ChebKernelStruct_free(void *inKernel){
    delete reinterpret_cast<ChebKernelStruct *>(inKernel)->matrix;
    int nb_threads = omp_get_max_threads();
    for(int idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        delete reinterpret_cast<ChebKernelStruct *>(inKernel)->kernel[idThreads];
    }
    delete [] reinterpret_cast<ChebKernelStruct *>(inKernel)->kernel;
    delete reinterpret_cast<ChebKernelStruct *>(inKernel);
}


extern "C" void ChebKernel_P2M(void * leafCell, FSize nbParticles, const FSize *particleIndexes, void * inKernel){
    //make temporary array of parts
    FP2PParticleContainerIndexed<double>* tempContainer = new FP2PParticleContainerIndexed<double>();
    tempContainer->reserve(nbParticles);
    FPoint<double> pos;
    for(int i=0 ; i<nbParticles ; ++i){
        pos = FPoint<double>(reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3  ],
                             reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3+1],
                             reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3+2]);
        double Phi = reinterpret_cast<UserData *>(inKernel)->myPhyValues[particleIndexes[i]];

        tempContainer->push(pos,particleIndexes[i],Phi);
    }
    //get the real cell struct
    ChebCellStruct * realCellStruct = reinterpret_cast<ChebCellStruct *>(leafCell);
    FChebCell<double,7> * realCell = realCellStruct->cell;

    //Identify thread number
    int id_thread = omp_get_thread_num();

    //get the real chebyshev struct
    UserData * userDataKernel = reinterpret_cast<UserData *>(inKernel);
    ChebKernelStruct * realKernel = userDataKernel->kernelStruct;

    realKernel->kernel[id_thread]->P2M(realCell, tempContainer);
    delete tempContainer;
}

extern "C" void  ChebKernel_M2M(int level, void* parentCell, int childPosition, void *childCell, void *inKernel){
    //Get our structures
    ChebCellStruct * parentCellStruct = reinterpret_cast<ChebCellStruct *>(parentCell);
    ChebCellStruct * childCellStruct = reinterpret_cast<ChebCellStruct *>(childCell);
    //get real cheb cell
    FChebCell<double,7>* parentChebCell = parentCellStruct->cell;
    FChebCell<double,7>* childChebCell = childCellStruct->cell;

    //Identify thread number
    int id_thread = omp_get_thread_num();

    //Get the kernel
    ChebKernelStruct * inKernelStruct = reinterpret_cast<UserData*>(inKernel)->kernelStruct;
    inKernelStruct->kernel[id_thread]->getPtrToInterpolator()->applyM2M(childPosition,
                                                                        childChebCell->getMultipole(0),
                                                                        parentChebCell->getMultipole(0));
}

extern "C" void ChebKernel_M2L(int level, void* targetCell,  void* sourceCell[343], void* inKernel){
    //Get our structures
    ChebCellStruct * targetCellStruct = reinterpret_cast<ChebCellStruct *>(targetCell);
    //get real cheb cell
    FChebCell<double,7>* const targetChebCell = targetCellStruct->cell;

    //copy to an array of FChebCell
    const FChebCell<double,7>* arrayOfChebCell[343];
    for(int i=0; i<343 ; ++i){
        if(sourceCell[i] != nullptr){
            arrayOfChebCell[i] = reinterpret_cast<ChebCellStruct*>(sourceCell[i])->cell;
        }
        else{
            arrayOfChebCell[i] = nullptr;
        }
    }

    //Identify thread number
    int id_thread = omp_get_thread_num();

    //Get the kernel
    ChebKernelStruct * inKernelStruct = reinterpret_cast<UserData*>(inKernel)->kernelStruct;
    inKernelStruct->kernel[id_thread]->M2L(targetChebCell,arrayOfChebCell,0,level);
}

extern "C" void ChebKernel_L2L(int level, void* parentCell, int childPosition, void* childCell, void* inKernel){
    //Get our structures
    ChebCellStruct * parentCellStruct = reinterpret_cast<ChebCellStruct *>(parentCell);
    ChebCellStruct * childCellStruct = reinterpret_cast<ChebCellStruct *>(childCell);
    //get real cheb cell
    FChebCell<double,7>* parentChebCell = parentCellStruct->cell;
    FChebCell<double,7>* childChebCell = childCellStruct->cell;

    //Identify thread number
    int id_thread = omp_get_thread_num();

    //Get the kernel
    ChebKernelStruct * inKernelStruct = reinterpret_cast<UserData*>(inKernel)->kernelStruct;
    inKernelStruct->kernel[id_thread]->getPtrToInterpolator()->applyL2L(childPosition,
                                                                        parentChebCell->getLocal(0),
                                                                        childChebCell->getLocal(0));
}

extern "C" void ChebKernel_L2P(void* leafCell, FSize nbParticles, const FSize* particleIndexes, void* inKernel){
    //Create temporary FSimpleLeaf
    FP2PParticleContainerIndexed<double>* tempContainer = new FP2PParticleContainerIndexed<double>();
    tempContainer->reserve(nbParticles);
    FPoint<double> pos;
    for(int i=0 ; i<nbParticles ; ++i){
        pos = FPoint<double>(reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3  ],
                             reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3+1],
                             reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3+2]);
        double Phi = reinterpret_cast<UserData *>(inKernel)->myPhyValues[particleIndexes[i]];
        tempContainer->push(pos,particleIndexes[i],Phi);
    }
    //Get our structures
    ChebCellStruct * leafCellStruct = reinterpret_cast<ChebCellStruct *>(leafCell);
    //get real cheb cell
    FChebCell<double,7>* leafChebCell = leafCellStruct->cell;

    //Identify thread number
    int id_thread = omp_get_thread_num();

    //Get the kernel
    ChebKernelStruct * inKernelStruct = reinterpret_cast<UserData*>(inKernel)->kernelStruct;

    inKernelStruct->kernel[id_thread]->L2P(leafChebCell,tempContainer);

    //Then retrieve the results
    double * forcesToFill     = reinterpret_cast<UserData *>(inKernel)->forcesComputed[id_thread];
    double * potentialsToFill = reinterpret_cast<UserData *>(inKernel)->potentials[id_thread];

    const FVector<FSize>& indexes = tempContainer->getIndexes();
    for(FSize idxPart = 0 ; idxPart<nbParticles ; ++idxPart){
        forcesToFill[indexes[idxPart]*3+0] += tempContainer->getForcesX()[idxPart];
        forcesToFill[indexes[idxPart]*3+1] += tempContainer->getForcesY()[idxPart];
        forcesToFill[indexes[idxPart]*3+2] += tempContainer->getForcesZ()[idxPart];
        potentialsToFill[indexes[idxPart]] += tempContainer->getPotentials()[idxPart];
    }

    delete tempContainer;
    tempContainer=nullptr;
}


void ChebKernel_P2P(FSize nbParticles, const FSize* particleIndexes,
                    const FSize * sourceParticleIndexes[27],FSize sourceNbPart[27],void* inKernel){

    //Create temporary FSimpleLeaf for target
    FP2PParticleContainerIndexed<double>* tempContTarget = new FP2PParticleContainerIndexed<double>();
    tempContTarget->reserve(nbParticles);
    for(int i=0 ; i<nbParticles ; ++i){
        FPoint<double> pos = FPoint<double>(reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3  ],
                                            reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3+1],
                                            reinterpret_cast<UserData *>(inKernel)->insertedPositions[particleIndexes[i]*3+2]);
        double Phi = reinterpret_cast<UserData *>(inKernel)->myPhyValues[particleIndexes[i]];
        tempContTarget->push(pos,particleIndexes[i],Phi);
    }

    //Create 27 FSimpleLeaf for 27 sources
    FP2PParticleContainerIndexed<double>* tempContSources[27];
    for(int idSource=0; idSource<27 ; ++idSource){
        if(sourceNbPart[idSource] != 0){
            //Create container
            tempContSources[idSource] = new FP2PParticleContainerIndexed<double>();
            //Allocate memory
            tempContSources[idSource]->reserve(sourceNbPart[idSource]);
            //Store a ptr to the indices of that source leaf
            const FSize * indSource = sourceParticleIndexes[idSource];
            //Then, for each part in this source
            for(int i=0 ; i<sourceNbPart[idSource] ; ++i){
                FPoint<double> pos = FPoint<double>(reinterpret_cast<UserData *>(inKernel)->insertedPositions[indSource[i]*3  ],
                                                    reinterpret_cast<UserData *>(inKernel)->insertedPositions[indSource[i]*3+1],
                                                    reinterpret_cast<UserData *>(inKernel)->insertedPositions[indSource[i]*3+2]);
                double Phi = reinterpret_cast<UserData *>(inKernel)->myPhyValues[indSource[i]];
                tempContSources[idSource]->push(pos,indSource[i],Phi);
            }
        }
        else{
            tempContSources[idSource] = nullptr;
        }
    }
    //Everything is fine, now, call Chebyshev P2P

    //Identify thread number
    int id_thread = omp_get_thread_num();

    //Get the kernel
    ChebKernelStruct * inKernelStruct = reinterpret_cast<UserData*>(inKernel)->kernelStruct;

    //Empty tree coordinate
    int coord[3] = {0,0,0};

    inKernelStruct->kernel[id_thread]->P2P(FTreeCoordinate(coord),tempContTarget,nullptr,tempContSources,0);

    //get back forces & potentials
    double * forcesToFill = reinterpret_cast<UserData *>(inKernel)->forcesComputed[id_thread];
    double * potentialsToFill = reinterpret_cast<UserData *>(inKernel)->potentials[id_thread];

    const FVector<FSize>& indexes = tempContTarget->getIndexes();
    for(FSize idxPart = 0 ; idxPart<nbParticles ; ++idxPart){
        forcesToFill[indexes[idxPart]*3+0] += tempContTarget->getForcesX()[idxPart];
        forcesToFill[indexes[idxPart]*3+1] += tempContTarget->getForcesY()[idxPart];
        forcesToFill[indexes[idxPart]*3+2] += tempContTarget->getForcesZ()[idxPart];
        potentialsToFill[indexes[idxPart]] += tempContTarget->getPotentials()[idxPart];
    }

    //Note that sources are also modified.
    //get back sources forces
    for(int idSource = 0 ; idSource < 27 ; ++idSource){
        const FVector<FSize>& indexesSource = tempContSources[idSource]->getIndexes();
        const FSize nbPartInSource = sourceNbPart[idSource];
        for(int idxSourcePart = 0; idxSourcePart < nbPartInSource ; ++idxSourcePart){
            forcesToFill[indexesSource[idxSourcePart]*3+0] += tempContSources[idSource]->getForcesX()[idxSourcePart];
            forcesToFill[indexesSource[idxSourcePart]*3+1] += tempContSources[idSource]->getForcesY()[idxSourcePart];
            forcesToFill[indexesSource[idxSourcePart]*3+2] += tempContSources[idSource]->getForcesZ()[idxSourcePart];
            potentialsToFill[indexesSource[idxSourcePart]] += tempContSources[idSource]->getPotentials()[idxSourcePart];
        }
    }

    //Release memory
    for(int idSource=0; idSource<27 ; ++idSource){
        if(tempContSources[idSource]) delete tempContSources[idSource];
    }
    delete tempContTarget;
}

void ChebCell_reset(int level, long long morton_index, int* tree_position, double* spatial_position, void * userCell,void * inKernel){
    ChebCellStruct *  cellStruct = reinterpret_cast<ChebCellStruct *>(userCell);
    FChebCell<double,7>* chebCell = cellStruct->cell;
    chebCell->resetToInitialState();
}

#endif
