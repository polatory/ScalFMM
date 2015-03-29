
// @SCALFMM_PRIVATE
#ifndef FSTARPUCPUWRAPPER_HPP
#define FSTARPUCPUWRAPPER_HPP


#include "../../Utils/FGlobal.hpp"
#include "../../Core/FCoreCommon.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FLog.hpp"
#include "../../Utils/FTic.hpp"
#include "../../Utils/FAssert.hpp"

#include "../Core/FOutOfBlockInteraction.hpp"

#ifdef SCALFMM_USE_MPI
#include "../../Utils/FMpi.hpp"
#endif

#include <vector>
#include <memory>

#include <omp.h>

//extern "C"{
#include <starpu.h>
//}

#ifdef STARPU_USE_MPI
//extern "C"{
#include <starpu_mpi.h>
//}
#endif

#include "FStarPUUtils.hpp"

template <class CellContainerClass, class CellClass, class KernelClass,
          class ParticleGroupClass, class ParticleContainerClass>
class FStarPUCpuWrapper {
protected:
    typedef FStarPUCpuWrapper<CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass> ThisClass;

    template <class OtherBlockClass>
    struct BlockInteractions{
        OtherBlockClass* otherBlock;
        int otherBlockId;
        std::vector<OutOfBlockInteraction> interactions;
    };

    const int treeHeight;
    KernelClass* kernels[STARPU_MAXCPUS];        //< The kernels

public:
    FStarPUCpuWrapper(const int inTreeHeight): treeHeight(inTreeHeight){
        memset(kernels, 0, sizeof(KernelClass*)*STARPU_MAXCPUS);
    }

    void initKernel(const int workerId, KernelClass* originalKernel){
        FAssertLF(kernels[workerId] == nullptr);
        kernels[workerId] = new KernelClass(*originalKernel);
    }

    void releaseKernel(const int workerId){
        delete kernels[workerId];
        kernels[workerId] = nullptr;
    }

    ~FStarPUCpuWrapper(){
        for(int idxKernel = 0 ; idxKernel < STARPU_MAXCPUS ; ++idxKernel ){
            FAssertLF(kernels[idxKernel] == nullptr);
        }
    }

    static void bottomPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass leafCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                            STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                            (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                            nullptr);
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                            STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                            nullptr);

        FStarPUPtrInterface* worker = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &intervalSize);
        worker->get<ThisClass>(FSTARPU_CPU_IDX)->bottomPassPerform(&leafCells, &containers);
    }

    void bottomPassPerform(CellContainerClass* leafCells, ParticleGroupClass* containers){
        const MortonIndex blockStartIdx = leafCells->getStartingIndex();
        const MortonIndex blockEndIdx = leafCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){            
            if(leafCells->exists(mindex)){
                CellClass cell = leafCells->getUpCell(mindex);
                FAssertLF(cell.getMortonIndex() == mindex);
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                FAssertLF(particles.isAttachedToSomething());
                kernel->P2M(&cell, &particles);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Upward Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void upwardPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                        nullptr);

        FStarPUPtrInterface* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel, &intervalSize);

        CellContainerClass* subCellGroups[9];
        memset(subCellGroups, 0, 9*sizeof(CellContainerClass*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroups[idxSubGroup] = new CellContainerClass(
                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+2]),
                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[(idxSubGroup*2)+2]),
                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+3]),
                        nullptr);
        }

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->upwardPassPerform(&currentCells, subCellGroups, nbSubCellGroups, idxLevel);

        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            delete subCellGroups[idxSubGroup];
        }
    }

    void upwardPassPerform(CellContainerClass*const currentCells,
                           CellContainerClass* subCellGroups[9],
                            const int nbSubCellGroups, const int idxLevel){
        FAssertLF(nbSubCellGroups != 0);
        const MortonIndex blockStartIdx = FMath::Max(currentCells->getStartingIndex(),
                                              subCellGroups[0]->getStartingIndex()>>3);
        const MortonIndex blockEndIdx   = FMath::Min(currentCells->getEndingIndex(),
                                              ((subCellGroups[nbSubCellGroups-1]->getEndingIndex()-1)>>3)+1);
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        int idxSubCellGroup = 0;

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
            if(currentCells->exists(mindex)){
                CellClass cell = currentCells->getUpCell(mindex);
                FAssertLF(cell.getMortonIndex() == mindex);
                CellClass childData[8];
                CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if( subCellGroups[idxSubCellGroup]->getEndingIndex() <= ((mindex<<3)+idxChild) ){
                        idxSubCellGroup += 1;
                    }
                    if( idxSubCellGroup == nbSubCellGroups ){
                        break;
                    }
                    if(subCellGroups[idxSubCellGroup]->exists((mindex<<3)+idxChild)){
                        childData[idxChild] = subCellGroups[idxSubCellGroup]->getUpCell((mindex<<3)+idxChild);
                        child[idxChild] = &childData[idxChild];
                        FAssertLF(child[idxChild] == nullptr || child[idxChild]->getMortonIndex() == ((mindex<<3)+idxChild));
                    }
                }

                kernel->M2M(&cell, child, idxLevel);
            }
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass Mpi
    /////////////////////////////////////////////////////////////////////////////////////
#ifdef STARPU_USE_MPI
    static void transferInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                        nullptr,
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]));
        CellContainerClass externalCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]),
                                        nullptr);

        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions, &intervalSize);

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->transferInoutPassPerformMpi(&currentCells, &externalCells, idxLevel, outsideInteractions);
    }


    void transferInoutPassPerformMpi(CellContainerClass*const currentCells,
                                  CellContainerClass*const cellsOther,
                                  const int idxLevel,
                                  const std::vector<OutOfBlockInteraction>* outsideInteractions){
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
            if(cellsOther->exists((*outsideInteractions)[outInterIdx].outIndex)){
                CellClass interCell = cellsOther->getDownCell((*outsideInteractions)[outInterIdx].outIndex);
                FAssertLF(interCell.getMortonIndex() == (*outsideInteractions)[outInterIdx].outIndex);
                CellClass cell = currentCells->getUpCell((*outsideInteractions)[outInterIdx].insideIndex);
                FAssertLF(cell.getMortonIndex() == (*outsideInteractions)[outInterIdx].insideIndex);

                const CellClass* interactions[343];
                memset(interactions, 0, 343*sizeof(CellClass*));
                interactions[(*outsideInteractions)[outInterIdx].outPosition] = &interCell;
                const int counter = 1;
                kernel->M2L( &cell , interactions, counter, idxLevel);
            }
        }
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void transferInPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]));

        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &intervalSize);

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->transferInPassPerform(&currentCells, idxLevel);
    }

    void transferInPassPerform(CellContainerClass*const currentCells, const int idxLevel){
        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
        const MortonIndex blockEndIdx = currentCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            if(currentCells->exists(mindex)){
                CellClass cell = currentCells->getCompleteCell(mindex);
                FAssertLF(cell.getMortonIndex() == mindex);
                MortonIndex interactionsIndexes[189];
                int interactionsPosition[189];
                const FTreeCoordinate coord(cell.getCoordinate());
                int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                CellClass interactionsData[343];
                const CellClass* interactions[343];
                memset(interactions, 0, 343*sizeof(CellClass*));
                int counterExistingCell = 0;

                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                        if(currentCells->exists(interactionsIndexes[idxInter])){
                            interactionsData[interactionsPosition[idxInter]] = currentCells->getCompleteCell(interactionsIndexes[idxInter]);
                            FAssertLF(interactionsData[interactionsPosition[idxInter]].getMortonIndex() == interactionsIndexes[idxInter]);
                            FAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                            interactions[interactionsPosition[idxInter]] = &interactionsData[interactionsPosition[idxInter]];
                            counterExistingCell += 1;
                        }
                    }
                }

                kernel->M2L( &cell , interactions, counterExistingCell, idxLevel);
            }
        }
    }

    static void transferInoutPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]));
        CellContainerClass externalCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[3]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[4]),
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[5]));

        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions, &intervalSize);

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->transferInoutPassPerform(&currentCells, &externalCells, idxLevel, outsideInteractions);
    }


    void transferInoutPassPerform(CellContainerClass*const currentCells,
                                  CellContainerClass*const cellsOther,
                                  const int idxLevel,
                                  const std::vector<OutOfBlockInteraction>* outsideInteractions){
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
            if(cellsOther->exists((*outsideInteractions)[outInterIdx].outIndex)){
                CellClass interCell = cellsOther->getCompleteCell((*outsideInteractions)[outInterIdx].outIndex);
                FAssertLF(interCell.getMortonIndex() == (*outsideInteractions)[outInterIdx].outIndex);
                CellClass cell = currentCells->getCompleteCell((*outsideInteractions)[outInterIdx].insideIndex);
                FAssertLF(cell.getMortonIndex() == (*outsideInteractions)[outInterIdx].insideIndex);

                const CellClass* interactions[343];
                memset(interactions, 0, 343*sizeof(CellClass*));
                interactions[(*outsideInteractions)[outInterIdx].outPosition] = &interCell;
                const int counter = 1;
                kernel->M2L( &cell , interactions, counter, idxLevel);

                interactions[(*outsideInteractions)[outInterIdx].outPosition] = nullptr;
                interactions[getOppositeInterIndex((*outsideInteractions)[outInterIdx].outPosition)] = &cell;
                kernel->M2L( &interCell , interactions, counter, idxLevel);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Downard Pass
    /////////////////////////////////////////////////////////////////////////////////////
    static void downardPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                        nullptr,
                                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel, &intervalSize);

        CellContainerClass* subCellGroups[9];
        memset(subCellGroups, 0, 9*sizeof(CellContainerClass*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroups[idxSubGroup] = new CellContainerClass(
                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+2]),
                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[(idxSubGroup*2)+2]),
                        nullptr,
                        (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[(idxSubGroup*2)+3]));
        }

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->downardPassPerform(&currentCells, subCellGroups, nbSubCellGroups, idxLevel);

        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            delete subCellGroups[idxSubGroup];
        }
    }

    void downardPassPerform(CellContainerClass*const currentCells,
                            CellContainerClass* subCellGroups[9],
                             const int nbSubCellGroups, const int idxLevel){
        FAssertLF(nbSubCellGroups != 0);
        const MortonIndex blockStartIdx = FMath::Max(currentCells->getStartingIndex(),
                                              subCellGroups[0]->getStartingIndex()>>3);
        const MortonIndex blockEndIdx   = FMath::Min(currentCells->getEndingIndex(),
                                              ((subCellGroups[nbSubCellGroups-1]->getEndingIndex()-1)>>3)+1);
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        int idxSubCellGroup = 0;

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
            if(currentCells->exists(mindex)){
                CellClass cell = currentCells->getDownCell(mindex);
                FAssertLF(cell.getMortonIndex() == mindex);
                CellClass childData[8];
                CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if( subCellGroups[idxSubCellGroup]->getEndingIndex() <= ((mindex<<3)+idxChild) ){
                        idxSubCellGroup += 1;
                    }
                    if( idxSubCellGroup == nbSubCellGroups ){
                        break;
                    }
                    if(subCellGroups[idxSubCellGroup]->exists((mindex<<3)+idxChild)){
                        childData[idxChild] = subCellGroups[idxSubCellGroup]->getDownCell((mindex<<3)+idxChild);
                        child[idxChild] = &childData[idxChild];
                        FAssertLF(child[idxChild] == nullptr || child[idxChild]->getMortonIndex() == ((mindex<<3)+idxChild));
                    }
                }

                kernel->L2L(&cell, child, idxLevel);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass MPI
    /////////////////////////////////////////////////////////////////////////////////////

#ifdef STARPU_USE_MPI
    static void directInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                      (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]));
        ParticleGroupClass externalContainers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                                      nullptr);

        FStarPUPtrInterface* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions, &intervalSize);

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->directInoutPassPerform(&containers, &externalContainers, outsideInteractions);
    }

    void directInoutPassPerformMpi(ParticleGroupClass* containers, ParticleGroupClass* containersOther,
                                const std::vector<OutOfBlockInteraction>* outsideInteractions){
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
            ParticleContainerClass interParticles = containersOther->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].outIndex);
            if(interParticles.isAttachedToSomething()){
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].insideIndex);
                FAssertLF(particles.isAttachedToSomething());
                ParticleContainerClass* interactions[27];
                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));
                interactions[(*outsideInteractions)[outInterIdx].outPosition] = &interParticles;
                const int counter = 1;
                kernel->P2PRemote( FTreeCoordinate((*outsideInteractions)[outInterIdx].insideIndex, treeHeight-1), &particles, &particles , interactions, counter);
            }
        }
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void directInPassCallback(void *buffers[], void *cl_arg){
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                      (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &intervalSize);
        worker->get<ThisClass>(FSTARPU_CPU_IDX)->directInPassPerform(&containers);
    }

    void directInPassPerform(ParticleGroupClass* containers){
        const MortonIndex blockStartIdx = containers->getStartingIndex();
        const MortonIndex blockEndIdx = containers->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
            if(particles.isAttachedToSomething()){
                MortonIndex interactionsIndexes[26];
                int interactionsPosition[26];
                FTreeCoordinate coord(mindex, treeHeight-1);
                int counter = coord.getNeighborsIndexes(treeHeight,interactionsIndexes,interactionsPosition);

                ParticleContainerClass interactionsObjects[27];
                ParticleContainerClass* interactions[27];
                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));
                int counterExistingCell = 0;

                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                        interactionsObjects[counterExistingCell] = containers->template getLeaf<ParticleContainerClass>(interactionsIndexes[idxInter]);
                        if(interactionsObjects[counterExistingCell].isAttachedToSomething()){
                            FAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                            interactions[interactionsPosition[idxInter]] = &interactionsObjects[counterExistingCell];
                            counterExistingCell += 1;
                        }
                    }
                }

                kernel->P2P( coord, &particles, &particles , interactions, counterExistingCell);
            }
        }
    }

    static void directInoutPassCallback(void *buffers[], void *cl_arg){
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                      (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]));
        ParticleGroupClass externalContainers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                                      (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]));

        FStarPUPtrInterface* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions, &intervalSize);

        worker->get<ThisClass>(FSTARPU_CPU_IDX)->directInoutPassPerform(&containers, &externalContainers, outsideInteractions);
    }

    void directInoutPassPerform(ParticleGroupClass* containers, ParticleGroupClass* containersOther,
                                const std::vector<OutOfBlockInteraction>* outsideInteractions){
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
            ParticleContainerClass interParticles = containersOther->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].outIndex);
            if(interParticles.isAttachedToSomething()){
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].insideIndex);
                FAssertLF(particles.isAttachedToSomething());
                ParticleContainerClass* interactions[27];
                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));
                interactions[(*outsideInteractions)[outInterIdx].outPosition] = &interParticles;
                const int counter = 1;
                kernel->P2PRemote( FTreeCoordinate((*outsideInteractions)[outInterIdx].insideIndex, treeHeight-1), &particles, &particles , interactions, counter);

                interactions[(*outsideInteractions)[outInterIdx].outPosition] = nullptr;
                interactions[getOppositeNeighIndex((*outsideInteractions)[outInterIdx].outPosition)] = &particles;
                kernel->P2PRemote( FTreeCoordinate((*outsideInteractions)[outInterIdx].outIndex, treeHeight-1), &interParticles, &interParticles , interactions, counter);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Merge Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void mergePassCallback(void *buffers[], void *cl_arg){
        CellContainerClass leafCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                     STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                                     nullptr,
                                     (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]));
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[2]),
                                     STARPU_VARIABLE_GET_ELEMSIZE(buffers[2]),
                                     (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[3]));

        FStarPUPtrInterface* worker = nullptr;
        int intervalSize;
        starpu_codelet_unpack_args(cl_arg, &worker, &intervalSize);
        worker->get<ThisClass>(FSTARPU_CPU_IDX)->mergePassPerform(&leafCells, &containers);
    }

    void mergePassPerform(CellContainerClass* leafCells, ParticleGroupClass* containers){
        const MortonIndex blockStartIdx = leafCells->getStartingIndex();
        const MortonIndex blockEndIdx = leafCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            if(leafCells->exists(mindex)){
                CellClass cell = leafCells->getDownCell(mindex);
                FAssertLF(cell.getMortonIndex() == mindex);
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                FAssertLF(particles.isAttachedToSomething());
                kernel->L2P(&cell, &particles);
            }
        }
    }

    static int getOppositeNeighIndex(const int index) {
        // ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1)
        return 27-index-1;
    }

    static int getOppositeInterIndex(const int index) {
        // ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3
        return 343-index-1;
    }
};

#endif // FSTARPUCPUWRAPPER_HPP

