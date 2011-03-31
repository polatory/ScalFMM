#ifndef FFMMALGORITHMTHREADED_HPP
#define FFMMALGORITHMTHREADED_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FOpenMPThread.hpp"

#include "../Containers/FOctree.hpp"

#include <omp.h>


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFMMAlgorithmThreaded
* @brief
* Please read the license
* This is a parallel FMM Algorithm
* Most of your code is unchanged excepted that your kerneks must have a copy operator.
*
* The parallel algorithm is simple, each thread is taking a value from the iterator (protected by a mutex)
*/
template<template< class ParticuleClass, class CellClass, int OctreeHeight> class KernelClass, class ParticuleClass, class CellClass, int OctreeHeight, int SubtreeHeight>
class FFMMAlgorithmThreaded : protected FAssertable{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator FOctreeIterator;

    static const int NbThreads = 4;                                             //< Number of threads (currently a static number)

    KernelClass<ParticuleClass, CellClass, OctreeHeight>* kernels[NbThreads];   //< The kernels (one by thread)
    Octree* const tree;                                                         //< The octree to work on

    FDEBUG(FTic counter);                                                  //< In case of debug count the time

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree
      * @param inKernels the kernels
      * an assert is launched if one of the arguments is null
      */
    FFMMAlgorithmThreaded(Octree* const inTree,
                  KernelClass<ParticuleClass, CellClass, OctreeHeight>* const inKernels)
                    : tree(inTree) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        for(int idxThread = 0 ; idxThread < NbThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticuleClass, CellClass, OctreeHeight>(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFMMAlgorithmThreaded\n" );
    }

    /** Default destructor */
    virtual ~FFMMAlgorithmThreaded(){
        for(int idxThread = 0 ; idxThread < NbThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
    }

    /** To execute the fmm algorithm
      * Call this function to run the complete algo
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        for(int idxThread = 0 ; idxThread < NbThreads ; ++idxThread){
            this->kernels[idxThread]->init();
        }

        bottomPass();
        upwardPass();

        downardPass();

        directPass();
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( counter.tic() );

        FOctreeIterator octreeIterator(tree);
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();

        omp_lock_t mutex;
        omp_init_lock(&mutex);
        bool stop = false;
        #pragma omp parallel shared(octreeIterator,mutex,stop) num_threads(NbThreads)
        {
            const int threadId = omp_get_thread_num();
            omp_set_lock(&mutex);
            while(!stop){
                CellClass*const cell = octreeIterator.getCurrentCell();
                const FList<ParticuleClass*>* const particules = octreeIterator.getCurrentList();
                if(!octreeIterator.moveRight()) stop = true;
                omp_unset_lock(&mutex);

                // We need the current cell that represent the leaf
                // and the list of particules
                kernels[threadId]->P2M( cell , particules);

                omp_set_lock(&mutex);
            }
            omp_unset_lock(&mutex);
        }
        omp_destroy_lock(&mutex);

        FDEBUG( counter.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counter.elapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG( counter.tic() );

        FOctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        omp_lock_t mutex;
        omp_init_lock(&mutex);
        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            bool stop = false;
            #pragma omp parallel shared(octreeIterator,mutex,stop) num_threads(NbThreads)
            {
                const int threadId = omp_get_thread_num();
                omp_set_lock(&mutex);
                // for each cells
                while(!stop){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    CellClass*const cell = octreeIterator.getCurrentCell();
                    const CellClass*const*const child = octreeIterator.getCurrentChild();
                    if(!octreeIterator.moveRight()) stop = true;
                    omp_unset_lock(&mutex);

                    kernels[threadId]->M2M( cell , child, idxLevel);

                    omp_set_lock(&mutex);
                }
                omp_unset_lock(&mutex);
            }
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
        }
        omp_destroy_lock(&mutex);

        FDEBUG( counter.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counter.elapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG( counter.tic() );

        { // first M2L
            FOctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            omp_lock_t mutex;
            omp_init_lock(&mutex);

            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                bool stop = false;
                #pragma omp parallel shared(octreeIterator,mutex,idxLevel,stop) num_threads(NbThreads)
                {
                    const int threadId = omp_get_thread_num();
                    CellClass* neighbors[208];
                    omp_set_lock(&mutex);
                    // for each cells
                    while(!stop){
                        CellClass* const cell = octreeIterator.getCurrentCell();
                        const MortonIndex cellIndex = octreeIterator.getCurrentGlobalIndex();
                        if(!octreeIterator.moveRight()) stop = true;
                        omp_unset_lock(&mutex);

                        const int counter = tree->getDistantNeighbors(neighbors, cellIndex,idxLevel);
                        if(counter) kernels[threadId]->M2L( cell, neighbors, counter, idxLevel);

                        omp_set_lock(&mutex);
                    }
                    omp_unset_lock(&mutex);
                }
                octreeIterator.gotoLeft();
                octreeIterator.moveDown();
            }
            omp_destroy_lock(&mutex);
        }
        FDEBUG( counter.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counter.elapsed() << "s)\n" );

        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG( counter.tic() );
        { // second L2L
            FOctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();
            const int heightMinusOne = OctreeHeight - 1;

            omp_lock_t mutex;
            omp_init_lock(&mutex);
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                bool stop = false;
                #pragma omp parallel shared(octreeIterator,mutex,idxLevel,stop) num_threads(NbThreads)
                {
                    const int threadId = omp_get_thread_num();
                    omp_set_lock(&mutex);
                    // for each cells
                    while(!stop){
                        const CellClass * const cell = octreeIterator.getCurrentCell();
                        CellClass ** const child = octreeIterator.getCurrentChild();
                        if(!octreeIterator.moveRight()) stop = true;
                        omp_unset_lock(&mutex);

                        kernels[threadId]->L2L( cell, child, idxLevel);

                        omp_set_lock(&mutex);
                    }
                    omp_unset_lock(&mutex);
                }

                octreeIterator.gotoLeft();
                octreeIterator.moveDown();
            }
            omp_destroy_lock(&mutex);
        }

        FDEBUG( counter.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counter.elapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2P */
    void directPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG( counter.tic() );

        FOctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();

        omp_lock_t mutex;
        omp_init_lock(&mutex);
        bool stop = false;
        #pragma omp parallel shared(octreeIterator,mutex,stop) num_threads(NbThreads)
        {
            const int threadId = omp_get_thread_num();
            // There is a maximum of 26 neighbors
            FList<ParticuleClass*>* neighbors[26];
            const int heightMinusOne = OctreeHeight - 1;

            omp_set_lock(&mutex);
            // for each leafs
            while(!stop){
                const CellClass * const cell = octreeIterator.getCurrentCell();
                FList<ParticuleClass*>* const particules = octreeIterator.getCurrentList();
                const MortonIndex cellIndex = octreeIterator.getCurrentGlobalIndex();
                if(!octreeIterator.moveRight()) stop = true;
                omp_unset_lock(&mutex);

                kernels[threadId]->L2P(cell, particules);
                // need the current particules and neighbors particules
                const int counter = tree->getLeafsNeighbors(neighbors, cellIndex,heightMinusOne);
                kernels[threadId]->P2P( particules , neighbors, counter);

                omp_set_lock(&mutex);
            }
            omp_unset_lock(&mutex);
        }
        omp_destroy_lock(&mutex);

        FDEBUG( counter.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counter.elapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};


#endif //FFMMALGORITHMTHREADED_HPP

// [--LICENSE--]
