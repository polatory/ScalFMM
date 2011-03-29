#ifndef FFMMALGORITHMTHREADEDTHREADED_HPP
#define FFMMALGORITHMTHREADEDTHREADED_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FOpenMPThread.hpp"

#include "../Containers/FOctree.hpp"

#include <omp.h>


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFMMAlgorithmThreadedInterval
* @brief
* Please read the license
*
*/
template<template< class ParticuleClass, class CellClass> class KernelClass, class ParticuleClass, class CellClass, int OctreeHeight, int SubtreeHeight>
class FFMMAlgorithmThreadedInterval : protected FAssertable{
    typedef typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator FOctreeIterator;

    static const int NbThreads = 1;
    static const int SizeInterval = 50;

    FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>* const tree;    //< The octree to work on
    KernelClass<ParticuleClass, CellClass>* kernels[NbThreads];                     //< The kernels

    FDEBUG_TIME(FTic counter);       //< In case of debug count the time

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree
      * @param inKernels the kernels
      * an assert is launched if one of the arguments is null
      */
    FFMMAlgorithmThreadedInterval(FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>* const inTree,
                  KernelClass<ParticuleClass, CellClass>* const inKernels)
        : tree(inTree) {
        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        for(int idxThread = 0 ; idxThread < NbThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticuleClass, CellClass>();
            *this->kernels[idxThread] = *inKernels;
        }

        FDEBUG_TRACE(FDebug::Controller.write("FFMMAlgorithmThreadedInterval\n"));
    }

    /** Default destructor */
    virtual ~FFMMAlgorithmThreadedInterval(){
        for(int idxThread = 0 ; idxThread < NbThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
    }

    /** To execute the fmm algorithm
      * Call this function to run the complete algo
      */
    void execute(){
        for(int idxThread = 0 ; idxThread < NbThreads ; ++idxThread){
            this->kernels[idxThread]->init();
        }

        bottomPass();
        upwardPass();

        downardPass();

        directPass();
    }

    /** P2M */
    void bottomPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG_TIME(counter.tic(););

        FOctreeIterator octreeIterator(tree);
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();

        omp_lock_t mutex;
        omp_init_lock(&mutex);
        bool stop = false;
        #pragma omp parallel shared(octreeIterator,mutex,stop) num_threads(NbThreads)
        {
            const int threadId = omp_get_thread_num();

            FOctreeIterator threadIter;
            int idxSizeInterval = 0;

            omp_set_lock(&mutex);
            while(!stop){
                threadIter = octreeIterator;
                for(idxSizeInterval = 1 ; idxSizeInterval < SizeInterval && octreeIterator.moveRight(); ++idxSizeInterval);
                if(idxSizeInterval != SizeInterval || !octreeIterator.moveRight()) stop = true;
                omp_unset_lock(&mutex);

                while(idxSizeInterval--){
                    // We need the current cell that represent the leaf
                    // and the list of particules
                    kernels[threadId]->P2M( threadIter.getCurrentCell() , threadIter.getCurrentList());
                    threadIter.moveRight();
                }
                omp_set_lock(&mutex);
            }
            omp_unset_lock(&mutex);
        }
        omp_destroy_lock(&mutex);

        FDEBUG_TIME(counter.tac(););
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (") FDEBUG_TIME(<< counter.elapsed() <<) FDEBUG_TRACE("s)\n"; )
    }

    /** M2M */
    void upwardPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG_TIME(counter.tic(););

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

                FOctreeIterator threadIter;
                int idxSizeInterval = 0;

                omp_set_lock(&mutex);
                // for each cells
                while(!stop){
                    threadIter = octreeIterator;
                    for(idxSizeInterval = 1 ; idxSizeInterval < SizeInterval && octreeIterator.moveRight(); ++idxSizeInterval);
                    if(idxSizeInterval != SizeInterval || !octreeIterator.moveRight()) stop = true;
                    omp_unset_lock(&mutex);

                    while(idxSizeInterval--){
                        // We need the current cell and the child
                        // child is an array (of 8 child) that may be null
                        kernels[threadId]->M2M( threadIter.getCurrentCell() , threadIter.getCurrentChild(), idxLevel);
                        threadIter.moveRight();
                    }

                    omp_set_lock(&mutex);
                }
                omp_unset_lock(&mutex);
            }
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
        }
        omp_destroy_lock(&mutex);

        FDEBUG_TIME(counter.tac(););
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (") FDEBUG_TIME(<< counter.elapsed() <<) FDEBUG_TRACE("s)\n"; )
    }

    /** M2L L2L */
    void downardPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG_TIME(counter.tic(););

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

                    FOctreeIterator threadIter;
                    int idxSizeInterval = 0;
                    CellClass* neighbors[208];

                    omp_set_lock(&mutex);
                    // for each cells
                    while(!stop){
                        threadIter = octreeIterator;
                        for(idxSizeInterval = 1 ; idxSizeInterval < SizeInterval && octreeIterator.moveRight(); ++idxSizeInterval);
                        if(idxSizeInterval != SizeInterval || !octreeIterator.moveRight()) stop = true;
                        omp_unset_lock(&mutex);

                        while(idxSizeInterval--){
                            const int counter = tree->getDistantNeighbors(neighbors, threadIter.getCurrentGlobalIndex(),idxLevel);
                            if(counter) kernels[threadId]->M2L( threadIter.getCurrentCell(), neighbors, counter, idxLevel);
                            threadIter.moveRight();
                        }

                        omp_set_lock(&mutex);
                    }
                    omp_unset_lock(&mutex);
                }
                octreeIterator.gotoLeft();
                octreeIterator.moveDown();
            }
            omp_destroy_lock(&mutex);
        }
        FDEBUG_TIME(counter.tac(););
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (") FDEBUG_TIME(<< counter.elapsed() <<) FDEBUG_TRACE("s)\n"; )

        FDEBUG_TRACE( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG_TIME(counter.tic(););
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

                    FOctreeIterator threadIter;
                    int idxSizeInterval = 0;

                    omp_set_lock(&mutex);
                    // for each cells
                    while(!stop){
                        threadIter = octreeIterator;
                        for(idxSizeInterval = 1 ; idxSizeInterval < SizeInterval && octreeIterator.moveRight(); ++idxSizeInterval);
                        if(idxSizeInterval != SizeInterval || !octreeIterator.moveRight()) stop = true;
                        omp_unset_lock(&mutex);

                        while(idxSizeInterval--){
                            kernels[threadId]->L2L( threadIter.getCurrentCell(),  threadIter.getCurrentChild(), idxLevel);
                            threadIter.moveRight();
                        }

                        omp_set_lock(&mutex);
                    }
                    omp_unset_lock(&mutex);
                }

                octreeIterator.gotoLeft();
                octreeIterator.moveDown();
            }
            omp_destroy_lock(&mutex);
        }

        FDEBUG_TIME(counter.tac(););
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (") FDEBUG_TIME(<< counter.elapsed() <<) FDEBUG_TRACE("s)\n"; )

    }

    /** P2P */
    void directPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG_TIME(counter.tic(););

        FOctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        const int heightMinusOne = OctreeHeight - 1;

        omp_lock_t mutex;
        omp_init_lock(&mutex);
        bool stop = false;
        #pragma omp parallel shared(octreeIterator,mutex,stop) num_threads(NbThreads)
        {
            const int threadId = omp_get_thread_num();
            // There is a maximum of 26 neighbors
            FList<ParticuleClass*>* neighbors[26];
            FOctreeIterator threadIter;
            int idxSizeInterval = 0;

            omp_set_lock(&mutex);
            // for each leafs
            while(!stop){
                threadIter = octreeIterator;
                for(idxSizeInterval = 1 ; idxSizeInterval < SizeInterval && octreeIterator.moveRight(); ++idxSizeInterval);
                if(idxSizeInterval != SizeInterval || !octreeIterator.moveRight()) stop = true;
                omp_unset_lock(&mutex);

                while(idxSizeInterval--){
                    kernels[threadId]->L2P(threadIter.getCurrentCell(), threadIter.getCurrentList());
                    // need the current particules and neighbors particules
                    const int counter = tree->getLeafsNeighbors(neighbors, threadIter.getCurrentGlobalIndex(),heightMinusOne);
                    kernels[threadId]->P2P( threadIter.getCurrentList() , neighbors, counter);
                    threadIter.moveRight();
                }

                omp_set_lock(&mutex);
            }
            omp_unset_lock(&mutex);
        }
        omp_destroy_lock(&mutex);

        FDEBUG_TIME(counter.tac(););
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (") FDEBUG_TIME(<< counter.elapsed() <<) FDEBUG_TRACE("s)\n"; )

    }

};


#endif //FFMMALGORITHMTHREADEDTHREADED_HPP

// [--LICENSE--]
