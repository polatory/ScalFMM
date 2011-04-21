#ifndef FFMMALGORITHMTHREADPROC_HPP
#define FFMMALGORITHMTHREADPROC_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"

#include "../Containers/FOctree.hpp"


//================================================================================================
#ifdef FUSE_MPI
// Compile by mpic++ testApplication.cpp ../Utils/FAssertable.cpp -o testApplication.exe
// run by mpirun -np 4 ./testApplication.exe
#include "../Utils/FMpiApplication.hpp"
typedef FMpiApplication ApplicationImplementation;
#else
// Compile by g++ testApplication.cpp ../Utils/FAssertable.cpp -o testApplication.exe
#include "../Utils/FSingleApplication.hpp"
typedef FSingleApplication ApplicationImplementation;
#endif
//================================================================================================

#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadProc
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
*
* Of course this class does not deallocate pointer given in arguements.
*
* Threaded & based on the inspector-executor model
* schedule(runtime)
*/
template<template< class ParticleClass, class CellClass, int OctreeHeight> class KernelClass,
        class ParticleClass, class CellClass,
        template<class ParticleClass> class LeafClass,
        int OctreeHeight, int SubtreeHeight>
class FFmmAlgorithmThreadProc : protected FAssertable, protected ApplicationImplementation{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator OctreeIterator;
    typedef KernelClass<ParticleClass, CellClass, OctreeHeight> Kernel;

    Octree* const tree;                  //< The octree to work on
    Kernel* kernels[FThreadNumbers];          //< The kernels

    FDEBUG(FTic counterTime);                //< In case of debug: to count the elapsed time
    FDEBUG(FTic computationCounter);     //< In case of debug: to  count computation time

    OctreeIterator* iterArray;

    static const int SizeShape = 3*3*3;
    int shapeLeaf[SizeShape];

    void run(){}

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadProc(Octree* const inTree, Kernel* const inKernels, const int inArgc, char ** const inArgv )
                      : ApplicationImplementation(inArgc,inArgv), tree(inTree) , iterArray(0) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass, OctreeHeight>(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadProc\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadProc(){
        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            delete this->kernels[idxThread];
        }
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            this->shapeLeaf[idxShape] = 0;
        }
        const int LeafIndex = OctreeHeight - 1;

        // Count leaf
        int leafs = 0;
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++leafs;
            const MortonIndex index = octreeIterator.getCurrentGlobalIndex();
            FTreeCoordinate coord;
            coord.setPositionFromMorton(index, LeafIndex);
            ++this->shapeLeaf[(coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3)];

        } while(octreeIterator.moveRight());
        iterArray = new OctreeIterator[leafs];
        assert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);

        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            this->kernels[idxThread]->init();
        }

        bottomPass();
        upwardPass();

        downardPass();

        directPass();

        delete [] iterArray;
        iterArray = 0;

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( counterTime.tic() );

        OctreeIterator octreeIterator(tree);
        const int nbProcess = processCount();
        const int idPorcess = processId();
        int leafs = 0;
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[leafs] = octreeIterator;
            ++leafs;
        } while(octreeIterator.moveRight());

        const float stepIdx = float(leafs) / nbProcess;
        const int startIdx = idPorcess*stepIdx;
        const int endIdx = (idPorcess+1)*stepIdx;

        FDEBUG(computationCounter.tic());
        #pragma omp parallel num_threads(FThreadNumbers)
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for
            for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSources());
            }

            #pragma omp single nowait
            {
                for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                    for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                        if(idxProc != idPorcess){
                            sendData(idxProc,sizeof(CellClass),iterArray[idxLeafs].getCurrentCell(),idxLeafs);
                        }
                    }
                }
            }

            #pragma omp single nowait
            {
                int needToReceive = leafs - (endIdx-startIdx);
                CellClass tempCell;
                int source = 0, tag = 0, filled = 0;

                while(needToReceive){
                    receiveData(sizeof(CellClass),&tempCell,&source,&tag,&filled);
                    if(filled){
                        *iterArray[tag].getCurrentCell() = tempCell;
                    }
                    --needToReceive;
                }
            }
        }
        FDEBUG(computationCounter.tac());

        processBarrier();

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        // Start from leal level - 1
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        OctreeIterator avoidGotoLeftIterator(octreeIterator);

        const int nbProcess = processCount();
        const int idPorcess = processId();

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            int leafs = 0;
            // for each cells
            do{
                iterArray[leafs] = octreeIterator;
                ++leafs;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();


            const float stepIdx = float(leafs) / nbProcess;
            const int startIdx = idPorcess*stepIdx;
            const int endIdx = (idPorcess+1)*stepIdx;

            FDEBUG(computationCounter.tic());

            #pragma omp parallel num_threads(FThreadNumbers)
            {
                Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for
                for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                    myThreadkernels->M2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                }

                #pragma omp single nowait
                {
                    FDEBUG(computationCounter.tac());
                    FDEBUG(totalComputation += computationCounter.elapsed());
                }

                // send computed data
                #pragma omp single nowait
                {
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                            if(idxProc != idPorcess){
                                sendData(idxProc,sizeof(CellClass),iterArray[idxLeafs].getCurrentCell(),idxLeafs);
                            }
                        }
                    }
                }

                // received computed data
                #pragma omp single nowait
                {
                    int needToReceive = leafs - (endIdx-startIdx);
                    CellClass tempCell;
                    int source = 0, tag = 0, filled = 0;

                    while(needToReceive){
                        receiveData(sizeof(CellClass),&tempCell,&source,&tag,&filled);
                        if(filled){
                            *iterArray[tag].getCurrentCell() = tempCell;
                        }
                        --needToReceive;
                    }
                }
            }
            processBarrier();
        }

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        { // first M2L
            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();
            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int nbProcess = processCount();
            const int idPorcess = processId();


            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                int leafs = 0;
                // for each cells
                do{
                    iterArray[leafs] = octreeIterator;
                    ++leafs;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const float stepIdx = float(leafs) / nbProcess;
                const int startIdx = idPorcess*stepIdx;
                const int endIdx = (idPorcess+1)*stepIdx;

                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    CellClass* neighbors[208];

                    #pragma omp for
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        const int counter = tree->getDistantNeighbors(neighbors,  iterArray[idxLeafs].getCurrentGlobalIndex(),idxLevel);
                        if(counter) myThreadkernels->M2L(  iterArray[idxLeafs].getCurrentCell() , neighbors, counter, idxLevel);
                    }

                    #pragma omp single nowait
                    {
                        for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                                if(idxProc != idPorcess){
                                    sendData(idxProc,sizeof(CellClass),iterArray[idxLeafs].getCurrentCell(),idxLeafs);
                                }
                            }
                        }
                    }

                    // received computed data
                    #pragma omp single nowait
                    {
                        int needToReceive = leafs - (endIdx-startIdx);
                        CellClass tempCell;
                        int source = 0, tag = 0, filled = 0;

                        while(needToReceive){
                            receiveData(sizeof(CellClass),&tempCell,&source,&tag,&filled);
                            if(filled){
                                *iterArray[tag].getCurrentCell() = tempCell;
                            }
                            --needToReceive;
                        }
                    }
                }
                processBarrier();
            }
        }
        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );

        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( totalComputation = 0 );
        { // second L2L
            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            const int nbProcess = processCount();
            const int idPorcess = processId();

            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                int leafs = 0;
                // for each cells
                do{
                    iterArray[leafs] = octreeIterator;
                    ++leafs;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const float stepIdx = float(leafs) / nbProcess;
                const int startIdx = idPorcess*stepIdx;
                const int endIdx = (idPorcess+1)*stepIdx;

                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        myThreadkernels->L2L( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                    }

                    #pragma omp single nowait
                    {
                        const int sizeBuffer = 8 * sizeof(CellClass);
                        char buffer[sizeBuffer];
                        for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                            int position=0;
                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                if(iterArray[idxLeafs].getCurrentChild()[idxChild]){
                                    memcpy(&buffer[position*sizeof(CellClass)],iterArray[idxLeafs].getCurrentChild()[idxChild],sizeof(CellClass));
                                    ++position;
                                }
                            }

                            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                                if(idxProc != idPorcess){
                                    sendData(idxProc,position*sizeof(CellClass),buffer,idxLeafs);
                                }
                            }
                        }
                    }

                    // received computed data
                    #pragma omp single nowait
                    {
                        int needToReceive = leafs - (endIdx-startIdx);
                        int source = 0, tag = 0, filled = 0;
                        const int sizeBuffer = 8 * sizeof(CellClass);
                        char buffer[sizeBuffer];

                        while(needToReceive){
                            receiveData(sizeBuffer,buffer,&source,&tag,&filled);
                            if(filled){
                                int position = 0;
                                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                    if(iterArray[tag].getCurrentChild()[idxChild]){
                                        memcpy(iterArray[tag].getCurrentChild()[idxChild],&buffer[position*sizeof(CellClass)],sizeof(CellClass));
                                        ++position;
                                    }
                                }
                            }
                            --needToReceive;
                        }
                    }
                }
                processBarrier();
            }
        }

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2P */
    void directPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );

        const int LeafIndex = OctreeHeight - 1;
        int leafs = 0;
        {
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leafs
            do{
                iterArray[leafs] = octreeIterator;
                ++leafs;
            } while(octreeIterator.moveRight());
        }

        FDEBUG(computationCounter.tic());

        #pragma omp parallel num_threads(FThreadNumbers)
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            FList<ParticleClass*>* neighbors[26];

            const float stepIdx = float(leafs) / processCount();
            const int startIdx = processId()*stepIdx;
            const int endIdx = (processId()+1)*stepIdx;

            #pragma omp for
            for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                myThreadkernels->L2P(iterArray[idxLeafs].getCurrentCell(), iterArray[idxLeafs].getCurrentListTargets());
                // need the current particles and neighbors particles
                const int counter = tree->getLeafsNeighbors(neighbors, iterArray[idxLeafs].getCurrentGlobalIndex(),LeafIndex);
                myThreadkernels->P2P( iterArray[idxLeafs].getCurrentListTargets(), iterArray[idxLeafs].getCurrentListSources() , neighbors, counter);
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** This function test the octree to be sure that the fmm algorithm
      * has worked completly.
      */
    void ValidateFMMAlgoProc(){
        std::cout << "Check Result\n";
        int NbPart = 0;
        int NbLeafs = 0;
        { // Check that each particle has been summed with all other
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSources()->getSize();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
        }
        {
            const float stepIdx = float(NbLeafs) / processCount();
            const int startIdx = processId()*stepIdx;
            const int endIdx = (processId()+1)*stepIdx;
            // Check that each particle has been summed with all other
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                typename FList<ParticleClass*>::BasicIterator iter(*octreeIterator.getCurrentListTargets());

                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSources());

                while( iter.isValide() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( (!isUsingTsm && iter.value()->getDataDown() != NbPart - 1) ||
                        (isUsingTsm && iter.value()->getDataDown() != NbPart) ){
                        std::cout << "Problem L2P + P2P, value on particle is : " << iter.value()->getDataDown() << "\n";
                    }
                    iter.progress();
                }
                octreeIterator.moveRight();
            }
        }

        std::cout << "Done\n";
    }

    void print(){
        OctreeIterator octreeIterator(tree);
        for(int idxLevel = OctreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            do{
                std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
            } while(octreeIterator.moveRight());
            std::cout << "\n";
            octreeIterator.gotoLeft();
            octreeIterator.moveDown();
        }
    }

};






#endif //FFMMALGORITHMTHREAD_HPP

// [--LICENSE--]
