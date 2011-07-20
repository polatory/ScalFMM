#ifndef FFMMALGORITHMTHREADPROC_HPP
#define FFMMALGORITHMTHREADPROC_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"

#include "../Containers/FBoolArray.hpp"
#include "../Containers/FOctree.hpp"
#include "../Containers/FBufferVector.hpp"

#include "../Utils/FMpi.hpp"

#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadProc
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm with mpi.
* It just iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
*
* Of course this class does not deallocate pointer given in arguements.
*
* Threaded & based on the inspector-executor model
* schedule(runtime) export OMP_NUM_THREADS=2
* export OMPI_CXX=`which g++-4.4`
* mpirun -np 2 valgrind --suppressions=/usr/share/openmpi/openmpi-valgrind.supp
* --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes
* ./Tests/testFmmAlgorithmProc ../Data/testLoaderSmall.fma.tmp
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThreadProc : protected FAssertable {

    FMpi& app;                               //< The app to communicate

    OctreeClass* const tree;                 //< The octree to work on
    KernelClass** kernels;                   //< The kernels

    typename OctreeClass::Iterator* iterArray;
    int numberOfLeafs;                          //< To store the size at the previous level

    const int MaxThreads;               //< the max number of thread allowed by openmp

    const int nbProcess;                //< Number of process
    const int idProcess;                //< Id of current process

    const static int BufferSize = 2000;      //< To know max of the buffer we receive
    FBufferVector<BufferSize> * sendBuffer;  //< To put data to send into a buffer

    const int OctreeHeight;

    struct Interval{
        MortonIndex min;
        MortonIndex max;
    };
    Interval*const intervals;
    Interval*const intervalsPerLevel;
    Interval*const realIntervalsPerLevel;


    static void mpiassert(const int test, const unsigned line, const char* const message = 0){
        if(test != MPI_SUCCESS){
            printf("[ERROR] Test failled at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadProc(FMpi& inApp, OctreeClass* const inTree, KernelClass* const inKernels)
            : app(inApp), tree(inTree) , kernels(0), numberOfLeafs(0),
            MaxThreads(omp_get_max_threads()), nbProcess(inApp.processCount()), idProcess(inApp.processId()),
            sendBuffer(0), OctreeHeight(tree->getHeight()),intervals(new Interval[inApp.processCount()]),
            intervalsPerLevel(new Interval[inApp.processCount() * tree->getHeight()]),
            realIntervalsPerLevel(new Interval[inApp.processCount() * tree->getHeight()]){

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);

        this->kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        this->sendBuffer = new FBufferVector<BufferSize>[nbProcess];

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadProc\n");
        FDEBUG(FDebug::Controller << "Max threads = "  << MaxThreads << ", Procs = " << nbProcess << ".\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadProc(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
        delete [] this->sendBuffer;
        delete [] intervals;
        delete [] intervalsPerLevel;
        delete [] realIntervalsPerLevel;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        // Count leaf
        this->numberOfLeafs = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        intervals[idProcess].min = octreeIterator.getCurrentGlobalIndex();
        do{
            ++this->numberOfLeafs;
        } while(octreeIterator.moveRight());
        intervals[idProcess].max = octreeIterator.getCurrentGlobalIndex();

        iterArray = new typename OctreeClass::Iterator[numberOfLeafs];
        fassert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);

        mpiassert( MPI_Allgather( &intervals[idProcess], sizeof(Interval), MPI_BYTE, intervals, sizeof(Interval), MPI_BYTE, MPI_COMM_WORLD),  __LINE__ );
        for(int idxLevel = 0 ; idxLevel < OctreeHeight ; ++idxLevel){
            const int offset = idxLevel * nbProcess;
            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                intervalsPerLevel[offset + idxProc].max = intervals[idxProc].max >> (3 * (OctreeHeight - idxLevel - 1));
                intervalsPerLevel[offset + idxProc].min = intervals[idxProc].min >> (3 * (OctreeHeight - idxLevel - 1));
            }

            realIntervalsPerLevel[offset + 0] = intervalsPerLevel[offset + 0];
            for(int idxProc = 1 ; idxProc < nbProcess ; ++idxProc){
                realIntervalsPerLevel[offset + idxProc].min = FMath::Max( intervalsPerLevel[offset + idxProc].min,
                                                                          intervalsPerLevel[offset + idxProc - 1].max);
                realIntervalsPerLevel[offset + idxProc].max = intervalsPerLevel[offset + idxProc].max;
            }
        }

        // run
        FBoolArray alreadySent(nbProcess);
        const long P2PSent = preP2P(alreadySent);

        bottomPass();

        upwardPass();

        downardPass();

        directPass();

        // delete array
        delete [] iterArray;
        iterArray = 0;

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    long preP2P(FBoolArray& alreadySent){
        const long limite = 1 << (this->OctreeHeight - 1);
        long sent = 0;
        FVector<MPI_Request> requests;

        for(MortonIndex idxMort = intervals[idProcess].min ; idxMort <= intervals[idProcess].max ; ++idxMort){
            FTreeCoordinate center;
            center.setPositionFromMorton(idxMort, OctreeHeight - 1);

            bool hasLeaf = true;
            ContainerClass* leaf = 0;
            alreadySent.setToZeros();

            // We test all cells around
            for(long idxX = -1 ; idxX <= 1 ; ++idxX){
                if(!FMath::Between(center.getX() + idxX,0l,limite)) continue;

                for(long idxY = -1 ; idxY <= 1 ; ++idxY){
                    if(!FMath::Between(center.getY() + idxY,0l,limite)) continue;

                    for(long idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                        if(!FMath::Between(center.getZ() + idxZ,0l,limite)) continue;

                        // if we are not on the current cell
                        if( !(!idxX && !idxY && !idxZ) ){
                            const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                            const MortonIndex mortonOther = other.getMortonIndex(this->OctreeHeight - 1);
                            if(mortonOther < intervals[idProcess].min || intervals[idProcess].max < mortonOther){
                                // find the proc that need this information
                                int procToReceive = idProcess;
                                while( procToReceive != 0 && mortonOther < intervals[procToReceive].min){
                                    --procToReceive;
                                }
                                while( procToReceive != nbProcess - 1 && intervals[procToReceive].max < mortonOther){
                                    ++procToReceive;
                                }

                                if(!alreadySent.get(procToReceive)){
                                    alreadySent.set(procToReceive,true);
                                    ++sent;

                                    // get cell only when needed
                                    if(hasLeaf && !leaf){
                                        leaf = this->tree->getLeafSrc(mortonOther);
                                        if(!leaf) hasLeaf = false;
                                    }

                                    // add to list if not null
                                    MPI_Request req;
                                    if(leaf){
                                        MPI_Isend( leaf->data(), leaf->getSize() * sizeof(ParticleClass) , MPI_BYTE , procToReceive, 0, MPI_COMM_WORLD, &req);
                                    }
                                    else{
                                        MPI_Isend( 0, 0, MPI_BYTE , procToReceive, 0, MPI_COMM_WORLD, &req);
                                    }
                                    requests.push(req);
                                }
                            }
                        }
                    }
                }
            }
        }
        return sent;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Utils functions
    /////////////////////////////////////////////////////////////////////////////

    int getLeft(const int inSize) const {
        const float step = (float(inSize) / nbProcess);
        return int(FMath::Ceil(step * idProcess));
    }

    int getRight(const int inSize) const {
        const float step = (float(inSize) / nbProcess);
        const int res = int(FMath::Ceil(step * (idProcess+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    int getOtherRight(const int inSize,const int other) const {
        const float step = (float(inSize) / nbProcess);
        const int res = int(FMath::Ceil(step * (other+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    int getProc(const int position, const int inSize) const {
        const float step = (float(inSize) / nbProcess);
        return int(position/step);
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M Bottom Pass */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        int leafs = 0;
        do{
            iterArray[leafs++] = octreeIterator;
        } while(octreeIterator.moveRight());


        FDEBUG(FTic computationCounter);
        #pragma omp parallel
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);
        FDEBUG(FTic sendCounter);
        FDEBUG(FTic receiveCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);


        int previousLeftProc = idProcess;
        bool getWork = true;

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 && getWork ; --idxLevel ){
            // We do not touche cells that we are not responsible
            while( octreeIterator.getCurrentGlobalIndex() < (realIntervalsPerLevel[(idxLevel + 1) * nbProcess + idProcess].min >>3) && (getWork = octreeIterator.moveRight()) ){}
            if(!getWork) continue;

            // copy cells to work with
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;

            int needToSend = 0;
            MPI_Request* requestsSend = 0;
            int iterRequests = 0;
            MortonIndex* cellsIndexs= 0;

            MPI_Request* requestsRecv = 0;
            int* needToReceive = 0;
            int firstProcThatSend = idProcess + 1;
            int endProcThatSend = firstProcThatSend;

            // find cells to send
            if(idProcess != 0){
                while( needToSend < numberOfCells && iterArray[needToSend].getCurrentGlobalIndex() < this->intervals[idxLevel * nbProcess + idProcess - 1].max ){
                    ++needToSend;
                }
                if(needToSend){
                    // find the proc to send the data to
                    while(previousLeftProc != 0 && this->intervals[idxLevel * nbProcess + previousLeftProc - 1].max <= iterArray[needToSend].getCurrentGlobalIndex() ){
                        --previousLeftProc;
                    }

                    requestsSend = new MPI_Request[9*needToSend + 1];
                    cellsIndexs = new MortonIndex[needToSend];

                    MPI_Isend(&needToSend, 1, MPI_INT, previousLeftProc, 8, MPI_COMM_WORLD, &requestsSend[iterRequests++]);

                    for(int idxCell = 0 ; idxCell < needToSend ; ++idxCell){
                        CellClass** const child = iterArray[idxCell].getCurrentChild();

                        cellsIndexs[idxCell] = iterArray[idxCell].getCurrentGlobalIndex();
                        MPI_Isend(&cellsIndexs[idxCell], sizeof(MortonIndex), MPI_BYTE, previousLeftProc, 9, MPI_COMM_WORLD, &requestsSend[iterRequests++]);

                        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                            MPI_Isend(child[idxChild], child[idxChild]? sizeof(CellClass) : 0, MPI_BYTE, previousLeftProc, idxChild, MPI_COMM_WORLD, &requestsSend[iterRequests++]);
                        }
                    }
                }
            }
            if(idProcess != nbProcess - 1){
                while(firstProcThatSend < nbProcess &&
                      (this->intervals[(idxLevel+1) * nbProcess + idProcess + firstProcThatSend + 1].max >> 3) <= this->intervals[idxLevel * nbProcess + idProcess].max ){
                    ++firstProcThatSend;
                }
                while( endProcThatSend < nbProcess &&
                        (this->intervals[(idxLevel+1) * nbProcess + idProcess + endProcThatSend + 1].min >> 3) <= this->intervals[idxLevel * nbProcess + idProcess].max){
                    ++endProcThatSend;
                }
                if(firstProcThatSend != endProcThatSend){
                    needToReceive = new int[endProcThatSend - firstProcThatSend];
                    requestsRecv = new MPI_Request[endProcThatSend - firstProcThatSend];
                    for(int idxProc = firstProcThatSend ; idxProc < endProcThatSend ; ++idxProc ){
                        MPI_Irecv(&needToReceive[idxProc - firstProcThatSend], 1, MPI_INT, idxProc, 8, MPI_COMM_WORLD, &requestsRecv[idxProc - firstProcThatSend]);
                    }
                }
            }

            KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
            for(int idxCell = needToSend ; idxCell < numberOfCells ; ++idxCell){
                myThreadkernels.M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
            }


            if(iterRequests) MPI_Waitall(iterRequests, requestsSend, 0);

            delete[] requestsSend;
            delete[] cellsIndexs;

            delete[] needToReceive;
            delete[] requestsRecv;
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downard
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        { // first M2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);
            FDEBUG(FTic sendCounter);
            FDEBUG(FTic receiveCounter);
            FDEBUG(FTic waitingToReceiveCounter);
            FDEBUG(FTic waitSendCounter);
            FDEBUG(FTic findCounter);



            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Wait data to Receive : " << waitingToReceiveCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t\tTotal time to send in mpi "  << waitSendCounter.cumulated() << " s.\n" );
            FDEBUG( FDebug::Controller << "\t\t\tTotal time to find "  << findCounter.cumulated() << " s.\n" );
        }


        { // second L2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);
            FDEBUG(FTic sendCounter);
            FDEBUG(FTic receiveCounter);



            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
        }

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        // init
        const int LeafIndex = OctreeHeight - 1;
        const int SizeShape = 3*3*3;

        int shapeLeaf[SizeShape];
        memset(shapeLeaf,0,SizeShape*sizeof(int));

        struct LeafData{
            MortonIndex index;
            CellClass* cell;
            ContainerClass* targets;
            ContainerClass* sources;
        };
        LeafData* const leafsDataArray = new LeafData[this->numberOfLeafs];

        // split data
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            // to store which shape for each leaf
            typename OctreeClass::Iterator* const myLeafs = new typename OctreeClass::Iterator[this->numberOfLeafs];
            int*const shapeType = new int[this->numberOfLeafs];

            for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf){
                myLeafs[idxLeaf] = octreeIterator;

                const FTreeCoordinate& coord = octreeIterator.getCurrentCell()->getCoordinate();
                const int shape = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);
                shapeType[idxLeaf] = shape;

                ++shapeLeaf[shape];

                octreeIterator.moveRight();
            }

            int startPosAtShape[SizeShape];
            startPosAtShape[0] = 0;
            for(int idxShape = 1 ; idxShape < SizeShape ; ++idxShape){
                startPosAtShape[idxShape] = startPosAtShape[idxShape-1] + shapeLeaf[idxShape-1];
            }

            int idxInArray = 0;
            for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf, ++idxInArray){
                const int shapePosition = shapeType[idxInArray];

                leafsDataArray[startPosAtShape[shapePosition]].index = myLeafs[idxInArray].getCurrentGlobalIndex();
                leafsDataArray[startPosAtShape[shapePosition]].cell = myLeafs[idxInArray].getCurrentCell();
                leafsDataArray[startPosAtShape[shapePosition]].targets = myLeafs[idxInArray].getCurrentListTargets();
                leafsDataArray[startPosAtShape[shapePosition]].sources = myLeafs[idxInArray].getCurrentListSrc();

                ++startPosAtShape[shapePosition];
            }

            delete[] shapeType;
            delete[] myLeafs;
        }

        FDEBUG(FTic computationCounter);

        #pragma omp parallel
        {
            KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[26];
            MortonIndex neighborsIndex[26];
            int previous = 0;

            for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                const int endAtThisShape = shapeLeaf[idxShape] + previous;

                #pragma omp for schedule(dynamic)
                for(int idxLeafs = previous ; idxLeafs < endAtThisShape ; ++idxLeafs){
                    LeafData& currentIter = leafsDataArray[idxLeafs];
                    myThreadkernels.L2P(currentIter.cell, currentIter.targets);
                    // need the current particles and neighbors particles
                    const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, currentIter.index,LeafIndex);
                    myThreadkernels.P2P( currentIter.index,currentIter.targets, currentIter.sources , neighbors, neighborsIndex, counter);
                }

                previous = endAtThisShape;
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};






#endif //FFMMALGORITHMTHREAD_HPP

// [--LICENSE--]
