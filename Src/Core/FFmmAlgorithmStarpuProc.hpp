#ifndef FFMMALGORITHM_HPP
#define FFMMALGORITHM_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"
#include "../Containers/FBoolArray.hpp"

#include <limits>
#include <starpu.h>
#include <starpu_mpi.h>
#include <mpi.h>
#include <map>

#include "../Utils/FMpi.hpp"


template <class RealCell>
class StarCell : public RealCell {
public:
    starpu_data_handle handle;

    void initHandle(){
        memset(&handle, 0, sizeof(starpu_data_handle));
        starpu_variable_data_register(&handle, 0, (uintptr_t)static_cast<RealCell*>(this), sizeof(RealCell));
    }
};


template < class RealContainer >
class StarContainer : public RealContainer {
public:
    starpu_data_handle handle;

    void initHandle(){
        memset(&handle, 0, sizeof(starpu_data_handle));
        starpu_vector_data_register(&handle, 0, (uintptr_t)&(*this)[0], RealContainer::getSize(), sizeof((*this)[0]));
    }
};

template<class T>
class StarVector {
protected:
    T* array;
    int size;

    StarVector& operator=(const StarVector&){ return *this;}
    StarVector(const StarVector&){}

public:
    StarVector(T*const inData, const unsigned int inSize) : array(inData), size(inSize) {
    }

    virtual ~StarVector(){
    }

    int getSize() const {
        return size;
    }

    /** To get the entire array
      * @return the array allocated by the vector
      */
    T* data(){
        return this->array;
    }

    /** To get the entire array
      * @return the array allocated by the vector
      */
    const T* data() const{
        return this->array;
    }


    class BasicIterator{
    protected:
        StarVector* const vector;  /**< the vector to work on*/
        int index;                 /**< the current node*/

    public:

        virtual ~BasicIterator(){}

        BasicIterator(StarVector<T>& inVector) : vector(&inVector), index(0){}

        void gotoNext(){
            ++this->index;
        }

        bool hasNotFinished() const{
            return this->index < this->vector->size;
        }

        T& data(){
            return this->vector->array[this->index];
        }

        const T& data() const{
            return this->vector->array[this->index];
        }

    };
    friend class BasicIterator;

    /** This class is a basic const iterator
      * it uses a const vector to work on
      */
    class ConstBasicIterator{
    protected:
        const StarVector* const vector;  /**< the vector to work on*/
        int index;              /**< the current node*/

    public:

        virtual ~ConstBasicIterator(){}

        ConstBasicIterator(const StarVector<T>& inVector) : vector(&inVector), index(0){}

        void gotoNext(){
            ++this->index;
        }

        bool hasNotFinished() const{
            return this->index < this->vector->size;
        }

        const T& data() const{
            return this->vector->array[this->index];
        }

    };
    friend class ConstBasicIterator;

};






/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmStarpu
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmStarpuProc : protected FAssertable{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;

    //////////////////////////////////////////////////////////////////
    // MPI Specific
    //////////////////////////////////////////////////////////////////

    int numberOfLeafs;
    int nbProcess;                //< Number of process
    int idProcess;                //< Id of current process

    struct Interval{
        MortonIndex min;
        MortonIndex max;
    };
    Interval* intervals;
    Interval* workingIntervalsPerLevel;

    MPI_Request* requestsP2P;
    MPI_Status* statusP2P;
    int iterRequestP2P;
    int nbMessagesToRecvP2P;

    ParticleClass** sendBufferP2P;
    ParticleClass** recvBufferP2P;

    int* globalReceiveMapP2P;

    FBoolArray* leafsNeedOtherP2P;

    OctreeClass* otherP2Ptree;


    CellClass* emptyCellsM2M;
    CellClass* recvCellsM2M;

    CellClass* cellsFromOtherM2L;
    std::map<MortonIndex, CellClass*>* tempTreeM2L;



    void initProc(const FMpi::FComm& comm){
        nbProcess = comm.processCount();
        idProcess = comm.processId();

        intervals = new Interval[comm.processCount()];
        workingIntervalsPerLevel = new Interval[comm.processCount() * tree->getHeight()];

        Interval myLastInterval;
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            myLastInterval.min = octreeIterator.getCurrentGlobalIndex();
            do{
                ++this->numberOfLeafs;
            } while(octreeIterator.moveRight());
            myLastInterval.max = octreeIterator.getCurrentGlobalIndex();
        }

        // We get the min/max indexes from each procs
        mpiassert( MPI_Allgather( &myLastInterval, sizeof(Interval), MPI_BYTE, intervals, sizeof(Interval), MPI_BYTE, MPI_COMM_WORLD),  __LINE__ );

        Interval myIntervals[OctreeHeight];
        myIntervals[OctreeHeight - 1] = myLastInterval;
        for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 0 ; --idxLevel){
            myIntervals[idxLevel].min = myIntervals[idxLevel+1].min >> 3;
            myIntervals[idxLevel].max = myIntervals[idxLevel+1].max >> 3;
        }
        if(idProcess != 0){
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            octreeIterator.moveUp();

            MortonIndex currentLimit = intervals[idProcess-1].max >> 3;

            for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 1 ; --idxLevel){
                while(octreeIterator.getCurrentGlobalIndex() <= currentLimit){
                    if( !octreeIterator.moveRight() ) break;
                }
                myIntervals[idxLevel].min = octreeIterator.getCurrentGlobalIndex();
                octreeIterator.moveUp();
                currentLimit >>= 3;
            }
        }

        // We get the min/max indexes from each procs
        mpiassert( MPI_Allgather( myIntervals, sizeof(Interval) * OctreeHeight, MPI_BYTE,
                                 workingIntervalsPerLevel, sizeof(Interval) * OctreeHeight, MPI_BYTE, MPI_COMM_WORLD),  __LINE__ );


        // P2P Send/Receive Data
        requestsP2P = new MPI_Request[2 * nbProcess];
        statusP2P = new MPI_Status[2 * nbProcess];
        iterRequestP2P = 0;
        nbMessagesToRecvP2P = 0;

        sendBufferP2P = new ParticleClass*[nbProcess];
        memset(sendBufferP2P, 0, sizeof(ParticleClass*) * nbProcess);

        recvBufferP2P = new ParticleClass*[nbProcess];
        memset(recvBufferP2P, 0, sizeof(ParticleClass*) * nbProcess);

        globalReceiveMapP2P = new int[nbProcess * nbProcess];
        memset(globalReceiveMapP2P, 0, sizeof(int) * nbProcess * nbProcess);

        leafsNeedOtherP2P = new FBoolArray(this->numberOfLeafs);

        otherP2Ptree = new OctreeClass( tree->getHeight(), tree->getSubHeight(), tree->getBoxWidth(), tree->getBoxCenter() );

        emptyCellsM2M = new CellClass[8 * OctreeHeight];
        recvCellsM2M = new CellClass[8 * OctreeHeight];

        tempTreeM2L = new std::map<MortonIndex, CellClass*>[OctreeHeight];

        cellsFromOtherM2L = 0;
    }

    void releaseProc(){
        delete [] intervals;
        delete [] workingIntervalsPerLevel;

        delete [] requestsP2P;
        delete [] statusP2P;

        delete [] sendBufferP2P;

        delete [] recvBufferP2P;

        delete [] globalReceiveMapP2P;

        delete leafsNeedOtherP2P;

        delete otherP2Ptree;

        delete [] emptyCellsM2M;
        delete [] recvCellsM2M;

        delete [] tempTreeM2L;

        delete [] cellsFromOtherM2L;
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

    //////////////////////////////////////////////////////////////////
    // Codelets
    //////////////////////////////////////////////////////////////////

    starpu_codelet p2m_cl;
    starpu_codelet p2p_cl[27];
    starpu_codelet m2m_cl[9];
    starpu_codelet m2l_cl[190];
    starpu_codelet l2l_cl[9];
    starpu_codelet l2p_cl;

    void initCodelets(){
        memset(&p2m_cl, 0, sizeof(p2m_cl));
        p2m_cl.where = STARPU_CPU;
        p2m_cl.cpu_func = p2m_cpu;
        p2m_cl.nbuffers = 2;

        memset(&p2p_cl, 0, 27 * sizeof(p2p_cl[0]));
        for(int idxNeigh = 0 ; idxNeigh <= 26 ; ++idxNeigh){
            p2p_cl[idxNeigh].where = STARPU_CPU;
            p2p_cl[idxNeigh].cpu_func = p2p_cpu;
            p2p_cl[idxNeigh].nbuffers = idxNeigh + 1;
        }

        memset(&m2l_cl, 0, 190 * sizeof(m2l_cl[0]));
        for(int idxNeigh = 0 ; idxNeigh <= 189 ; ++idxNeigh){
            m2l_cl[idxNeigh].where = STARPU_CPU;
            m2l_cl[idxNeigh].cpu_func = m2l_cpu;
            m2l_cl[idxNeigh].nbuffers = idxNeigh + 1;
        }

        memset(&l2p_cl, 0, sizeof(l2p_cl));
        l2p_cl.where = STARPU_CPU;
        l2p_cl.cpu_func = l2p_cpu;
        l2p_cl.nbuffers = 2;

        memset(&m2m_cl, 0, 9 * sizeof(m2m_cl[0]));
        memset(&l2l_cl, 0, 9 * sizeof(m2m_cl[0]));
        for(int idxChild = 0 ; idxChild <= 8 ; ++idxChild){
            m2m_cl[idxChild].where = STARPU_CPU;
            m2m_cl[idxChild].cpu_func = m2m_cpu;
            m2m_cl[idxChild].nbuffers = idxChild + 1;

            l2l_cl[idxChild].where = STARPU_CPU;
            l2l_cl[idxChild].cpu_func = l2l_cpu;
            l2l_cl[idxChild].nbuffers = idxChild + 1;
        }
    }

    void releaseCodelets(){
    }

    //////////////////////////////////////////////////////////////////
    // Codelets Args
    //////////////////////////////////////////////////////////////////

    struct M2lProperties {
        int level;
        int bufferSize;
    };


    int* levels;
    int p2pNeighbors[27];
    M2lProperties** m2lNeighbors;
    int* bitsDescriptors;

    void initArgs(){
        levels = new int[OctreeHeight];
        for(int idxLevel = 0 ; idxLevel < OctreeHeight ; ++idxLevel ){
            levels[idxLevel] = idxLevel;
        }

        for(int idxP2pNeighbors = 0 ; idxP2pNeighbors <= 26 ; ++idxP2pNeighbors ){
            p2pNeighbors[idxP2pNeighbors] = idxP2pNeighbors;
        }

        m2lNeighbors = new M2lProperties*[OctreeHeight];
        for(int idxLevel = 0 ; idxLevel < OctreeHeight ; ++idxLevel){
            m2lNeighbors[idxLevel] = new M2lProperties[190];
            for(int idxm2lNeighbors = 0 ; idxm2lNeighbors <= 189 ; ++idxm2lNeighbors){
                m2lNeighbors[idxLevel][idxm2lNeighbors].level = idxLevel;
                m2lNeighbors[idxLevel][idxm2lNeighbors].bufferSize = idxm2lNeighbors;
            }
        }

        bitsDescriptors = new int[OctreeHeight * 256];
        for(int idxLevel = 0 ; idxLevel < OctreeHeight ; ++idxLevel){
            for(int idxBit = 0 ; idxBit < 256 ; ++idxBit){
                bitsDescriptors[idxLevel * 256 + idxBit] = (idxBit<<8) | idxLevel;
            }
        }
    }

    void releaseArgs(){
        delete[] levels;

        for(int idxLevel = 0 ; idxLevel < OctreeHeight ; ++idxLevel){
            delete[] m2lNeighbors[idxLevel];
        }
        delete[] m2lNeighbors;

        delete[] bitsDescriptors;
    }

    //////////////////////////////////////////////////////////////////
    // Init handles
    //////////////////////////////////////////////////////////////////

    void initHandles(){
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
        // init leaf handle
        do{
            octreeIterator.getCurrentLeaf()->getSrc()->initHandle();
            if(octreeIterator.getCurrentLeaf()->getSrc() != octreeIterator.getCurrentLeaf()->getTargets()){
                octreeIterator.getCurrentLeaf()->getTargets()->initHandle();
            }
        } while(octreeIterator.moveRight());

        octreeIterator = avoidGotoLeftIterator;

        // init cells handle
        for(int idxLevel = OctreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            do{
                octreeIterator.getCurrentCell()->initHandle();
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }
    }

    void releaseHandles(){
    }

    //////////////////////////////////////////////////////////////////
    // Init Kernels
    //////////////////////////////////////////////////////////////////

    void initKernels(){
        globalKernels = new KernelClass*[starpu_worker_get_count()];
        memset(globalKernels, 0, sizeof(KernelClass*) * starpu_worker_get_count());

        for(unsigned int workerid = 0; workerid < starpu_worker_get_count(); ++workerid){
            if( starpu_worker_get_type(workerid) == STARPU_CPU_WORKER ){
                globalKernels[workerid] = new KernelClass(*kernels);
            }
        }
    }

    void releaseKernels(){
        for(unsigned int workerid = 0; workerid < starpu_worker_get_count(); ++workerid){
            delete globalKernels[workerid];
        }
        delete[] globalKernels;
    }


    static void mpiassert(const int test, const unsigned line, const char* const message = 0){
        if(test != MPI_SUCCESS){
            printf("[ERROR] Test failled at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

    Interval& getWorkingInterval(const int level, const int proc){
        return workingIntervalsPerLevel[OctreeHeight * proc + level];
    }


public:

    Interval& getWorkingInterval(const int level){
        return getWorkingInterval(level, idProcess);
    }

    bool hasWorkAtLevel(const int level){
        return idProcess == 0 || getWorkingInterval(level, idProcess - 1).max < getWorkingInterval(level, idProcess).max;
    }

    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmStarpuProc(const FMpi::FComm& comm, OctreeClass* const inTree, KernelClass* const inKernels)
                      : tree(inTree) , kernels(inKernels), OctreeHeight(tree->getHeight()) {
        FDEBUG(FDebug::Controller << "FFmmAlgorithmStarpu\n");

        starpu_init(NULL);
        starpu_mpi_initialize();

        initProc(comm);

        initCodelets();
        initArgs();
        initHandles();
        initKernels();
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmStarpuProc(){
        starpu_mpi_shutdown();
        starpu_shutdown();

        releaseCodelets();
        releaseArgs();
        releaseHandles();
        releaseKernels();
        releaseProc();
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        P2MTasks();
        P2PTasks();
        M2MTasks();
        M2LTasks();
        L2LTasks();
        L2PTasks();

        P2PTasksOther();

        starpu_task_wait_for_all();

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void P2MTasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart P2M\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            //kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
            starpu_insert_task( &p2m_cl, STARPU_RW, octreeIterator.getCurrentCell()->handle, STARPU_R, octreeIterator.getCurrentLeaf()->getSrc()->handle, 0);
        } while(octreeIterator.moveRight());

        FDEBUG( FDebug::Controller << "\tFinished (@P2M = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2P
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void P2PTasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart P2P\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        { // Now send and receive P2P
            FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );

            // Box limite
            const long limite = 1 << (this->OctreeHeight - 1);
            // pointer to send
            typename OctreeClass::Iterator* toSend[nbProcess];
            memset(toSend, 0, sizeof(typename OctreeClass::Iterator*) * nbProcess );

            int sizeToSend[nbProcess];
            memset(sizeToSend, 0, sizeof(int) * nbProcess);
            // index
            int indexToSend[nbProcess];
            memset(indexToSend, 0, sizeof(int) * nbProcess);
            // index
            int partsToSend[nbProcess];
            memset(partsToSend, 0, sizeof(int) * nbProcess);

            // To know if a leaf has been already sent to a proc
            int alreadySent[nbProcess];

            MortonIndex indexesNeighbors[26];

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            int idxLeaf = 0;

            do {
                FTreeCoordinate center;
                center.setPositionFromMorton(octreeIterator.getCurrentGlobalIndex(), OctreeHeight - 1);

                memset(alreadySent, 0, sizeof(int) * nbProcess);
                bool needOther = false;

                const int neighCount = getNeighborsIndexes(octreeIterator.getCurrentGlobalIndex(), limite, indexesNeighbors);

                for(int idxNeigh = 0 ; idxNeigh < neighCount ; ++idxNeigh){
                    if(indexesNeighbors[idxNeigh] < intervals[idProcess].min || intervals[idProcess].max < indexesNeighbors[idxNeigh]){
                        needOther = true;

                        // find the proc that need this information
                        int procToReceive = idProcess;
                        while( procToReceive != 0 && indexesNeighbors[idxNeigh] < intervals[procToReceive].min){
                            --procToReceive;
                        }

                        while( procToReceive != nbProcess - 1 && intervals[procToReceive].max < indexesNeighbors[idxNeigh]){
                            ++procToReceive;
                        }

                        if( !alreadySent[procToReceive] && intervals[procToReceive].min <= indexesNeighbors[idxNeigh] && indexesNeighbors[idxNeigh] <= intervals[procToReceive].max){

                            alreadySent[procToReceive] = 1;
                            if(indexToSend[procToReceive] ==  sizeToSend[procToReceive]){
                                const int previousSize = sizeToSend[procToReceive];
                                sizeToSend[procToReceive] = FMath::Max(10*int(sizeof(typename OctreeClass::Iterator)), int(sizeToSend[procToReceive] * 1.5));
                                typename OctreeClass::Iterator* temp = toSend[procToReceive];
                                toSend[procToReceive] = reinterpret_cast<typename OctreeClass::Iterator*>(new char[sizeof(typename OctreeClass::Iterator) * sizeToSend[procToReceive]]);
                                memcpy(toSend[procToReceive], temp, previousSize * sizeof(typename OctreeClass::Iterator));
                                delete[] reinterpret_cast<char*>(temp);
                            }
                            toSend[procToReceive][indexToSend[procToReceive]++] = octreeIterator;
                            partsToSend[procToReceive] += octreeIterator.getCurrentListSrc()->getSize();
                        }
                    }
                }

                if(needOther){
                    leafsNeedOtherP2P->set(idxLeaf,true);
                }
                else {
                    ContainerClass* neighbors[26];
                    MortonIndex neighborsIndex[26];

                    const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, octreeIterator.getCurrentGlobalIndex(), OctreeHeight - 1);

                    //kernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                    struct starpu_task* const task = starpu_task_create();

                    task->buffers[0].handle = octreeIterator.getCurrentLeaf()->getTargets()->handle;
                    task->buffers[0].mode = STARPU_RW;

                    for(int idxCounterNeigh = 0 ; idxCounterNeigh < counter ; ++idxCounterNeigh){
                        task->buffers[idxCounterNeigh + 1].handle = neighbors[idxCounterNeigh]->handle;
                        task->buffers[idxCounterNeigh + 1].mode   = STARPU_R;
                    }

                    task->cl = &p2p_cl[counter];

                    task->cl_arg = &p2pNeighbors[counter];
                    task->cl_arg_size = sizeof(int);

                    starpu_task_submit(task);
                }
                ++idxLeaf;
            } while(octreeIterator.moveRight());

            mpiassert( MPI_Allgather( partsToSend, nbProcess, MPI_INT, globalReceiveMapP2P, nbProcess, MPI_INT, MPI_COMM_WORLD),  __LINE__ );

            // Prepare receive
            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                if(globalReceiveMapP2P[idxProc * nbProcess + idProcess]){
                    recvBufferP2P[idxProc] = reinterpret_cast<ParticleClass*>(new char[sizeof(ParticleClass) * globalReceiveMapP2P[idxProc * nbProcess + idProcess]]);

                    mpiassert( MPI_Irecv(recvBufferP2P[idxProc], globalReceiveMapP2P[idxProc * nbProcess + idProcess]*sizeof(ParticleClass), MPI_BYTE,
                                         idxProc, FMpi::TagFmmP2P, MPI_COMM_WORLD, &requestsP2P[iterRequestP2P++]) , __LINE__ );
                }
            }
            nbMessagesToRecvP2P = iterRequestP2P;
            // Prepare send
            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                if(indexToSend[idxProc] != 0){
                    sendBufferP2P[idxProc] = reinterpret_cast<ParticleClass*>(new char[sizeof(ParticleClass) * partsToSend[idxProc]]);

                    int currentIndex = 0;
                    for(int idxLeaf = 0 ; idxLeaf < indexToSend[idxProc] ; ++idxLeaf){
                        memcpy(&sendBufferP2P[idxProc][currentIndex], toSend[idxProc][idxLeaf].getCurrentListSrc()->data(),
                               sizeof(ParticleClass) * toSend[idxProc][idxLeaf].getCurrentListSrc()->getSize() );
                        currentIndex += toSend[idxProc][idxLeaf].getCurrentListSrc()->getSize();
                    }

                    mpiassert( MPI_Isend( sendBufferP2P[idxProc], sizeof(ParticleClass) * partsToSend[idxProc] , MPI_BYTE ,
                                         idxProc, FMpi::TagFmmP2P, MPI_COMM_WORLD, &requestsP2P[iterRequestP2P++]) , __LINE__ );

                }
            }


            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                delete [] reinterpret_cast<char*>(toSend[idxProc]);
                toSend[idxProc] = 0;
            }
        }

        FDEBUG( FDebug::Controller << "\tFinished (@P2P = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    void P2PTasksOther(){
        int complete = 0;
        while( complete != iterRequestP2P){

            int indexMessage[nbProcess * 2];
            memset(indexMessage, 0, sizeof(int) * nbProcess * 2);
            int countMessages = 0;
            // Wait data
            MPI_Waitsome(iterRequestP2P, requestsP2P, &countMessages, indexMessage, statusP2P);
            complete += countMessages;


            for(int idxRcv = 0 ; idxRcv < countMessages ; ++idxRcv){
                if( indexMessage[idxRcv] < nbMessagesToRecvP2P ){
                    const int idxProc = statusP2P[idxRcv].MPI_SOURCE;
                    for(int idxPart = 0 ; idxPart < globalReceiveMapP2P[idxProc * nbProcess + idProcess] ; ++idxPart){
                        otherP2Ptree->insert(recvBufferP2P[idxProc][idxPart]);
                    }
                    delete [] reinterpret_cast<char*>(recvBufferP2P[idxProc]);
                    recvBufferP2P[idxProc] = 0;
                }
            }
        }


        for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
            delete [] reinterpret_cast<char*>(sendBufferP2P[idxProc]);
            sendBufferP2P[idxProc] = 0;
        }

        {
            typename OctreeClass::Iterator octreeIterator(otherP2Ptree);
            octreeIterator.gotoBottomLeft();
            // init leaf handle
            do{
                octreeIterator.getCurrentLeaf()->getSrc()->initHandle();
            } while(octreeIterator.moveRight());
        }


        { // Exec all P2P
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            int idxLeaf = 0;
            const long limite = 1 << (this->OctreeHeight - 1);
            // for each leafs
            do{
                if( leafsNeedOtherP2P->get(idxLeaf)){
                    ContainerClass* neighbors[26];
                    MortonIndex neighborsIndex[26];
                    int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, octreeIterator.getCurrentGlobalIndex(), OctreeHeight - 1);


                    MortonIndex indexesNeighbors[26];
                    const int neighCount = getNeighborsIndexes(octreeIterator.getCurrentGlobalIndex(), limite, indexesNeighbors);

                    for(int idxNeigh = 0 ; idxNeigh < neighCount ; ++idxNeigh){
                        if(indexesNeighbors[idxNeigh] < intervals[idProcess].min || intervals[idProcess].max < indexesNeighbors[idxNeigh]){
                            ContainerClass*const hypotheticNeighbor = otherP2Ptree->getLeafSrc(indexesNeighbors[idxNeigh]);
                            if(hypotheticNeighbor){
                                neighbors[counter] = hypotheticNeighbor;
                                neighborsIndex[counter] = indexesNeighbors[idxNeigh];
                                ++counter;
                            }
                        }
                    }

                    //kernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                    struct starpu_task* const task = starpu_task_create();

                    task->buffers[0].handle = octreeIterator.getCurrentLeaf()->getTargets()->handle;
                    task->buffers[0].mode = STARPU_RW;

                    for(int idxCounterNeigh = 0 ; idxCounterNeigh < counter ; ++idxCounterNeigh){
                        task->buffers[idxCounterNeigh + 1].handle = neighbors[idxCounterNeigh]->handle;
                        task->buffers[idxCounterNeigh + 1].mode   = STARPU_R;
                    }

                    task->cl = &p2p_cl[counter];

                    task->cl_arg = &p2pNeighbors[counter];
                    task->cl_arg_size = sizeof(int);

                    starpu_task_submit(task);
                }

                ++idxLeaf;
            } while(octreeIterator.moveRight());
        }

    }

    /////////////////////////////////////////////////////////////////////////////
    // M2M
    /////////////////////////////////////////////////////////////////////////////


    /** M2M */
    void M2MTasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart M2M\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        int firstProcThatSend = idProcess + 1;

        int sendToProc = idProcess;

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            if(idProcess != 0
                    && getWorkingInterval((idxLevel+1), idProcess).max <= getWorkingInterval((idxLevel+1), idProcess - 1).max){
                break;
            }

            typename OctreeClass::Iterator previous = octreeIterator;
            while(octreeIterator.getCurrentGlobalIndex() < getWorkingInterval(idxLevel, idProcess).min){
                previous = octreeIterator;
                octreeIterator.moveRight();
            }

            if(idProcess != 0
                    && (getWorkingInterval((idxLevel+1), idProcess).min >>3) <= (getWorkingInterval((idxLevel+1), idProcess - 1).max >>3)){

                const CellClass* const* const child = previous.getCurrentChild();
                const int startIndex = getWorkingInterval((idxLevel+1), idProcess).min & 0x7;
                int endIndex = 8;
                if(idProcess != nbProcess - 1
                        && previous.getCurrentGlobalIndex() == (getWorkingInterval((idxLevel+1), idProcess + 1).min >>3)){
                    endIndex = getWorkingInterval((idxLevel+1), idProcess + 1).min & 0x7;
                }

                while( sendToProc && previous.getCurrentGlobalIndex() == getWorkingInterval(idxLevel , sendToProc - 1).max){
                    --sendToProc;
                }

                for(int idxChild = startIndex ; idxChild < endIndex ; ++idxChild){
                    if( child[idxChild] && getWorkingInterval((idxLevel+1), idProcess).min <= child[idxChild]->getMortonIndex() ){
                        starpu_mpi_isend_detached(child[idxChild]->handle ,
                                                  sendToProc, (((((previous.getCurrentGlobalIndex() << 3) | idxChild) << 5) | idxLevel) << 5) | FMpi::TagFmmM2M , MPI_COMM_WORLD , 0, 0 );
                    }
                    else {
                        emptyCellsM2M[idxLevel * 8 + idxChild].initHandle();
                        starpu_mpi_isend_detached(emptyCellsM2M[idxLevel * 8 + idxChild].handle ,
                                                  sendToProc, (((((previous.getCurrentGlobalIndex() << 3) | idxChild) << 5) | idxLevel) << 5) | FMpi::TagFmmM2M , MPI_COMM_WORLD , 0, 0 );
                    }
                }
            }

            do {
                //kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                struct starpu_task* const task = starpu_task_create();

                task->buffers[0].handle = octreeIterator.getCurrentCell()->handle;
                task->buffers[0].mode = STARPU_RW;

                int idxRealNeig   = 1;
                int bitNeigIndex  = 1;
                int bitDescriptor = 0;

                CellClass*const*const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if(child[idxChild]){
                        task->buffers[idxRealNeig].handle = child[idxChild]->handle;
                        task->buffers[idxRealNeig].mode = STARPU_R;
                        ++idxRealNeig;
                        bitDescriptor |= bitNeigIndex;
                    }
                    bitNeigIndex <<= 1;
                }

                task->cl = &m2m_cl[idxRealNeig - 1];

                task->cl_arg = &bitsDescriptors[bitDescriptor + idxLevel * 256];
                task->cl_arg_size = sizeof(int);

                starpu_task_submit(task);
            } while( octreeIterator.moveRight() );


            bool needToRecv = false;
            int endProcThatSend = firstProcThatSend;

            if(idProcess != nbProcess - 1){
                while(firstProcThatSend < nbProcess
                      && getWorkingInterval((idxLevel+1), firstProcThatSend).max < getWorkingInterval((idxLevel+1), idProcess).max){
                    ++firstProcThatSend;
                }

                if(firstProcThatSend < nbProcess &&
                        (getWorkingInterval((idxLevel+1), firstProcThatSend).min >>3) <= (getWorkingInterval((idxLevel+1) , idProcess).max>>3) ){

                    endProcThatSend = firstProcThatSend;

                    while( endProcThatSend < nbProcess &&
                          (getWorkingInterval((idxLevel+1) ,endProcThatSend).min >>3) <= (getWorkingInterval((idxLevel+1) , idProcess).max>>3)){
                        ++endProcThatSend;
                    }
                    if(firstProcThatSend != endProcThatSend){
                        needToRecv = true;
                        const int startIndex = getWorkingInterval((idxLevel+1) ,firstProcThatSend).min & 0x7;

                        for(int idxChild = startIndex ; idxChild < 8 ; ++idxChild){
                            recvCellsM2M[idxLevel * 8 + idxChild].initHandle();
                            starpu_mpi_irecv_detached(recvCellsM2M[idxLevel * 8 + idxChild].handle,
                                                      MPI_ANY_SOURCE, (((((octreeIterator.getCurrentGlobalIndex() << 3) | idxChild) << 5) | idxLevel) << 5) | FMpi::TagFmmM2M , MPI_COMM_WORLD , 0, 0);
                        }
                    }
                }
            }


            struct starpu_task* const task = starpu_task_create();

            task->buffers[0].handle = octreeIterator.getCurrentCell()->handle;
            task->buffers[0].mode = STARPU_RW;

            int idxRealNeig   = 1;
            int bitNeigIndex  = 1;
            int bitDescriptor = 0;

            CellClass*const*const child = octreeIterator.getCurrentChild();
            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if(child[idxChild]){
                    task->buffers[idxRealNeig].handle = child[idxChild]->handle;
                    task->buffers[idxRealNeig].mode = STARPU_R;
                    ++idxRealNeig;
                    bitDescriptor |= bitNeigIndex;
                }
                bitNeigIndex <<= 1;
            }
            if(needToRecv){
                const int startIndex = getWorkingInterval((idxLevel+1) ,firstProcThatSend).min & 0x7;
                bitNeigIndex = 1 << startIndex;

                for(int idxChild = startIndex ; idxChild < 8 ; ++idxChild){
                    task->buffers[idxRealNeig].handle = recvCellsM2M[idxLevel * 8 + idxChild].handle;
                    task->buffers[idxRealNeig].mode = STARPU_R;
                    ++idxRealNeig;
                    bitDescriptor |= bitNeigIndex;
                    bitNeigIndex <<= 1;
                }

                firstProcThatSend = endProcThatSend - 1;
            }

            task->cl = &m2m_cl[idxRealNeig - 1];

            task->cl_arg = &bitsDescriptors[bitDescriptor + idxLevel * 256];
            task->cl_arg_size = sizeof(int);

            starpu_task_submit(task);


            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }

        FDEBUG( FDebug::Controller << "\tFinished (@M2M = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }


    /////////////////////////////////////////////////////////////////////////////
    // M2L
    /////////////////////////////////////////////////////////////////////////////

    /** M2L */
    void M2LTasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart M2L\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        FVector<MortonIndex> indexes[nbProcess];
        int nbToSend[OctreeHeight][nbProcess];
        memset(nbToSend, 0, sizeof(int) * OctreeHeight * nbProcess);

        FBoolArray* cellsNeedOther[OctreeHeight];
        memset( cellsNeedOther, 0, sizeof(FBoolArray*) * OctreeHeight);

        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const CellClass* neighbors[189];

            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                if(idProcess != 0
                        && getWorkingInterval(idxLevel, idProcess).max <= getWorkingInterval(idxLevel, idProcess - 1).max){
                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    continue;
                }

                FVector< CellClass* > cellsToSend[nbProcess];
                cellsNeedOther[idxLevel] = new FBoolArray(numberOfLeafs);

                while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).min){
                    octreeIterator.moveRight();
                }

                int idLeaf = 0;
                do{
                    // Add exec
                    const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);

                    if(counter){
                        struct starpu_task* const task = starpu_task_create();

                        task->buffers[0].handle = octreeIterator.getCurrentCell()->handle;
                        task->buffers[0].mode = STARPU_RW;

                        for(int idxNeigh = 0 ; idxNeigh < counter ; ++idxNeigh){
                            task->buffers[idxNeigh + 1].handle = neighbors[idxNeigh]->handle;
                            task->buffers[idxNeigh + 1].mode = STARPU_R;
                        }

                        task->cl = &m2l_cl[counter];

                        task->cl_arg = &m2lNeighbors[idxLevel][counter];
                        task->cl_arg_size = sizeof(M2lProperties);

                        starpu_task_submit(task);
                    }

                    // Look if we need to send
                    bool alreadySent[nbProcess];
                    memset(alreadySent, false, sizeof(bool) * nbProcess);

                    MortonIndex neighborsIndexes[189];
                    const int counterDistant = getDistantNeighbors(octreeIterator.getCurrentGlobalCoordinate(),idxLevel,neighborsIndexes);
                    bool needOther = false;

                    // Test each negibors to know which one do not belong to us
                    for(int idxNeigh = 0 ; idxNeigh < counterDistant ; ++idxNeigh){
                        if(neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , idProcess).min
                                || getWorkingInterval(idxLevel , idProcess).max < neighborsIndexes[idxNeigh]){
                            int procToReceive = idProcess;
                            while( 0 != procToReceive && neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , procToReceive).min ){
                                --procToReceive;
                            }
                            while( procToReceive != nbProcess -1 && getWorkingInterval(idxLevel , procToReceive).max < neighborsIndexes[idxNeigh]){
                                ++procToReceive;
                            }
                            // Maybe already sent to that proc?
                            if( !alreadySent[procToReceive]
                                && getWorkingInterval(idxLevel , procToReceive).min <= neighborsIndexes[idxNeigh]
                                && neighborsIndexes[idxNeigh] <= getWorkingInterval(idxLevel , procToReceive).max){

                                alreadySent[procToReceive] = true;
                                ++nbToSend[idxLevel][procToReceive];
                                indexes[procToReceive].push(octreeIterator.getCurrentGlobalIndex());

                                needOther = true;

                                starpu_mpi_isend_detached(octreeIterator.getCurrentCell()->handle ,
                                    procToReceive, (((octreeIterator.getCurrentGlobalIndex() << 5) | idxLevel ) << 5) | FMpi::TagFmmM2L, MPI_COMM_WORLD , 0, 0 );
                            }
                        }
                    }

                    if( needOther ){
                        cellsNeedOther[idxLevel]->set(idLeaf, true);
                    }
                    ++idLeaf;

                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }
        }

        int allSendRecv[nbProcess][OctreeHeight][nbProcess];
        MPI_Allgather( nbToSend, nbProcess * OctreeHeight, MPI_INT, allSendRecv, nbProcess * OctreeHeight, MPI_INT, MPI_COMM_WORLD);

        MPI_Request requests[nbProcess * 2];
        int iterRequest = 0;

        MortonIndex* indexesFromOther[nbProcess];
        memset( indexesFromOther, 0, nbProcess * sizeof(MortonIndex*) );

        int totalToRecv = 0;

        for( int idxProc = 0 ; idxProc < nbProcess ; ++idxProc ){
            if( indexes[idxProc].getSize() ){
                MPI_Isend( &indexes[idxProc][0], indexes[idxProc].getSize() , MPI_LONG_LONG ,
                          idxProc, FMpi::TagFmmM2L, MPI_COMM_WORLD, &requests[iterRequest++]);
            }
            int hasToReceive = 0;
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                hasToReceive += allSendRecv[idxProc][idxLevel][idProcess];
            }
            if(hasToReceive){
                indexesFromOther[idxProc] = new MortonIndex[hasToReceive];
                MPI_Irecv(indexesFromOther[idxProc], hasToReceive, MPI_LONG_LONG,
                        idxProc, FMpi::TagFmmM2L, MPI_COMM_WORLD, &requests[iterRequest++]);
            }
            totalToRecv += hasToReceive;
        }

        MPI_Waitall( iterRequest, requests, MPI_STATUSES_IGNORE);

        int idxCellsBuffer = 0;
        cellsFromOtherM2L = new CellClass[totalToRecv];

        for( int idxProc = 0 ; idxProc < nbProcess ; ++idxProc ){
            if(idxProc != idProcess){
                int idxMortonBuffer = 0;
                for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){

                    for(int idxMorton = 0 ; idxMorton < allSendRecv[idxProc][idxLevel][idProcess] ; ++idxMorton){
                        cellsFromOtherM2L[idxCellsBuffer].initHandle();

                        starpu_mpi_irecv_detached(cellsFromOtherM2L[idxCellsBuffer].handle,
                            idxProc, (((indexesFromOther[idxProc][idxMortonBuffer]  << 5) | idxLevel ) << 5) | FMpi::TagFmmM2L, MPI_COMM_WORLD , 0, 0);

                        tempTreeM2L[idxLevel][ indexesFromOther[idxProc][idxMortonBuffer] ] = &cellsFromOtherM2L[idxCellsBuffer];

                        ++idxMortonBuffer;
                        ++idxCellsBuffer;
                    }
                }
            }
            delete[] indexesFromOther[idxProc];
        }

        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const CellClass* neighbors[189];

            int nbToSend[OctreeHeight][nbProcess];
            memset(nbToSend, 0, sizeof(int) * OctreeHeight * nbProcess);

            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                if(idProcess != 0
                        && getWorkingInterval(idxLevel, idProcess).max <= getWorkingInterval(idxLevel, idProcess - 1).max){

                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    continue;
                }

                while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).min){
                    octreeIterator.moveRight();
                }

                int realCellId = 0;

                do{
                    if( cellsNeedOther[idxLevel]->get(realCellId++) ){
                        // Add exec
                        MortonIndex neighborsIndex[189];
                        const int counterNeighbors = getDistantNeighbors(octreeIterator.getCurrentGlobalCoordinate(), idxLevel, neighborsIndex);

                        int counter = 0;
                        // does we receive this index from someone?
                        for(int idxNeig = 0 ;idxNeig < counterNeighbors ; ++idxNeig){
                            typename std::map<MortonIndex,CellClass*>::iterator it = tempTreeM2L[idxLevel].find(neighborsIndex[idxNeig]);
                            if(it != tempTreeM2L[idxLevel].end() ){
                                neighbors[counter++] = it->second;
                            }
                        }

                        if(counter){
                            struct starpu_task* const task = starpu_task_create();

                            task->buffers[0].handle = octreeIterator.getCurrentCell()->handle;
                            task->buffers[0].mode = STARPU_RW;

                            for(int idxNeigh = 0 ; idxNeigh < counter ; ++idxNeigh){
                                task->buffers[idxNeigh + 1].handle = neighbors[idxNeigh]->handle;
                                task->buffers[idxNeigh + 1].mode = STARPU_R;
                            }

                            task->cl = &m2l_cl[counter];

                            task->cl_arg = &m2lNeighbors[idxLevel][counter];
                            task->cl_arg_size = sizeof(M2lProperties);

                            starpu_task_submit(task);
                        }
                    }
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                delete cellsNeedOther[idxLevel];
                cellsNeedOther[idxLevel] = 0;
            }
        }

        FDEBUG( FDebug::Controller << "\tFinished (@M2L = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // L2L
    /////////////////////////////////////////////////////////////////////////////

    /** L2L */
    void L2LTasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart L2L\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        for(int idxLevel = 2 ; idxLevel < OctreeHeight - 1 ; ++idxLevel ){
            if(idProcess != 0
                    && getWorkingInterval((idxLevel+1) , idProcess).max <= getWorkingInterval((idxLevel+1) , idProcess - 1).max){

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                continue;
            }

            typename OctreeClass::Iterator previous = octreeIterator;
            while(octreeIterator.getCurrentGlobalIndex() < getWorkingInterval(idxLevel , idProcess).min){
                previous = octreeIterator;
                octreeIterator.moveRight();
            }

            bool needToRecv = false;
            if(idProcess != 0
                    && (getWorkingInterval((idxLevel + 1) , idProcess).min >> 3 ) <= (getWorkingInterval((idxLevel+1) , idProcess - 1).max >> 3 ) ){
                needToRecv = true;

                starpu_mpi_irecv_detached(previous.getCurrentCell()->handle,
                                          MPI_ANY_SOURCE, (((octreeIterator.getCurrentGlobalIndex()  << 5) | idxLevel) << 5) | FMpi::TagFmmL2L , MPI_COMM_WORLD , 0, 0);
                previous = octreeIterator;
            }

            do{
                //kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                struct starpu_task* const task = starpu_task_create();

                task->buffers[0].handle = octreeIterator.getCurrentCell()->handle;
                task->buffers[0].mode = STARPU_R;

                int idxRealNeig   = 1;
                int bitNeigIndex  = 1;
                int bitDescriptor = 0;

                CellClass*const*const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if(child[idxChild]){
                        task->buffers[idxRealNeig].handle = child[idxChild]->handle;
                        task->buffers[idxRealNeig].mode = STARPU_RW;
                        ++idxRealNeig;
                        bitDescriptor |= bitNeigIndex;
                    }
                    bitNeigIndex <<= 1;
                }

                task->cl = &l2l_cl[idxRealNeig - 1];

                task->cl_arg = &bitsDescriptors[bitDescriptor + idxLevel * 256];
                task->cl_arg_size = sizeof(int);

                starpu_task_submit(task);
            } while(octreeIterator.moveRight());

            if(idProcess != nbProcess - 1){
                int firstProcThatRecv = idProcess + 1;
                while( firstProcThatRecv < nbProcess &&
                      getWorkingInterval((idxLevel + 1) , firstProcThatRecv).max <= getWorkingInterval((idxLevel+1) , idProcess).max){
                    ++firstProcThatRecv;
                }

                int endProcThatRecv = firstProcThatRecv;
                while( endProcThatRecv < nbProcess &&
                      (getWorkingInterval((idxLevel + 1) , endProcThatRecv).min >> 3) <= (getWorkingInterval((idxLevel+1) , idProcess).max >> 3) ){
                    ++endProcThatRecv;
                }

                if(firstProcThatRecv != endProcThatRecv){
                    for(int idxProc = firstProcThatRecv ; idxProc < endProcThatRecv ; ++idxProc ){
                        starpu_mpi_isend_detached(octreeIterator.getCurrentCell()->handle ,
                                                  idxProc, (((octreeIterator.getCurrentGlobalIndex()  << 5) | idxLevel) << 5) | FMpi::TagFmmL2L, MPI_COMM_WORLD , 0 , 0);
                    }
                }
            }

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }

        FDEBUG( FDebug::Controller << "\tFinished (@L2L = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // L2P
    /////////////////////////////////////////////////////////////////////////////

    /** L2P */
    void L2PTasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart L2P\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            //kernels->L2P( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
            starpu_insert_task(&l2p_cl, STARPU_R, octreeIterator.getCurrentCell()->handle, STARPU_RW, octreeIterator.getCurrentLeaf()->getTargets()->handle, 0);
        } while(octreeIterator.moveRight());

        FDEBUG( FDebug::Controller << "\tFinished (@L2P = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Callback
    /////////////////////////////////////////////////////////////////////////////

    static KernelClass** globalKernels;

    static void p2m_cpu(void *descr[], void *)
    {
        CellClass* const currentCell = (CellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);

        starpu_vector_interface_t* const data = (starpu_vector_interface_t*)descr[1];
        StarVector<ParticleClass> particles((ParticleClass*)STARPU_VECTOR_GET_PTR(data), STARPU_VECTOR_GET_NX(data));

        globalKernels[starpu_worker_get_id()]->P2M( currentCell , &particles );
    }

    static void m2m_cpu(void *descr[], void *cl_arg)
    {
        int descriptor = (*static_cast<int*>(cl_arg));
        const int level = descriptor & 0xFF;
        descriptor >>= 8;

        CellClass* const currentCell = (CellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);

        const CellClass* child[8];
        memset(child, 0, sizeof(CellClass*)*8);

        int bufferIndex = 1;
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(descriptor & 1){
                child[idxChild] = ((const CellClass*)STARPU_VARIABLE_GET_PTR(descr[bufferIndex++]));
            }
            descriptor >>= 1;
        }

        globalKernels[starpu_worker_get_id()]->M2M( currentCell , child , level);
    }

    static void m2l_cpu(void *descr[], void *cl_arg)
    {
        const M2lProperties*const properties = static_cast<M2lProperties*>(cl_arg);
        CellClass* const currentCell = (CellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);

        const CellClass* neighbor[189];
        memset(neighbor, 0, 189 * sizeof(CellClass*));

        for(int idxNeig = 0 ; idxNeig < properties->bufferSize ; ++idxNeig){
            neighbor[idxNeig] =  ((CellClass*)STARPU_VARIABLE_GET_PTR(descr[idxNeig + 1]));
        }

        globalKernels[starpu_worker_get_id()]->M2L( currentCell , neighbor, properties->bufferSize, properties->level );

    }


    static void l2l_cpu(void *descr[], void * cl_arg)
    {
        int descriptor = (*static_cast<int*>(cl_arg));
        const int level = descriptor & 0xFF;
        descriptor >>= 8;

        const CellClass* const currentCell = (const CellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);

        CellClass* child[8];
        memset(child, 0, sizeof(CellClass*)*8);

        int bufferIndex = 1;
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(descriptor & 1){
                child[idxChild] = ((CellClass*)STARPU_VARIABLE_GET_PTR(descr[bufferIndex++]));
            }
            descriptor >>= 1;
        }

        globalKernels[starpu_worker_get_id()]->L2L( currentCell , child , level);
    }

    static void l2p_cpu(void *descr[], void *)
    {
        const CellClass* const currentCell = (const CellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);

        starpu_vector_interface_t* const data = (starpu_vector_interface_t*)descr[1];
        StarVector<ParticleClass> particles((ParticleClass*)STARPU_VECTOR_GET_PTR(data), STARPU_VECTOR_GET_NX(data));

        globalKernels[starpu_worker_get_id()]->L2P( currentCell , &particles );
    }

    static void p2p_cpu(void *descr[], void *cl_arg) {

        const int nbNeighbors = (*static_cast<int*>(cl_arg));

        starpu_vector_interface_t* const data0 = (starpu_vector_interface_t*)descr[0];
        StarVector<ParticleClass> currentContainer((ParticleClass*)STARPU_VECTOR_GET_PTR(data0), STARPU_VECTOR_GET_NX(data0));

        StarVector<ParticleClass>* neighors[26];
        memset(neighors, 0, 26 * sizeof(ContainerClass*) );

        for(int idxNeig = 0 ; idxNeig < nbNeighbors; ++idxNeig){
            starpu_vector_interface_t* const data = (starpu_vector_interface_t*)descr[idxNeig + 1];
            neighors[idxNeig] = new StarVector<ParticleClass>((ParticleClass*)STARPU_VECTOR_GET_PTR(data), STARPU_VECTOR_GET_NX(data));
        }

        globalKernels[starpu_worker_get_id()]->P2P( &currentContainer, &currentContainer , neighors, nbNeighbors );

        for(int idxNeig = 0 ; idxNeig < nbNeighbors; ++idxNeig){
            delete neighors[idxNeig];
        }
    }


    int getNeighborsIndexes(const MortonIndex centerIndex, const long limite, MortonIndex indexes[26]) const{
        FTreeCoordinate center;
        center.setPositionFromMorton(centerIndex, OctreeHeight - 1);
        int idxNeig = 0;
        // We test all cells around
        for(long idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(center.getX() + idxX,0l, limite)) continue;

            for(long idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(center.getY() + idxY,0l, limite)) continue;

                for(long idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(center.getZ() + idxZ,0l, limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                        indexes[idxNeig++] = other.getMortonIndex(this->OctreeHeight - 1);
                    }
                }
            }
        }
        return idxNeig;
    }

    int getDistantNeighbors(const FTreeCoordinate& workingCell,const int inLevel, MortonIndex inNeighbors[189]) const{

        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(workingCell.getX()>>1,workingCell.getY()>>1,workingCell.getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const long limite = FMath::pow(2,inLevel-1);

        int idxNeighbors = 0;
        // We test all cells around
        for(long idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(parentCell.getX() + idxX,0l,limite)) continue;

            for(long idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(parentCell.getY() + idxY,0l,limite)) continue;

                for(long idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(parentCell.getZ() + idxZ,0l,limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        const MortonIndex mortonOther = other.getMortonIndex(inLevel-1);

                        // For each child
                        for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                            const FTreeCoordinate potentialNeighbor((other.getX()<<1) | (idxCousin>>2 & 1),
                                                                   (other.getY()<<1) | (idxCousin>>1 & 1),
                                                                   (other.getZ()<<1) | (idxCousin&1));
                            // Test if it is a direct neighbor
                            if(FMath::Abs(workingCell.getX() - potentialNeighbor.getX()) > 1 ||
                                    FMath::Abs(workingCell.getY() - potentialNeighbor.getY()) > 1 ||
                                    FMath::Abs(workingCell.getZ() - potentialNeighbor.getZ()) > 1){
                                // add to neighbors
                                inNeighbors[idxNeighbors++] = (mortonOther << 3) | idxCousin;
                            }
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }
};


#endif //FFMMALGORITHM_HPP

// [--LICENSE--]
