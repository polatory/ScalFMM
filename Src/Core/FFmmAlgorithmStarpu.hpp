#ifndef FFMMALGORITHMSTARPU_HPP
#define FFMMALGORITHMSTARPU_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FNoCopyable.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"

#include <limits>
#include <starpu.h>


struct StarHandle : public FNoCopyable, public FNoAssignement {
    starpu_data_handle_t handle;

    StarHandle(){
        memset(&handle, 0, sizeof(starpu_data_handle_t));
    }

    ~StarHandle(){
        if( handle ){
            starpu_data_unregister(handle);
        }
    }

    template <class ObjectType>
    void registerVariable(ObjectType*const inData){
        starpu_variable_data_register(&handle, 0, (uintptr_t)inData, sizeof(ObjectType));
    }

    template <class ObjectType>
    void registerVariable(ObjectType*const inData, const int sizeOfBlock){
        starpu_variable_data_register(&handle, 0, (uintptr_t)inData, sizeOfBlock);
    }

    template <class ObjectType>
    void registerVector(ObjectType*const inData, const int inSize){
        starpu_vector_data_register(&handle, 0, (uintptr_t)inData, inSize, sizeof(ObjectType));
    }

    void unregisterData(){
        if( handle ){
            starpu_data_unregister(handle);
            memset(&handle, 0, sizeof(starpu_data_handle_t));
        }
    }
};


class AbstractStarCell {
public:
    StarHandle handleUp;
    StarHandle handleDown;

    virtual void initHandle() = 0;
};


template < class ElementClass >
class StarVector : public FVector<ElementClass> {
public:
    StarHandle handle;

    void initHandle() {
        handle.registerVector( FVector<ElementClass>::data(), FVector<ElementClass>::getSize());
    }
};


template<class T>
class DataVector {
protected:
    T* array;
    int size;

    DataVector& operator=(const DataVector&){ return *this;}
    DataVector(const DataVector&){}

public:
    DataVector(T*const inData, const unsigned int inSize) : array(inData), size(inSize) {
    }

    virtual ~DataVector(){
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


    class BasicIterator {
    protected:
        DataVector* const vector;  /**< the vector to work on*/
        int index;              /**< the current node*/

    public:
        /** Empty destructor */
        virtual ~BasicIterator(){}

        /** Constructor need a vector */
        explicit BasicIterator(DataVector<T>& inVector) : vector(&inVector), index(0){}

        /** Go to next vector element */
        void gotoNext() {
            ++index;
        }

        /** is it over
          * @return true if we are over the vector
          */
        bool hasNotFinished() const {
            return index < vector->size;
        }

        /** Get current data */
        T& data(){
            return vector->array[index];
        }

        /** Get current data */
        const T& data() const{
            return vector->array[index];
        }

        /** Set the data */
        void setData(const T& inData){
            vector->array[index] = inData;
        }

        /** Remove current data
          * It will move all the data after to their previous position
          */
        void remove(){
            if( hasNotFinished() ){
                for(int idxMove = index + 1; idxMove < vector->index ; ++idxMove){
                    vector->array[idxMove - 1] = vector->array[idxMove];
                }
                vector->index -= 1;
            }
        }

    };
    friend class BasicIterator;

    /** This class is a basic const iterator
      * it uses a const vector to work on
      */
    class ConstBasicIterator {
    protected:
        const DataVector* const vector;  /**< the vector to work on*/
        int index;                    /**< the current node*/

    public:
        /** Empty destructor */
        virtual ~ConstBasicIterator(){}

        /** Constructor need a vector */
        explicit ConstBasicIterator(const DataVector<T>& inVector) : vector(&inVector), index(0){}

        /** Go to next vector element */
        void gotoNext(){
            ++index;
        }

        /** is it over
          * @return true if we are over the vector
          */
        bool hasNotFinished() const{
            return index < vector->size;
        }

        /** Get current data */
        const T& data() const{
            return vector->array[index];
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
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass, class CellType>
class FFmmAlgorithmStarpu : protected FAssertable{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;

    //////////////////////////////////////////////////////////////////
    // Codelets
    //////////////////////////////////////////////////////////////////

    static unsigned int EmptyValue;
    static StarHandle EmptyHandle;


    starpu_codelet p2m_cl;
    starpu_codelet p2p_cl;
    starpu_codelet m2m_cl;
    starpu_codelet m2l_cl;
    starpu_codelet l2l_cl;
    starpu_codelet l2p_cl;

    void initCodelets(){
        memset(&p2m_cl, 0, sizeof(p2m_cl));
        p2m_cl.where = STARPU_CPU;
        p2m_cl.cpu_funcs[0] = p2m_cpu;
        p2m_cl.nbuffers = 2;

        memset(&p2p_cl, 0, sizeof(p2p_cl));
        p2p_cl.where = STARPU_CPU;
        p2p_cl.cpu_funcs[0] = p2p_cpu;
        p2p_cl.nbuffers = 28;

        memset(&m2l_cl, 0, sizeof(m2l_cl));
        m2l_cl.where = STARPU_CPU;
        m2l_cl.cpu_funcs[0] = m2l_cpu;
        m2l_cl.nbuffers = 344;

        memset(&l2p_cl, 0, sizeof(l2p_cl));
        l2p_cl.where = STARPU_CPU;
        l2p_cl.cpu_funcs[0] = l2p_cpu;
        l2p_cl.nbuffers = 2;

        memset(&m2m_cl, 0, sizeof(m2m_cl));
        memset(&l2l_cl, 0, sizeof(l2l_cl));

        m2m_cl.where = STARPU_CPU;
        m2m_cl.cpu_funcs[0] = m2m_cpu;
        m2m_cl.nbuffers = 9;

        l2l_cl.where = STARPU_CPU;
        l2l_cl.cpu_funcs[0] = l2l_cpu;
        l2l_cl.nbuffers = 9;
    }

    void releaseCodelets(){
    }

    //////////////////////////////////////////////////////////////////
    // Manage args
    //////////////////////////////////////////////////////////////////

    int* argsLevels;

    void initArgs(){
        argsLevels = new int[OctreeHeight];
        for( int idx = 0 ;idx  < OctreeHeight ; ++idx){
            argsLevels[idx] = idx;
        }
    }

    void releaseArgs(){
        delete[] argsLevels;
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


public:	
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmStarpu(OctreeClass* const inTree, KernelClass* const inKernels)
                      : tree(inTree) , kernels(inKernels), OctreeHeight(tree->getHeight()) {
        FDEBUG(FDebug::Controller << "FFmmAlgorithmStarpu\n");

        starpu_init(NULL);

        EmptyHandle.registerVariable( &EmptyValue );

        initCodelets();
        initArgs();
        initHandles();
        initKernels();
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmStarpu(){
        EmptyHandle.unregisterData();

        releaseCodelets();
        releaseArgs();
        releaseHandles();
        releaseKernels();

        starpu_shutdown();
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        P2M_P2P_Tasks();

        M2M_M2L_Tasks();

        L2LTasks();
        L2PTasks();

        FDEBUG(FDebug::Controller << "Waiting...\n");
        starpu_task_wait_for_all();

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void P2M_P2P_Tasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart P2M && P2P\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        ContainerClass* neighbors[27];
        const int heightMinusOne = OctreeHeight - 1;

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            // P2M
            {
                //kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
                starpu_insert_task( &p2m_cl, STARPU_RW, octreeIterator.getCurrentCell()->handleUp.handle,
                                    STARPU_R, octreeIterator.getCurrentLeaf()->getSrc()->handle.handle, 0);
            }
            // P2P
            {
                /*const int counter =*/ tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), heightMinusOne);
                //kernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                struct starpu_task* const task = starpu_task_create();

                task->handles[0] = octreeIterator.getCurrentLeaf()->getTargets()->handle.handle;
                //task->handles[0].mode = STARPU_RW;

                for(int idxCounterNeigh = 0 ; idxCounterNeigh < 27 ; ++idxCounterNeigh){
                    if( neighbors[idxCounterNeigh] ){
                        task->handles[idxCounterNeigh + 1] = neighbors[idxCounterNeigh]->handle.handle;
                        //task->handles[idxCounterNeigh + 1].mode   = STARPU_RW;
                    }
                    else {
                        task->handles[idxCounterNeigh + 1] = EmptyHandle.handle;
                        //task->handles[idxCounterNeigh + 1].mode   = STARPU_R;
                    }
                }

                task->cl = &p2p_cl;

                task->cl_arg = const_cast<FTreeCoordinate*>(&octreeIterator.getCurrentGlobalCoordinate());
                task->cl_arg_size = sizeof(FTreeCoordinate);

                task->priority = starpu_sched_get_min_priority();

                starpu_task_submit(task);

            }
        } while(octreeIterator.moveRight());

        FDEBUG( FDebug::Controller << "\tFinished (@P2M P2P = "  << counterTime.tacAndElapsed() << "s)\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // M2M M2L
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void M2M_M2L_Tasks(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart M2M M2L\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const CellClass* neighbors[343];

        // M2L at leaf level
        {
            const int idxLevel = OctreeHeight - 1;
            do{
                 const int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);

                if(counter){
                    struct starpu_task* const task = starpu_task_create();

                    task->handles[0] = octreeIterator.getCurrentCell()->handleDown.handle;
                    //task->handles[0].mode = STARPU_RW;

                    for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                        if( neighbors[idxNeigh] ){
                            task->handles[idxNeigh+1] = neighbors[idxNeigh]->handleUp.handle;
                            //task->handles[idxNeigh+1].mode = STARPU_R;
                        }
                        else {
                            task->handles[idxNeigh+1] = EmptyHandle.handle;
                            //task->handles[idxNeigh+1].mode = STARPU_R;
                        }
                    }

                    task->cl = &m2l_cl;

                    task->cl_arg = &argsLevels[idxLevel];
                    task->cl_arg_size = sizeof(int);

                    starpu_task_submit(task);
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            // M2M
            do{
                //kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                struct starpu_task* const task = starpu_task_create();

                task->handles[0] = octreeIterator.getCurrentCell()->handleUp.handle;
                //task->handles[0].mode = STARPU_RW;

                CellClass*const*const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if(child[idxChild]){
                        task->handles[idxChild+1] = child[idxChild]->handleUp.handle;
                        //task->handles[idxChild+1].mode = STARPU_R;
                    }
                    else{
                        task->handles[idxChild+1] = EmptyHandle.handle;
                        //task->handles[idxChild+1].mode = STARPU_R;
                    }
                }

                task->cl = &m2m_cl;

                task->cl_arg = &argsLevels[idxLevel];
                task->cl_arg_size = sizeof(int);

                starpu_task_submit(task);
            } while(octreeIterator.moveRight());

            octreeIterator = avoidGotoLeftIterator;

            // M2L
            do{
                const int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);

                if(counter){
                    struct starpu_task* const task = starpu_task_create();

                    task->handles[0] = octreeIterator.getCurrentCell()->handleDown.handle;
                    //task->handles[0].mode = STARPU_RW;

                    for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                        if( neighbors[idxNeigh] ){
                            task->handles[idxNeigh+1] = neighbors[idxNeigh]->handleUp.handle;
                            //task->handles[idxNeigh+1].mode = STARPU_R;
                        }
                        else {
                            task->handles[idxNeigh+1] = EmptyHandle.handle;
                            //task->handles[idxNeigh+1].mode = STARPU_R;
                        }
                    }

                    task->cl = &m2l_cl;

                    task->cl_arg = &argsLevels[idxLevel];
                    task->cl_arg_size = sizeof(int);

                    starpu_task_submit(task);
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;

        }

        FDEBUG( FDebug::Controller << "\tFinished (@M2M = "  << counterTime.tacAndElapsed() << "s)\n" );
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
            do{
                //kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                struct starpu_task* const task = starpu_task_create();

                task->handles[0] = octreeIterator.getCurrentCell()->handleDown.handle;
                //task->handles[0].mode = STARPU_R;

                CellClass*const*const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if(child[idxChild]){
                        task->handles[idxChild+1] = child[idxChild]->handleDown.handle;
                        //task->handles[idxChild+1].mode = STARPU_RW;
                    }
                    else{
                        task->handles[idxChild+1] = EmptyHandle.handle;
                        //task->handles[idxChild+1].mode = STARPU_R;
                    }
                }

                task->cl = &l2l_cl;

                task->cl_arg = &argsLevels[idxLevel];
                task->cl_arg_size = sizeof(int);

                starpu_task_submit(task);

            } while(octreeIterator.moveRight());

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
            starpu_insert_task(&l2p_cl, STARPU_R, octreeIterator.getCurrentCell()->handleDown.handle,
                               STARPU_RW, octreeIterator.getCurrentLeaf()->getTargets()->handle.handle, 0);
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
        CellType* const currentCell = (CellType*)STARPU_VARIABLE_GET_PTR(descr[0]);

        DataVector<ParticleClass> particles((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[1]), STARPU_VECTOR_GET_NX(descr[1]));

        globalKernels[starpu_worker_get_id()]->P2M( currentCell , &particles );
    }

    static void m2m_cpu(void *descr[], void *cl_arg)
    {
        const CellType* child[8];
        memset(child, 0, sizeof(CellType*)*8);

        const int level = (*static_cast<int*>(cl_arg));
        CellType* const currentCell = (CellType*)STARPU_VARIABLE_GET_PTR(descr[0]);

        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            const unsigned int empty = *((const unsigned int*)STARPU_VARIABLE_GET_PTR(descr[idxChild+1]));
            if(empty != EmptyValue){
                child[idxChild] = ((const CellType*)STARPU_VARIABLE_GET_PTR(descr[idxChild+1]));
            }
        }

        globalKernels[starpu_worker_get_id()]->M2M( currentCell , child , level);
    }

    static void m2l_cpu(void *descr[], void *cl_arg)
    {return;
        const int level = (*static_cast<int*>(cl_arg));
        CellType* const currentCell = (CellType*)STARPU_VARIABLE_GET_PTR(descr[0]);

        const CellType* neighbor[343];
        memset(neighbor, 0, 343 * sizeof(CellType*));

        int counter = 0;
        for(int idxNeig = 0 ; idxNeig < 343 ; ++idxNeig){
            const unsigned int empty = *((const unsigned int*)STARPU_VARIABLE_GET_PTR(descr[idxNeig+1]));
            if(empty != EmptyValue){
                neighbor[idxNeig] = ((const CellType*)STARPU_VARIABLE_GET_PTR(descr[idxNeig+1]));
                ++counter;
            }
        }

        globalKernels[starpu_worker_get_id()]->M2L( currentCell , neighbor, counter, level );

    }


    static void l2l_cpu(void *descr[], void * cl_arg)
    {
        CellType* child[8];
        memset(child, 0, sizeof(CellType*)*8);

        const int level = (*static_cast<int*>(cl_arg));
        CellType* const currentCell = (CellType*)STARPU_VARIABLE_GET_PTR(descr[0]);

        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            const unsigned int empty = *((const unsigned int*)STARPU_VARIABLE_GET_PTR(descr[idxChild+1]));
            if(empty != EmptyValue){
                child[idxChild] = ((CellType*)STARPU_VARIABLE_GET_PTR(descr[idxChild+1]));
            }

        }

        globalKernels[starpu_worker_get_id()]->L2L( currentCell , child , level);
    }

    static void l2p_cpu(void *descr[], void *)
    {
        const CellType* const currentCell = (const CellType*)STARPU_VARIABLE_GET_PTR(descr[0]);

        DataVector<ParticleClass> particles((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[1]), STARPU_VECTOR_GET_NX(descr[1]));

        globalKernels[starpu_worker_get_id()]->L2P( currentCell , &particles );
    }

    static void p2p_cpu(void *descr[], void* cl_arg) {
        DataVector<ParticleClass> currentContainer((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[0]), STARPU_VECTOR_GET_NX(descr[0]));
        const FTreeCoordinate& coordinate = *(static_cast<const FTreeCoordinate*>(cl_arg));

        DataVector<ParticleClass>* neighors[27];
        memset(neighors, 0, 27 * sizeof(DataVector<ParticleClass>*) );

        int counter = 0;
        for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
            const unsigned int empty = *((const unsigned int*)STARPU_VARIABLE_GET_PTR(descr[idxNeig+1]));
            if(empty != EmptyValue){
                neighors[idxNeig] = new DataVector<ParticleClass>((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[idxNeig+1]), STARPU_VECTOR_GET_NX(descr[idxNeig+1]));
                ++counter;
            }
        }

        globalKernels[starpu_worker_get_id()]->P2P(coordinate, &currentContainer, &currentContainer , neighors, counter);

        for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
            delete neighors[idxNeig];
        }
    }

};

template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass, class CellType>
unsigned int FFmmAlgorithmStarpu<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass,CellType>::EmptyValue = 0xFFFFFFFF;
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass, class CellType>
StarHandle FFmmAlgorithmStarpu<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass,CellType>::EmptyHandle;

#endif //FFMMALGORITHMSTARPU_HPP

// [--LICENSE--]
