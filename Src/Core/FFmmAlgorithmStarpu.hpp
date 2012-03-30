// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
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


////////////////////////////////////////////////////////
// Utils
////////////////////////////////////////////////////////

/** This class is a c++ wrapper to starpu data handle
  * @warning Delete is not working yet
  */
struct StarHandle : public FNoCopyable, public FNoAssignement {
    /** The stapu handle */
    starpu_data_handle_t handle;

    /** Init with 0 */
    StarHandle(){
        memset(&handle, 0, sizeof(starpu_data_handle_t));
    }

    /** Release the handle */
    ~StarHandle(){
        if( handle != starpu_data_handle_t(0) ){
            starpu_data_unregister(handle);
        }
    }

    /** Register a data from its size as a variable */
    template <class ObjectType>
    void registerVariable(ObjectType*const inData){
        starpu_variable_data_register(&handle, 0, (uintptr_t)inData, sizeof(ObjectType));
    }

    /** Register a data from a given size as a variable */
    template <class ObjectType>
    void registerVariable(ObjectType*const inData, const int sizeOfBlock){
        starpu_variable_data_register(&handle, 0, (uintptr_t)inData, sizeOfBlock);
    }

    /** Register a data vector from its data-type size */
    template <class ObjectType>
    void registerVector(ObjectType*const inData, const int inSize){
        starpu_vector_data_register(&handle, 0, (uintptr_t)inData, inSize, sizeof(ObjectType));
    }

    /** Release data */
    void unregisterData(){
        if( handle != ((void *)0) ){
            starpu_data_unregister(handle);
            memset(&handle, 0, sizeof(starpu_data_handle_t));
        }
    }
};


/** This has to be used to make a cell
  * starpu enabled
  */
template <class CellClass>
class FStarCell : public CellClass{
public:
    /** The handle to register the data */
    StarHandle handleUp;
    StarHandle handleDown;
    /** Called by fmm starpu to register data */
    void initHandle(){
        handleUp.registerVariable( static_cast<CellClass*>(this) );
        handleDown.registerVariable( static_cast<CellClass*>(this) );
    }
};

/** This has to be used to make a vector
  * starpu enabled
  */
template < class ElementClass >
class StarVector : public FVector<ElementClass> {
public:
    /** The handle to register the data */
    StarHandle handle;
    /** Called by fmm starpu to register data */
    void initHandle() {
        handle.registerVector( FVector<ElementClass>::data(), FVector<ElementClass>::getSize());
    }
};

/** Wrapper vector to wrappe data from a real vector
  * without copying again the data
  * This vector do not release the data
  */
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

////////////////////////////////////////////////////////
// Core
////////////////////////////////////////////////////////


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
template<class OctreeClass, class ParticleClass, class CellClass, class RealCellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmStarpu : protected FAssertable{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;
    const bool putNameInTask;

    //////////////////////////////////////////////////////////////////
    // Codelets
    //////////////////////////////////////////////////////////////////

    // All the posible codelet
    starpu_codelet p2m_cl;
    starpu_codelet p2p_cl[28];
    starpu_codelet m2m_cl[8];
    starpu_codelet m2l_cl[189];
    starpu_codelet l2l_cl[8];
    starpu_codelet l2p_cl;

    starpu_perfmodel p2p_model;
    starpu_perfmodel p2m_model;
    starpu_perfmodel m2m_model;
    starpu_perfmodel m2l_model;
    starpu_perfmodel l2l_model;
    starpu_perfmodel l2p_model;

    // Init the codelet
    void initCodelets(){
        memset(&p2p_model, 0, sizeof(p2p_model));
        p2p_model.type = STARPU_HISTORY_BASED;
        p2p_model.symbol = "P2P";
        memset(&p2m_model, 0, sizeof(p2m_model));
        p2m_model.type = STARPU_HISTORY_BASED;
        p2m_model.symbol = "P2M";
        memset(&m2l_model, 0, sizeof(m2l_model));
        m2l_model.type = STARPU_HISTORY_BASED;
        m2l_model.symbol = "M2L";
        memset(&l2p_model, 0, sizeof(l2p_model));
        l2p_model.type = STARPU_HISTORY_BASED;
        l2p_model.symbol = "L2P";
        memset(&l2l_model, 0, sizeof(l2l_model));
        l2l_model.type = STARPU_HISTORY_BASED;
        l2l_model.symbol = "L2L";
        memset(&m2m_model, 0, sizeof(m2m_model));
        m2m_model.type = STARPU_HISTORY_BASED;
        m2m_model.symbol = "M2M";

        // P2M
        memset(&p2m_cl, 0, sizeof(p2m_cl));
        p2m_cl.where = STARPU_CPU;
        p2m_cl.cpu_funcs[0] = p2m_cpu;
        p2m_cl.nbuffers = 2;
        p2m_cl.modes[0] = STARPU_W;
        p2m_cl.modes[1] = STARPU_R;
        if(putNameInTask) p2m_cl.model = &p2m_model;

        // P2P
        memset(p2p_cl, 0, sizeof(starpu_codelet) * 28);
        for(int idxNeig = 0 ; idxNeig <= 27 ; ++idxNeig){
            p2p_cl[idxNeig].where = STARPU_CPU;
            p2p_cl[idxNeig].cpu_funcs[0] = p2p_cpu;
            p2p_cl[idxNeig].nbuffers = idxNeig + 1;

            if( putNameInTask ) p2p_cl[idxNeig].model = &p2p_model;

            for( int idxMode = 0 ; idxMode <= idxNeig ; ++idxMode){
                p2p_cl[idxNeig].modes[idxMode] = STARPU_RW;
            }
        }
        // M2L
        memset(m2l_cl, 0, sizeof(starpu_codelet) * 189);
        for(int idxNeig = 0 ; idxNeig < 189 ; ++idxNeig){
            m2l_cl[idxNeig].where = STARPU_CPU;
            m2l_cl[idxNeig].cpu_funcs[0] = m2l_cpu;
            m2l_cl[idxNeig].nbuffers = idxNeig + 2;

            if( putNameInTask ) m2l_cl[idxNeig].model = &m2l_model;

            m2l_cl[idxNeig].modes[0] = STARPU_RW;
            for( int idxMode = 0 ; idxMode <= idxNeig ; ++idxMode){
                m2l_cl[idxNeig].modes[idxMode+1] = STARPU_R;
            }
        }
        // L2P
        memset(&l2p_cl, 0, sizeof(l2p_cl));
        l2p_cl.where = STARPU_CPU;
        l2p_cl.cpu_funcs[0] = l2p_cpu;
        l2p_cl.nbuffers = 2;
        l2p_cl.modes[0] = STARPU_R;
        l2p_cl.modes[1] = STARPU_RW;
        if(putNameInTask)  l2p_cl.model = &l2p_model;

        // M2M & L2L
        memset(m2m_cl, 0, sizeof(starpu_codelet) * 8);
        memset(l2l_cl, 0, sizeof(starpu_codelet) * 8);
        for( int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            m2m_cl[idxChild].where = STARPU_CPU;
            m2m_cl[idxChild].cpu_funcs[0] = m2m_cpu;
            m2m_cl[idxChild].nbuffers = idxChild + 2;
            m2m_cl[idxChild].modes[0] = STARPU_W;
            if( putNameInTask)m2m_cl[idxChild].model = &m2m_model;

            l2l_cl[idxChild].where = STARPU_CPU;
            l2l_cl[idxChild].cpu_funcs[0] = l2l_cpu;
            l2l_cl[idxChild].nbuffers = idxChild + 2;
            l2l_cl[idxChild].modes[0] = STARPU_R;
            if( putNameInTask)l2l_cl[idxChild].model = &l2l_model;

            for( int idxMode = 0 ; idxMode <= idxChild ; ++idxMode){
                m2m_cl[idxChild].modes[idxMode+1] = STARPU_R;
                l2l_cl[idxChild].modes[idxMode+1] = STARPU_RW;
            }
        }
    }

    // Release the codelet
    void releaseCodelets(){
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
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
        // init leaf handle
        do{
            octreeIterator.getCurrentLeaf()->getSrc()->handle.unregisterData();
            if(octreeIterator.getCurrentLeaf()->getSrc() != octreeIterator.getCurrentLeaf()->getTargets()){
                octreeIterator.getCurrentLeaf()->getTargets()->handle.unregisterData();
            }
        } while(octreeIterator.moveRight());

        octreeIterator = avoidGotoLeftIterator;

        // init cells handle
        for(int idxLevel = OctreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            do{
                octreeIterator.getCurrentCell()->handleUp.unregisterData();
                octreeIterator.getCurrentCell()->handleDown.unregisterData();
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }
    }

    //////////////////////////////////////////////////////////////////
    // Init Kernels
    //////////////////////////////////////////////////////////////////

    // Init the fmm kernel (1 per thread)
    void initKernels(){
        globalKernels = new KernelClass*[starpu_worker_get_count()];
        memset(globalKernels, 0, sizeof(KernelClass*) * starpu_worker_get_count());

        for(unsigned int workerid = 0; workerid < starpu_worker_get_count(); ++workerid){
            if( starpu_worker_get_type(workerid) == STARPU_CPU_WORKER ){
                globalKernels[workerid] = new KernelClass(*kernels);
            }
        }
    }

    // Delete kernels
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
    FFmmAlgorithmStarpu(OctreeClass* const inTree, KernelClass* const inKernels, const bool inPutNameInTask = false)
        : tree(inTree) , kernels(inKernels), OctreeHeight(tree->getHeight()), putNameInTask(inPutNameInTask) {

        FDEBUG(FDebug::Controller << "FFmmAlgorithmStarpu\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmStarpu(){
    }

    /** Run starpu */
    void initStarpu(){
        // Run starpu
        starpu_init(NULL);
        FDEBUG(FDebug::Controller << "Init starpu, there are " << starpu_worker_get_count() << " workers\n");

        // Init
        initCodelets();
        initHandles();
        initKernels();
    }

    /** Release starpu */
    void releaseStarpu(){
        // Release stuff
        releaseCodelets();
        releaseHandles();
        releaseKernels();
        // Shutdown
        starpu_shutdown();
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        // Insert tasks
        P2M_P2P_Tasks();

        M2M_M2L_Tasks();

        L2LTasks();
        L2PTasks();
        // Wait
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

        // Neeed data
        ContainerClass* neighbors[27];
        const int heightMinusOne = OctreeHeight - 1;

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            // P2M
            {
                //kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
                starpu_insert_task( &p2m_cl, STARPU_W, octreeIterator.getCurrentCell()->handleUp.handle,
                                    STARPU_R, octreeIterator.getCurrentLeaf()->getSrc()->handle.handle, 0);
            }
            // P2P
            {
                const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), heightMinusOne);
                //kernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                struct starpu_task* const task = starpu_task_create();

                // handles 0 is the current leaf
                task->handles[0] = octreeIterator.getCurrentLeaf()->getTargets()->handle.handle;

                // then we insert neighbors with a mask system
                unsigned int mask = 0;
                int idxInsert = 1;
                for(int idxCounterNeigh = 0 ; idxCounterNeigh < 27 ; ++idxCounterNeigh){
                    if( neighbors[idxCounterNeigh] ){
                        task->handles[idxInsert++] = neighbors[idxCounterNeigh]->handle.handle;
                        mask = mask | (1 << idxCounterNeigh);
                    }
                }
                // Put the right codelet
                task->cl = &p2p_cl[counter];

                // Add buffer data
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                        STARPU_VALUE, &mask, sizeof(mask),
                        STARPU_VALUE, const_cast<FTreeCoordinate*>(&octreeIterator.getCurrentGlobalCoordinate()), sizeof(FTreeCoordinate),
                        0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                // Set priority to min
                task->priority = starpu_sched_get_min_priority();
                // submit task
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
        // Start from leaf level for M2L
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // Needed data
        const CellClass* neighbors[343];
        unsigned int mask_m2l[12];

        // M2L at leaf level
        {
            const int idxLevel = OctreeHeight - 1;
            // for all leafs
            do{
                // get interaction list
                const int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                // if not empty
                if(counter){
                    // create task
                    struct starpu_task* const task = starpu_task_create();
                    // buffer 0 is current leaf
                    task->handles[0] = octreeIterator.getCurrentCell()->handleDown.handle;

                    // insert other with a mask
                    memset(mask_m2l, 0, sizeof(unsigned int) * 12);
                    int idxInsert = 1;
                    for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                        if( neighbors[idxNeigh] ){
                            task->handles[idxInsert++] = neighbors[idxNeigh]->handleUp.handle;
                            mask_m2l[ idxNeigh >> 5 ] = mask_m2l[ idxNeigh >> 5 ] | (1 << (idxNeigh & 0x1F));
                        }
                    }
                    // put the right codelet
                    task->cl = &m2l_cl[counter-1];
                    // put args values
                    char *arg_buffer;
                    size_t arg_buffer_size;
                    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, mask_m2l, sizeof(unsigned int) * 12,
                            0);
                    task->cl_arg = arg_buffer;
                    task->cl_arg_size = arg_buffer_size;
                    // submit task
                    starpu_task_submit(task);
                }

            } while(octreeIterator.moveRight());
            // move up
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }

        // for all level above leaf-level
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            // M2M
            do{
                //kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                struct starpu_task* const task = starpu_task_create();
                // buffer 0 is parent cell
                task->handles[0] = octreeIterator.getCurrentCell()->handleUp.handle;
                // add child with mask
                unsigned int mask = 0;
                int idxInsert = 1;
                CellClass*const*const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if(child[idxChild]){
                        task->handles[idxInsert++] = child[idxChild]->handleUp.handle;
                        mask = mask | (1 << idxChild);
                    }
                }
                // put right codelet
                task->cl = &m2m_cl[idxInsert-2];
                // put args data
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                        STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                        STARPU_VALUE, &mask, sizeof(mask),
                        0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                // insert task
                starpu_task_submit(task);
            } while(octreeIterator.moveRight());

            // restart this leve
            octreeIterator = avoidGotoLeftIterator;

            // M2L
            do{
                // get interaction list
                const int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                // if not empty
                if(counter){
                    // create task
                    struct starpu_task* const task = starpu_task_create();
                    // buffer 0 is current leaf
                    task->handles[0] = octreeIterator.getCurrentCell()->handleDown.handle;

                    // insert other with a mask
                    memset(mask_m2l, 0, sizeof(unsigned int) * 12);
                    int idxInsert = 1;
                    for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                        if( neighbors[idxNeigh] ){
                            task->handles[idxInsert++] = neighbors[idxNeigh]->handleUp.handle;
                            mask_m2l[ idxNeigh >> 5 ] = mask_m2l[ idxNeigh >> 5 ] | (1 << (idxNeigh & 0x1F));
                        }
                    }
                    // put the right codelet
                    task->cl = &m2l_cl[counter-1];
                    // put args values
                    char *arg_buffer;
                    size_t arg_buffer_size;
                    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, mask_m2l, sizeof(unsigned int) * 12,
                            0);
                    task->cl_arg = arg_buffer;
                    task->cl_arg_size = arg_buffer_size;
                    // submit task
                    starpu_task_submit(task);
                }

            } while(octreeIterator.moveRight());

            // move up
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

        // start from level 2
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for all level until leaf-level
        for(int idxLevel = 2 ; idxLevel < OctreeHeight - 1 ; ++idxLevel ){
            do{
                //kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                struct starpu_task* const task = starpu_task_create();
                // buffer 0 is parent cell
                task->handles[0] = octreeIterator.getCurrentCell()->handleDown.handle;
                // insert children with mask
                unsigned int mask = 0;
                int idxInsert = 1;
                CellClass*const*const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if(child[idxChild]){
                        task->handles[idxInsert++] = child[idxChild]->handleDown.handle;
                        mask = mask | (1 << idxChild);
                    }
                }
                // put right codelet
                task->cl = &l2l_cl[idxInsert-2];
                // add buffer data
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                        STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                        STARPU_VALUE, &mask, sizeof(mask),
                        0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                // submit
                starpu_task_submit(task);

            } while(octreeIterator.moveRight());
            // go down
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

        // for all leaves
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

    // The current solution uses a global array of kernel
    static KernelClass** globalKernels;

    // P2M
    static void p2m_cpu(void *descr[], void *){
        // get leaf
        RealCellClass* const currentCell = (RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);
        // get particles
        DataVector<ParticleClass> particles((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[1]), STARPU_VECTOR_GET_NX(descr[1]));
        // compute
        globalKernels[starpu_worker_get_id()]->P2M( currentCell , &particles );
    }

    // M2M
    static void m2m_cpu(void *descr[], void *cl_arg)
    {
        // init pointers array
        const RealCellClass* child[8];
        memset(child, 0, sizeof(RealCellClass*)*8);
        // get args
        int level;
        unsigned int mask;
        starpu_codelet_unpack_args(cl_arg, &level, &mask);
        int counter = 1;
        // get child
        for(int idxChild = 0 ; idxChild < 8 && mask; ++idxChild){
            if( mask & 1 ){
                child[idxChild] = ((const RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[counter++]));
            }
            mask >>= 1;
        }
        // get parent
        RealCellClass* const currentCell = (RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);
        // compute
        globalKernels[starpu_worker_get_id()]->M2M( currentCell , child , level);
    }

    // M2L
    static void m2l_cpu(void *descr[], void *cl_arg)
    {
        // get args
        int level;
        unsigned int mask[12];
        starpu_codelet_unpack_args(cl_arg, &level, mask);
        // init pointers array
        const RealCellClass* neighbor[343];
        memset(neighbor, 0, 343 * sizeof(RealCellClass*));
        // get interaction list from mask
        int counter = 0;
        for(int idxNeig = 0 ; idxNeig < 343 ; ++idxNeig){
            if(mask[idxNeig >> 5] & (1 << (idxNeig & 0x1F))){
                ++counter;
                neighbor[idxNeig] = ((const RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[counter]));
            }
        }
        // get current cell
        RealCellClass* const currentCell = (RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);
        // compute
        globalKernels[starpu_worker_get_id()]->M2L( currentCell , neighbor, counter, level );
    }

    // L2L
    static void l2l_cpu(void *descr[], void * cl_arg)
    {
        // init pointers data
        RealCellClass* child[8];
        memset(child, 0, sizeof(RealCellClass*)*8);
        // get args
        int level;
        unsigned int mask;
        starpu_codelet_unpack_args(cl_arg, &level, &mask);
        int counter = 1;
        // retrieve child
        for(int idxChild = 0 ; idxChild < 8 && mask; ++idxChild){
            if( mask & 1){
                child[idxChild] = ((RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[counter++]));
            }
            mask >>= 1;
        }
        // get parents
        const RealCellClass* const currentCell = (RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);
        // compute
        globalKernels[starpu_worker_get_id()]->L2L( currentCell , child , level);
    }

    // L2L
    static void l2p_cpu(void *descr[], void *)
    {
        // get leaf
        const RealCellClass* const currentCell = (const RealCellClass*)STARPU_VARIABLE_GET_PTR(descr[0]);
        // get particles
        DataVector<ParticleClass> particles((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[1]), STARPU_VECTOR_GET_NX(descr[1]));
        // compute
        globalKernels[starpu_worker_get_id()]->L2P( currentCell , &particles );
    }

    // P2P
    static void p2p_cpu(void *descr[], void* cl_arg) {
        // retrieve args
        unsigned int mask;
        FTreeCoordinate coordinate;
        starpu_codelet_unpack_args(cl_arg, &mask, &coordinate);
        // pointers data
        DataVector<ParticleClass>* neighors[27];
        memset(neighors, 0, 27 * sizeof(DataVector<ParticleClass>*) );

        // get neigbors by allocating a vector wrapper
        int counter = 0;
        for(int idxNeig = 0 ; idxNeig < 27 && mask; ++idxNeig){
            if(mask & 1){
                ++counter;
                neighors[idxNeig] = new DataVector<ParticleClass>((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[counter]), STARPU_VECTOR_GET_NX(descr[counter]));
            }
            mask >>= 1;
        }
        // get current particles
        DataVector<ParticleClass> currentContainer((ParticleClass*)STARPU_VECTOR_GET_PTR(descr[0]), STARPU_VECTOR_GET_NX(descr[0]));
        // compute
        globalKernels[starpu_worker_get_id()]->P2P(coordinate, &currentContainer, &currentContainer , neighors, counter);
        // dealloc wrapper
        for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
            if(neighors[idxNeig]) delete neighors[idxNeig];
        }
    }

};

template<class OctreeClass, class ParticleClass, class CellClass, class RealCellClass, class ContainerClass, class KernelClass, class LeafClass>
KernelClass** FFmmAlgorithmStarpu<OctreeClass,ParticleClass,CellClass,RealCellClass,ContainerClass,KernelClass,LeafClass>::globalKernels = 0;

#endif //FFMMALGORITHMSTARPU_HPP

