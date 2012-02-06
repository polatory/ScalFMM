// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#ifndef FFMMALGORITHMTASK_HPP
#define FFMMALGORITHMTASK_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmTask
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmTask : protected FAssertable{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass** kernels;    //< The kernels

    const int MaxThreads;

    const int OctreeHeight;

public:	
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmTask(OctreeClass* const inTree, KernelClass* const inKernels)
                      : tree(inTree) , kernels(0),
                        MaxThreads(omp_get_max_threads()), OctreeHeight(tree->getHeight()) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(inKernels, "kernels cannot be null", __LINE__, __FILE__);

        this->kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFmmAlgorithmTask\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmTask(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        bottomPass();

        upwardPass();

        downardPass();

        directPass();         
    }

private:
    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        #pragma omp parallel
        {
            KernelClass*const myThreadkernels = kernels[omp_get_thread_num()];

            #pragma omp single nowait
            {
                typename OctreeClass::Iterator octreeIterator(tree);

                // Iterate on leafs
                octreeIterator.gotoBottomLeft();
                do{
                    // We need the current cell that represent the leaf
                    // and the list of particles
                    #pragma omp task
                    {
                        myThreadkernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
                    }

                } while(octreeIterator.moveRight());

                #pragma omp taskwait
            }
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        #pragma omp parallel
        {
            KernelClass*const myThreadkernels = kernels[omp_get_thread_num()];

            #pragma omp single nowait
            {
                // Start from leal level - 1
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();
                octreeIterator.moveUp();

                typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

                // for each levels
                for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
                    // for each cells
                    do{
                        // We need the current cell and the child
                        // child is an array (of 8 child) that may be null
                        #pragma omp task
                        {
                            myThreadkernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                        }
                    } while(octreeIterator.moveRight());

                    avoidGotoLeftIterator.moveUp();
                    octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

                    #pragma omp taskwait
                }
            }
        }


        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        { // first M2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            #pragma omp parallel
            {
                KernelClass*const myThreadkernels = kernels[omp_get_thread_num()];
                const CellClass* neighbors[189];

                #pragma omp single nowait
                {
                    typename OctreeClass::Iterator octreeIterator(tree);
                    octreeIterator.moveDown();

                    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);


                    // for each levels
                    for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                        // for each cells
                        do{
                            const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                            if(counter){
                                #pragma omp task
                                {
                                    myThreadkernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);
                                }
                            }
                        } while(octreeIterator.moveRight());

                        avoidGotoLeftIterator.moveDown();
                        octreeIterator = avoidGotoLeftIterator;

                        #pragma omp taskwait
                    }
                }
            }
            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        }

        { // second L2L            
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);

            #pragma omp parallel
            {
                KernelClass*const myThreadkernels = kernels[omp_get_thread_num()];

                #pragma omp single nowait
                {
                    typename OctreeClass::Iterator octreeIterator(tree);
                    octreeIterator.moveDown();

                    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

                    const int heightMinusOne = OctreeHeight - 1;
                    // for each levels exepted leaf level
                    for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                        // for each cells
                        do{
                            #pragma omp task
                            {
                                myThreadkernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                            }

                        } while(octreeIterator.moveRight());

                        avoidGotoLeftIterator.moveDown();
                        octreeIterator = avoidGotoLeftIterator;

                        #pragma omp taskwait
                    }
                }
            }

            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        }


    }

    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        const int heightMinusOne = OctreeHeight - 1;

        #pragma omp parallel
        {
            KernelClass*const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[26];
            MortonIndex neighborsIndex[26];

            #pragma omp single nowait
            {
                const int SizeShape = 3*3*3;
                FVector<typename OctreeClass::Iterator> shapes[SizeShape];

                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();

                // for each leafs
                do{
                    #pragma omp task
                    {
                        myThreadkernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
                    }

                    const FTreeCoordinate& coord = octreeIterator.getCurrentGlobalCoordinate();
                    const int shapePosition = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);

                    if( shapePosition == 0){
                        #pragma omp task
                        {
                            // need the current particles and neighbors particles
                            const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, octreeIterator.getCurrentGlobalIndex(),heightMinusOne);
                            myThreadkernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                        }
                    }
                    else{
                        shapes[shapePosition].push(octreeIterator);
                    }

                } while(octreeIterator.moveRight());

                #pragma omp taskwait

                for( int idxShape = 1 ; idxShape < SizeShape ; ++idxShape){
                    int iterLeaf = shapes[idxShape].getSize();
                    while( iterLeaf-- ){
                        typename OctreeClass::Iterator toWork = shapes[idxShape][iterLeaf];
                        #pragma omp task
                        {
                            const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, toWork.getCurrentGlobalIndex(),heightMinusOne);
                            myThreadkernels->P2P(toWork.getCurrentGlobalIndex(),toWork.getCurrentListTargets(), toWork.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                        }
                    }

                    #pragma omp taskwait
                }
            }
        }


        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

};


#endif //FFMMALGORITHMTASK_HPP


