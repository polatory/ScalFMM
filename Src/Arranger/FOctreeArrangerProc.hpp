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
#ifndef FOCTREEARRANGERPROC_HPP
#define FOCTREEARRANGERPROC_HPP

#include "../Utils/FGlobal.hpp"
#include "../Containers/FVector.hpp"
#include "../Utils/FAssertable.hpp"
#include "../Utils/FMpi.hpp"

/** This class is an arranger, it move the particles that need
  * to be hosted in a different leaf
  * This is the parallel version.
  */
template <class OctreeClass, class ContainerClass, class ParticleClass>
class FOctreeArrangerProc : FAssertable {
    /** Interval is the min/max morton index
      * for a proc
      */
    struct Interval{
        MortonIndex min;
        MortonIndex max;
    };


    /** Find the interval that contains mindex */
    int getInterval(const MortonIndex mindex, const int size, const Interval intervals[]) const{
        for(int idxProc = 0 ; idxProc < size ; ++idxProc){
            // does it contains the index?
            if( intervals[idxProc].min <= mindex && mindex < intervals[idxProc].max){
                return idxProc;
            }
        }
        // if no interval found return the lastest one
        return size - 1;
    }

    OctreeClass* const tree;


public:
    /** Basic constructor */
    FOctreeArrangerProc(OctreeClass* const inTree) : tree(inTree) {
        fassert(tree, "Tree cannot be null", __LINE__ , __FILE__ );
    }

    /** return false if the tree is empty after processing */
    bool rearrange(const FMpi::FComm& comm, const bool isPeriodic = false){
        // interval of each procs
        Interval*const intervals = new Interval[comm.processCount()];
        memset(intervals, 0, sizeof(Interval) * comm.processCount());

        {   // We need to exchange interval of each process, this interval
            // will be based on the current morton min max
            Interval myLastInterval;

            // take fist index
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            myLastInterval.min = octreeIterator.getCurrentGlobalIndex();
            // take last index
            octreeIterator.gotoRight();
            myLastInterval.max = octreeIterator.getCurrentGlobalIndex();

            // We get the min/max indexes from each procs
            FMpi::MpiAssert( MPI_Allgather( &myLastInterval, sizeof(Interval), MPI_BYTE, intervals, sizeof(Interval), MPI_BYTE, comm.getComm()),  __LINE__ );

            // increase interval in the empty morton index
            intervals[0].min = 0;
            for(int idxProc = 1 ; idxProc < comm.processCount() ; ++idxProc){
                intervals[idxProc].min = ((intervals[idxProc].min - intervals[idxProc-1].max)/2) + intervals[idxProc-1].max;
                intervals[idxProc-1].max = intervals[idxProc].min;
            }

            intervals[comm.processCount() - 1].max = ((1 << (3*(tree->getHeight()-1))) - 1);
        }

        // Particles that move
        FVector<ParticleClass>*const toMove = new FVector<ParticleClass>[comm.processCount()];

        { // iterate on the leafs and found particle to remove or to send
            // For periodic
            const FReal boxWidth = tree->getBoxWidth();
            const F3DPosition min(tree->getBoxCenter(),-boxWidth/2);
            const F3DPosition max(tree->getBoxCenter(),boxWidth/2);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();

                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                while( iter.hasNotFinished() ){
                    if(isPeriodic){
                        F3DPosition partPos = iter.data().getPosition();
                        while(partPos.getX() < min.getX()){
                            partPos.incX(boxWidth);
                        }
                        while(partPos.getX() >= max.getX()){
                            partPos.incX(-boxWidth);
                        }
                        while(partPos.getY() < min.getY()){
                            partPos.incY(boxWidth);
                        }
                        while(partPos.getY() >= max.getY()){
                            partPos.incY(-boxWidth);
                        }
                        while(partPos.getZ() < min.getZ()){
                            partPos.incZ(boxWidth);
                        }
                        while(partPos.getZ() >= max.getZ()){
                            partPos.incZ(-boxWidth);
                        }
                        iter.data().setPosition(partPos);
                    }
                    const MortonIndex particuleIndex = tree->getMortonFromPosition(iter.data().getPosition());
                    // is this particle need to be changed from its leaf
                    if(particuleIndex != currentIndex){
                        // find the right interval
                        const int procConcerned = getInterval( particuleIndex, comm.processCount(), intervals);
                        toMove[procConcerned].push(iter.data());
                        iter.remove();
                    }
                    else {
                        iter.gotoNext();
                    }
                }
            } while(octreeIterator.moveRight());
        }

        // To send and recv
        ParticleClass* toReceive = 0;
        MPI_Request*const requests = new MPI_Request[comm.processCount()*2];
        memset(requests, 0, sizeof(MPI_Request) * comm.processCount() * 2);
        long long int*const indexToReceive = new long long int[comm.processCount() + 1];
        memset(indexToReceive, 0, sizeof(long long int) * comm.processCount() + 1);

        int iterRequests = 0;
        int limitRecvSend = 0;
        int hasToRecvFrom = 0;

        { // gather what to send to who + isend data
            int*const counter = new int[comm.processCount()];
            memset(counter, 0, sizeof(int) * comm.processCount());

            for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
                counter[idxProc] = toMove[idxProc].getSize();
            }
            // say who send to who
            int*const allcounter = new int[comm.processCount()*comm.processCount()];
            FMpi::MpiAssert( MPI_Allgather( counter, comm.processCount(), MPI_INT, allcounter, comm.processCount(), MPI_INT, comm.getComm()),  __LINE__ );

            // prepare buffer to receive
            long long int sumToRecv = 0;
            indexToReceive[0] = 0;
            for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
                if( idxProc != comm.processId()){
                    sumToRecv += allcounter[idxProc * comm.processCount() + comm.processId()];
                }
                indexToReceive[idxProc + 1] = sumToRecv;
            }
            toReceive = new ParticleClass[sumToRecv];

            // send
            for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
                if(idxProc != comm.processId() && allcounter[idxProc * comm.processCount() + comm.processId()]){
                    FMpi::MpiAssert( MPI_Irecv(&toReceive[indexToReceive[idxProc]], allcounter[idxProc * comm.processCount() + comm.processId()] * int(sizeof(ParticleClass)), MPI_BYTE,
                              idxProc, 0, comm.getComm(), &requests[iterRequests++]),  __LINE__ );
                    hasToRecvFrom += 1;
                }
            }

            limitRecvSend = iterRequests;

            // recv
            for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
                if(idxProc != comm.processId() && toMove[idxProc].getSize()){
                    FMpi::MpiAssert( MPI_Isend(toMove[idxProc].data(), toMove[idxProc].getSize() * int(sizeof(ParticleClass)), MPI_BYTE,
                              idxProc, 0, comm.getComm(), &requests[iterRequests++]),  __LINE__ );
                }
            }

            delete[] allcounter;
            delete[] counter;
        }

        { // insert particles that moved
            for(int idxPart = 0 ; idxPart < toMove[comm.processId()].getSize() ; ++idxPart){
                tree->insert(toMove[comm.processId()][idxPart]);
            }
        }

        {   // wait any recv or send
            // if it is a recv then insert particles
            MPI_Status status;
            while( hasToRecvFrom ){
                int done = 0;
                FMpi::MpiAssert( MPI_Waitany( iterRequests, requests, &done, &status ),  __LINE__ );
                if( done < limitRecvSend ){
                    const int source = status.MPI_SOURCE;
                    for(long long int idxPart = indexToReceive[source] ; idxPart < indexToReceive[source+1] ; ++idxPart){
                        tree->insert(toReceive[idxPart]);
                    }
                    hasToRecvFrom -= 1;
                }
            }            
        }

        int counterLeavesAlive = 0;
        { // Remove empty leaves
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            bool workOnNext = true;

            do{
                // Empty leaf
                if( octreeIterator.getCurrentListTargets()->getSize() == 0 ){
                    const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();
                    workOnNext = octreeIterator.moveRight();
                    tree->removeLeaf( currentIndex );
                }
                // Not empty, just continue
                else {
                    // todo remove
                    if( octreeIterator.getCurrentGlobalIndex() < intervals[comm.processId()].min
                            || intervals[comm.processId()].max <= octreeIterator.getCurrentGlobalIndex()){
                        std::cout << comm.processId() << " Must delete this leaf at " << octreeIterator.getCurrentGlobalIndex()
                                  <<  " size " << octreeIterator.getCurrentListTargets()->getSize()  <<std::endl;
                    }

                    workOnNext = octreeIterator.moveRight();
                    counterLeavesAlive += 1;                    
                }
            } while( workOnNext );
        }

        // wait all send
        FMpi::MpiAssert( MPI_Waitall( iterRequests, requests, MPI_STATUSES_IGNORE),  __LINE__ );

        delete[] intervals;
        delete[] toMove;
        delete[] requests;
        delete[] toReceive;
        delete[] indexToReceive;

        // return false if tree is empty
        return counterLeavesAlive != 0;
    }

};

#endif // FOCTREEARRANGERPROC_HPP
