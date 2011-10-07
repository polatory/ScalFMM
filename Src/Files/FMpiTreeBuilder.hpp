#ifndef FMPITREEBUILDER_H
#define FMPITREEBUILDER_H

#include "../Utils/FMpi.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Utils/FMemUtils.hpp"
#include "../Utils/FTrace.hpp"

/** This class manage the loading of particles
  * for the mpi version.
  * It use a BinLoader and then sort the data
  * with a parallel quick sort.
  */
template<class ParticleClass>
class FMpiTreeBuilder{
    /** This method has been tacken from the octree
      * it computes a tree coordinate (x or y or z) from real position
      */
    static long getTreeCoordinate(const FReal inRelativePosition, const FReal boxWidthAtLeafLevel) {
            const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
            const long index = long(FMath::dfloor(indexFReal));
            if( index && FMath::LookEqual(inRelativePosition, boxWidthAtLeafLevel * FReal(index) ) ){
                    return index - 1;
            }
            return index;
    }

    /** get current rank */
    static int MpiGetRank(MPI_Comm comm = MPI_COMM_WORLD){
        int rank(0);
        MPI_Comm_rank(comm, &rank);
        return rank;
    }

    /** get current nb procs */
    static int MpiGetNbProcs(MPI_Comm comm = MPI_COMM_WORLD){
        int nb(0);
        MPI_Comm_size(comm, &nb);
        return nb;
    }

    /** receive data from a tag function */
    static void receiveDataFromTag(const int inSize, const int inTag, void* const inData, int* const inSource = 0, int* const inFilledSize = 0){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, inTag, MPI_COMM_WORLD, &status);
        if(inSource) *inSource = status.MPI_SOURCE;
        if(inFilledSize) MPI_Get_count(&status,MPI_CHAR,inFilledSize);
    }

    template< class T >
    static T GetLeft(const T inSize) {
        const double step = (double(inSize) / double(MpiGetNbProcs()));
        return T(FMath::Ceil(step * double(MpiGetRank())));
    }

    template< class T >
    static T GetRight(const T inSize) {
        const double step = (double(inSize) / double(MpiGetNbProcs()));
        const T res = T(FMath::Ceil(step * double(MpiGetRank()+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    template< class T >
    static T GetOtherRight(const T inSize, const int other) {
        const double step = (double(inSize) / MpiGetNbProcs());
        const T res = T(FMath::Ceil(step * (other+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    template< class T >
    static T GetOtherLeft(const T inSize, const int other) {
        const double step = (double(inSize) / double(MpiGetNbProcs()));
        return T(FMath::Ceil(step * double(other)));
    }

    template <class T1, class T2>
    static T1 Min( const T1 v1, const T2 v2){
        return (v1 < v2 ? v1 : v2);
    }

    /** This struct is used to represent a particles group to
      * sort them easily
      */
    struct ParticlesGroup {
        int number;
        int positionInArray;
        MortonIndex index;
        ParticlesGroup(const int inNumber = 0 , const int inPositionInArray = 0, const MortonIndex inIndex = 0)
            : number(inNumber), positionInArray(inPositionInArray), index(inIndex) {
        }
    };


    /** A particle may not have a MortonIndex Method
      * but they are sorted on their morton index
      * so this struct store a particle and its index
      */
    struct IndexedParticle{
        MortonIndex index;
        ParticleClass particle;

        operator MortonIndex(){
            return this->index;
        }
    };

    char* intervals;
    FSize nbLeavesInIntervals;

private:
    // Forbid copy
    FMpiTreeBuilder(const FMpiTreeBuilder&){}
    FMpiTreeBuilder& operator=(const FMpiTreeBuilder&){return *this;}

public:
    /** Constructor */
    FMpiTreeBuilder()
        :  intervals(0), nbLeavesInIntervals(0) {
    }

    /** Destructor */
    virtual ~FMpiTreeBuilder(){
        delete [] intervals;
    }

    /** Split and sort the file */
    template <class LoaderClass>
    bool splitAndSortFile(LoaderClass& loader, const int NbLevels){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader" , __FILE__ , __LINE__) );
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();

        // First we create the particles that belong to us (our proc)
        // we compute their morton index to be able to sort them
        //

        IndexedParticle* outputArray = 0;
        FSize outputSize = 0;
        {
            FTRACE( FTrace::FRegion regionTrace("Insert Particles", __FUNCTION__ , __FILE__ , __LINE__) );
            // create particles
            IndexedParticle*const realParticlesIndexed = new IndexedParticle[loader.getNumberOfParticles()];
            FMemUtils::memset(realParticlesIndexed, 0, sizeof(IndexedParticle) * loader.getNumberOfParticles());
            F3DPosition boxCorner(loader.getCenterOfBox() - (loader.getBoxWidth()/2));
            FTreeCoordinate host;

            const FReal boxWidthAtLeafLevel = loader.getBoxWidth() / FReal(1 << (NbLevels - 1) );
            for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(realParticlesIndexed[idxPart].particle);
                host.setX( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
                host.setY( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
                host.setZ( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));

                realParticlesIndexed[idxPart].index = host.getMortonIndex(NbLevels - 1);
            }

            // sort particles
            FQuickSort::QsMpi<IndexedParticle,MortonIndex>(realParticlesIndexed, loader.getNumberOfParticles(),outputArray,outputSize);
            delete [] (realParticlesIndexed);
        }
        // be sure there is no splited leaves
        // to do that we exchange the first index with the left proc
        {
            FTRACE( FTrace::FRegion regionTrace("Remove Splited leaves", __FUNCTION__ , __FILE__ , __LINE__) );

            MortonIndex otherFirstIndex = -1;
            {
                FMpi::Request req[2];
                int reqiter = 0;
                // can I send my first index? == I am not rank 0 & I have data
                if( 0 < rank && outputSize){
                    MPI_Isend( &outputArray[0].index, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs, MPI_COMM_WORLD, &req[reqiter++]);
                }
                if( rank != nbProcs - 1){
                    MPI_Irecv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs, MPI_COMM_WORLD, &req[reqiter++]);
                }

                MPI_Waitall(reqiter,req,MPI_STATUSES_IGNORE);

                // I could not send because I do not have data, so I transmit the data coming
                // from my right neigbors
                if( 0 < rank && !outputSize){
                    MPI_Send( &otherFirstIndex, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs, MPI_COMM_WORLD);
                }
            }

            MPI_Request req[2];
            int reqiter = 0;

            // at this point every one know the first index of his right neighbors
            const bool needToRecvBeforeSend = (rank != 0 && ((outputSize && otherFirstIndex == outputArray[0].index ) || !outputSize));
            if( needToRecvBeforeSend || (rank == nbProcs - 1) ){
                // Here we want to send data we do not have
                // so we first receive other data and put then into the array
                // this case happens only if I have one leaf with index MX
                // and if proc[rank - 1].last.index == MX && proc[rank + 1].first.index == MX

                int sendByOther = 0;

                MPI_Status probStatus;
                MPI_Probe(rank - 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &probStatus);
                MPI_Get_count( &probStatus,  MPI_BYTE, &sendByOther);

                if(sendByOther){
                    sendByOther /= sizeof(IndexedParticle);
                    const IndexedParticle* const reallocOutputArray = outputArray;
                    const FSize reallocOutputSize = outputSize;

                    outputSize += sendByOther;
                    outputArray = new IndexedParticle[outputSize];
                    FMemUtils::memcpy(&outputArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                    delete[] reallocOutputArray;

                    MPI_Recv(outputArray, sizeof(IndexedParticle) * sendByOther, MPI_BYTE, rank - 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &probStatus);
                }
                else{
                    MPI_Irecv(0, 0, MPI_BYTE, rank - 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &req[reqiter++]);
                }
            }

            if(rank != nbProcs - 1){

                FSize idxPart = outputSize - 1 ;
                while(idxPart >= 0 && outputArray[idxPart].index == otherFirstIndex){
                    --idxPart;
                }
                const int toSend = int(outputSize - 1 - idxPart);
                MPI_Isend( &outputArray[idxPart + 1], toSend * sizeof(IndexedParticle), MPI_BYTE, rank + 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &req[reqiter++]);

                if( rank != 0 && !needToRecvBeforeSend && (rank != nbProcs - 1)){
                    int sendByOther = 0;

                    MPI_Status probStatus;
                    MPI_Probe(rank - 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sendByOther);

                    if(sendByOther){
                        sendByOther /= sizeof(IndexedParticle);
                        char* const tempBuffer = new char[sizeof(IndexedParticle) * sendByOther];

                        MPI_Irecv(tempBuffer, sizeof(IndexedParticle) * sendByOther, MPI_BYTE, rank - 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &req[reqiter++]);

                        MPI_Waitall(reqiter,req, MPI_STATUSES_IGNORE);
                        reqiter = 0;

                        const IndexedParticle* const reallocOutputArray = outputArray;
                        const FSize reallocOutputSize = outputSize;

                        outputSize += sendByOther;
                        outputArray = new IndexedParticle[outputSize];
                        FMemUtils::memcpy(&outputArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                        delete[] reallocOutputArray;
                        memcpy(outputArray, tempBuffer, sendByOther * sizeof(IndexedParticle));
                        delete[] tempBuffer;
                    }
                    else{
                        MPI_Irecv( 0, 0, MPI_BYTE, rank - 1, FMpi::TagSplittedLeaf, MPI_COMM_WORLD, &req[reqiter++]);
                    }
                }
            }

            MPI_Waitall(reqiter,req,MPI_STATUSES_IGNORE);
        }

        // We now copy the data from a sorted type into real particles array + counter

        nbLeavesInIntervals = 0;
        if(outputSize){
            intervals = new char[outputSize * (sizeof(ParticleClass) + sizeof(int))];

            MortonIndex previousIndex = -1;
            char* writeIndex = intervals;
            int* writeCounter = 0;

            for( FSize idxPart = 0; idxPart < outputSize ; ++idxPart){
                if( outputArray[idxPart].index != previousIndex ){
                    previousIndex = outputArray[idxPart].index;
                    ++nbLeavesInIntervals;

                    writeCounter = reinterpret_cast<int*>( writeIndex );
                    writeIndex += sizeof(int);

                    (*writeCounter) = 0;
                }

                memcpy(writeIndex, &outputArray[idxPart].particle, sizeof(ParticleClass));

                writeIndex += sizeof(ParticleClass);
                ++(*writeCounter);
            }
        }
        delete [] outputArray;

        return true;
    }

    /** Put the interval into a tree */
    template <class OctreeClass>
    bool intervalsToTree(OctreeClass& realTree){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader" , __FILE__ , __LINE__) );
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();
        const FSize currentNbLeafs = nbLeavesInIntervals;

        int leavesPerProcs[nbProcs];
        memset(leavesPerProcs, 0, sizeof(int) * nbProcs);
        MPI_Allgather(&nbLeavesInIntervals, 1, MPI_INT, leavesPerProcs, 1, MPI_INT, MPI_COMM_WORLD);

        FSize currentLeafsOnMyLeft  = 0;
        FSize currentLeafsOnMyRight = 0;
        for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc ){
            if(idxProc < rank) currentLeafsOnMyLeft  += leavesPerProcs[idxProc];
            if(rank < idxProc) currentLeafsOnMyRight += leavesPerProcs[idxProc];
        }

        const FSize totalNbLeafs     = (currentLeafsOnMyLeft + currentNbLeafs + currentLeafsOnMyRight);
        const FSize nbLeafsOnMyLeft  = GetLeft(totalNbLeafs);
        const FSize nbLeafsOnMyRight = GetRight(totalNbLeafs);


        int leavesToSend[nbProcs];
        memset(leavesPerProcs, 0, sizeof(int) * nbProcs);
        int bytesToSend[nbProcs];
        memset(bytesToSend, 0, sizeof(int) * nbProcs);
        int bytesOffset[nbProcs];
        memset(bytesToSend, 0, sizeof(int) * nbProcs);

        FSize currentIntervalPosition = 0;

        const bool iNeedToSendToLeft  = currentLeafsOnMyLeft < nbLeafsOnMyLeft;
        const bool iNeedToSendToRight = nbLeafsOnMyRight < currentLeafsOnMyLeft + currentNbLeafs;

        int leftProcToStartSend = rank;
        if(iNeedToSendToLeft){
            int idxProc = rank - 1;
            while( idxProc >= 0 ){
                const FSize thisProcRight = GetOtherRight(totalNbLeafs, idxProc);
                if( currentLeafsOnMyLeft < thisProcRight ){
                    break;
                }
                --idxProc;
            }

            leftProcToStartSend = idxProc;
            int hasToGive = currentNbLeafs;
            leavesToSend[idxProc] = Min(GetOtherRight(totalNbLeafs, idxProc) - GetOtherLeft(totalNbLeafs, idxProc) , hasToGive);
            {
                bytesOffset[idxProc] = 0;
                for(FSize idxLeaf = 0 ; idxLeaf < leavesToSend[idxProc] ; ++idxLeaf){
                    currentIntervalPosition += ((*(int*)&intervals[currentIntervalPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                bytesToSend[idxProc] = currentIntervalPosition - bytesOffset[idxProc];
            }
            hasToGive -= leavesToSend[idxProc];
            ++idxProc;

            while(idxProc < rank && hasToGive){
                leavesToSend[idxProc] = Min( GetOtherRight(totalNbLeafs, idxProc) - GetOtherLeft(totalNbLeafs, idxProc), hasToGive);

                bytesOffset[idxProc] = currentIntervalPosition;
                for(FSize idxLeaf = 0 ; idxLeaf < leavesToSend[idxProc] ; ++idxLeaf){
                    currentIntervalPosition += ((*(int*)&intervals[currentIntervalPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                bytesToSend[idxProc] = currentIntervalPosition - bytesOffset[idxProc];

                hasToGive -= leavesToSend[idxProc];
                ++idxProc;
            }
        }

        {
            const FSize iNeedToSendLeftCount = nbLeafsOnMyLeft - currentLeafsOnMyLeft;
            const FSize iNeedToSendRightCount = currentLeafsOnMyLeft + currentNbLeafs - nbLeafsOnMyRight;
            FSize endForMe = currentNbLeafs;
            if(iNeedToSendToRight) endForMe -= iNeedToSendRightCount;

            for(FSize idxLeaf = iNeedToSendLeftCount ; idxLeaf < endForMe ; ++idxLeaf){
                const int nbPartInLeaf = (*(int*)&intervals[currentIntervalPosition]);
                ParticleClass* const particles = reinterpret_cast<ParticleClass*>(&intervals[currentIntervalPosition] + sizeof(int));

                for(int idxPart = 0 ; idxPart < nbPartInLeaf ; ++idxPart){
                    realTree.insert(particles[idxPart]);
                }
                currentIntervalPosition += (nbPartInLeaf * sizeof(ParticleClass)) + sizeof(int);
            }
        }

        int rightProcToStartSend = rank;
        if(iNeedToSendToRight){
            int idxProc = rank + 1;
            while( idxProc + 1 < nbProcs ){
                const FSize thisProcLeft = GetOtherLeft(totalNbLeafs, idxProc + 1);
                if( totalNbLeafs - currentLeafsOnMyRight < thisProcLeft ){
                    break;
                }
                ++idxProc;
            }

            rightProcToStartSend = idxProc;
            int hasToGive = currentNbLeafs;
            leavesToSend[idxProc] = Min((totalNbLeafs - currentLeafsOnMyRight) - GetOtherLeft(totalNbLeafs, idxProc) , hasToGive);

            {
                bytesOffset[idxProc] = currentIntervalPosition;
                for(FSize idxLeaf = 0 ; idxLeaf < leavesToSend[idxProc] ; ++idxLeaf){
                    currentIntervalPosition += ((*(int*)&intervals[currentIntervalPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                bytesToSend[idxProc] = currentIntervalPosition - bytesOffset[idxProc];
            }
            hasToGive -= leavesToSend[idxProc];
            --idxProc;

            while(rank < idxProc && hasToGive){
                leavesToSend[idxProc] = Min( GetOtherRight(totalNbLeafs, idxProc) - GetOtherLeft(totalNbLeafs, idxProc), hasToGive);

                bytesOffset[idxProc] = currentIntervalPosition;
                for(FSize idxLeaf = 0 ; idxLeaf < leavesToSend[idxProc] ; ++idxLeaf){
                    currentIntervalPosition += ((*(int*)&intervals[currentIntervalPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                bytesToSend[idxProc] = currentIntervalPosition - bytesOffset[idxProc];

                hasToGive -= leavesToSend[idxProc];
                --idxProc;
            }
        }

        int bytesToSendRecv[nbProcs * nbProcs];
        memset(bytesToSendRecv, 0, sizeof(int) * nbProcs * nbProcs);
        MPI_Allgather(bytesToSend, nbProcs, MPI_INT, bytesToSendRecv, nbProcs, MPI_INT, MPI_COMM_WORLD);

        int bytesToRecv[nbProcs];
        memset(bytesToRecv, 0, sizeof(int) * nbProcs);
        int bytesOffsetToRecv[nbProcs];
        memset(bytesOffsetToRecv, 0, sizeof(int) * nbProcs);

        FSize sumBytesToRecv = 0;
        for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
            if( bytesToSendRecv[idxProc * nbProcs + rank] ){
                bytesOffsetToRecv[idxProc] = sumBytesToRecv;
                sumBytesToRecv += bytesToSendRecv[idxProc * nbProcs + rank];
                bytesToRecv[idxProc] = bytesToSendRecv[idxProc * nbProcs + rank];
            }
        }

        char* const recvbuf = new char[sumBytesToRecv];

        MPI_Alltoallv(intervals, bytesToSend, bytesOffset, MPI_BYTE,
                      recvbuf, bytesToRecv, bytesOffsetToRecv, MPI_BYTE,
                      MPI_COMM_WORLD);

        FSize recvBufferPosition = 0;
        while( recvBufferPosition < sumBytesToRecv){
            const int nbPartInLeaf = (*reinterpret_cast<int*>(&recvbuf[recvBufferPosition]));
            ParticleClass* const particles = reinterpret_cast<ParticleClass*>(&recvbuf[recvBufferPosition] + sizeof(int));

            for(int idxPart = 0 ; idxPart < nbPartInLeaf ; ++idxPart){
                realTree.insert(particles[idxPart]);
            }
            recvBufferPosition += (nbPartInLeaf * sizeof(ParticleClass)) + sizeof(int);
        }

        delete[] recvbuf;

        return true;
    }
};

#endif // FMPITREEBUILDER_H
