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

        //////////////////////////////////////////////////////////////////////////////////
        // We inform the master proc about the data we have
        //////////////////////////////////////////////////////////////////////////////////
        FTRACE( FTrace::FRegion preprocessTrace("Preprocess", __FUNCTION__ , __FILE__ , __LINE__) );

        FSize nbLeafs = nbLeavesInIntervals;
        FSize leftLeafs = 0;
        FSize rightLeafs = 0;

        // receive from left and right
        if((rank == 0)){
            MPI_Request requests[2];
            MPI_Isend(&nbLeafs, sizeof(FSize), MPI_BYTE , 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(&rightLeafs, sizeof(FSize), MPI_BYTE, 1 , FMpi::TagExchangeNbLeafs , MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
        }
        else if(rank == nbProcs - 1){
            MPI_Request requests[2];
            MPI_Isend(&nbLeafs, sizeof(FSize), MPI_BYTE , rank - 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(&leftLeafs, sizeof(FSize), MPI_BYTE, rank - 1 , FMpi::TagExchangeNbLeafs , MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
        }
        else { //rank != 0) && rank != nbProcs - 1
            for(int idxToReceive = 0 ; idxToReceive < 2 ; ++idxToReceive){
                int source(0);
                FSize temp = 0;
                receiveDataFromTag(sizeof(FSize), FMpi::TagExchangeNbLeafs, &temp, &source);
                if(source < rank){ // come from left
                    leftLeafs = temp;
                    temp += nbLeafs;
                    MPI_Send(&temp, sizeof(FSize), MPI_BYTE , rank + 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD);
                }
                else { // come from right
                    rightLeafs = temp;
                    temp += nbLeafs;
                    MPI_Send(&temp, sizeof(FSize), MPI_BYTE , rank - 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD);
                }
            }
        }
        FTRACE( preprocessTrace.end() );
        //////////////////////////////////////////////////////////////////////////////////
        // We balance the data
        //////////////////////////////////////////////////////////////////////////////////

        const FSize totalNbLeafs = (leftLeafs + nbLeafs + rightLeafs);
        const FSize myLeftLeaf = GetLeft(totalNbLeafs);
        const FSize myRightLeaf = GetRight(totalNbLeafs);

        const bool iNeedToSendToLeft = leftLeafs < myLeftLeaf;
        const bool iNeedToSendToRight = myRightLeaf < leftLeafs + nbLeafs;

        const bool iWillReceiveFromRight = leftLeafs + nbLeafs < myRightLeaf;
        const bool iWillReceiveFromLeft = leftLeafs > myLeftLeaf;

        const bool iDoNotHaveEnoughtToSendRight = myRightLeaf < leftLeafs;
        const bool iDoNotHaveEnoughtToSendLeft = leftLeafs + nbLeafs < myLeftLeaf;


        const FSize iNeedToSendLeftCount = myLeftLeaf - leftLeafs;
        const FSize iCanSendToLeft = nbLeafs;

        const FSize iNeedToSendRightCount = leftLeafs + nbLeafs - myRightLeaf;
        const FSize iCanSendToRight = nbLeafs;

        MPI_Request requests[2];
        int iterRequest = 0;

        FSize hasBeenSentToLeft = 0;
        FSize hasBeenSentToRight = 0;

        char* particlesToSend = 0;

        ///////////////////////////////
        // Manage data we already have
        ///////////////////////////////
        FTRACE( FTrace::FRegion step1Trace("Step1", __FUNCTION__ , __FILE__ , __LINE__) );

        if(nbLeafs){
            particlesToSend = intervals;

            FSize currentLeafPosition = 0;

            //Send to Left (the first leaves
            if(iNeedToSendToLeft){
                for(FSize idxLeaf = 0 ; idxLeaf < iNeedToSendLeftCount && idxLeaf < iCanSendToLeft ; ++idxLeaf){
                    currentLeafPosition += ((*(int*)&particlesToSend[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                hasBeenSentToLeft = FMath::Min(iNeedToSendLeftCount, iCanSendToLeft);
                MPI_Isend(particlesToSend, int(currentLeafPosition), MPI_BYTE , rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);
            }

            // Insert the particles I host and that belong to me
            const FSize beginForMe = (iNeedToSendToLeft ? FMath::Min(iNeedToSendLeftCount,iCanSendToLeft) : 0);
            const FSize endForMe = nbLeafs - (iNeedToSendToRight ? FMath::Min(iNeedToSendRightCount,iCanSendToRight) : 0);

            for(FSize idxLeaf = beginForMe ; idxLeaf < endForMe ; ++idxLeaf){

                const int nbPartInLeaf = (*(int*)&particlesToSend[currentLeafPosition]);
                ParticleClass* const particles = reinterpret_cast<ParticleClass*>(&particlesToSend[currentLeafPosition] + sizeof(int));

                for(int idxPart = 0 ; idxPart < nbPartInLeaf ; ++idxPart){
                    realTree.insert(particles[idxPart]);
                }

                currentLeafPosition += (nbPartInLeaf * sizeof(ParticleClass)) + sizeof(int);
            }

            //Send to Right (the right-est leaves
            if(iNeedToSendToRight){
                const FSize beginWriteIndex = currentLeafPosition;

                for(int idxLeaf = 0 ; idxLeaf < iNeedToSendRightCount && idxLeaf < iCanSendToRight ; ++idxLeaf){
                    currentLeafPosition += (*(int*)&particlesToSend[currentLeafPosition]* sizeof(ParticleClass)) + sizeof(int);
                }

                hasBeenSentToRight = FMath::Min(iNeedToSendRightCount, iCanSendToRight);
                MPI_Isend( &particlesToSend[beginWriteIndex], int(currentLeafPosition - beginWriteIndex), MPI_BYTE , rank + 1, FMpi::TagSandSettling,
                          MPI_COMM_WORLD, &requests[iterRequest++]);
            }
        }

        char* toRecvFromLeft = 0;
        char* toRecvFromRight = 0;
        int countReceive = int(iWillReceiveFromLeft) + int(iWillReceiveFromRight);
        int sizeOfLeftBuffer = 0;
        int sizeOfRightBuffer = 0;
        int sizeOfRightData = 0;
        int sizeOfLeftData = 0;

        int sourceToWhileRecv = MPI_ANY_SOURCE;

        // Now prepare to receive data
        while(countReceive--){
            MPI_Status recvStatus;
            MPI_Probe(sourceToWhileRecv, FMpi::TagSandSettling, MPI_COMM_WORLD, &recvStatus);
            // receive from left
            if(recvStatus.MPI_SOURCE == rank - 1){
                MPI_Get_count( &recvStatus,  MPI_BYTE, &sizeOfLeftBuffer);
                toRecvFromLeft = new char[sizeOfLeftBuffer];
                sizeOfLeftData = sizeOfLeftBuffer;
                MPI_Irecv(toRecvFromLeft, sizeOfLeftBuffer, MPI_BYTE, rank - 1 , FMpi::TagSandSettling , MPI_COMM_WORLD, &requests[iterRequest++]);
                sourceToWhileRecv = rank + 1;
            }
            // receive from right
            else{
                MPI_Get_count( &recvStatus,  MPI_BYTE, &sizeOfRightBuffer);
                toRecvFromRight = new char[sizeOfRightBuffer];
                sizeOfRightData = sizeOfRightBuffer;
                MPI_Irecv(toRecvFromRight, sizeOfRightBuffer, MPI_BYTE, rank + 1 , FMpi::TagSandSettling , MPI_COMM_WORLD, &requests[iterRequest++]);
                sourceToWhileRecv = rank - 1;
            }
        }

        ///////////////////////////////
        // Wait send receive
        ///////////////////////////////
        MPI_Waitall(iterRequest, requests, MPI_STATUSES_IGNORE);
        // We can delete the buffer use to send our particles only
        FTRACE( step1Trace.end() );
        FTRACE( FTrace::FRegion step2Trace("Step2", __FUNCTION__ , __FILE__ , __LINE__) );
        ///////////////////////////////
        // Process received data
        // and transfer if needed
        ///////////////////////////////
        // We have to receive from right and transfere to left
        int hasToBeReceivedFromLeft = int(leftLeafs - myLeftLeaf);
        int hasToBeReceivedFromRight = int(myRightLeaf - (leftLeafs + nbLeafs));
        int arrayIdxRight = 0;
        int arrayIdxLeft = 0;

        if(iDoNotHaveEnoughtToSendLeft){

            do{
                arrayIdxRight = 0;
                while(arrayIdxRight < sizeOfRightData && hasBeenSentToLeft < iNeedToSendLeftCount){
                    const int particlesInThisLeaf = *(int*)&toRecvFromRight[arrayIdxRight];
                    arrayIdxRight += sizeof(int) + sizeof(ParticleClass) * particlesInThisLeaf;
                    --hasToBeReceivedFromRight;
                    ++hasBeenSentToLeft;
                }

                MPI_Send(toRecvFromRight, arrayIdxRight, MPI_BYTE , rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD);
                if(hasBeenSentToLeft < iNeedToSendLeftCount){
                    MPI_Status probStatus;
                    MPI_Probe(MPI_ANY_SOURCE, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfRightData);
                    if(sizeOfRightBuffer < sizeOfRightData){
                        sizeOfRightBuffer = sizeOfRightData;
                        delete[] toRecvFromRight;
                        toRecvFromRight = new char[sizeOfRightData];
                    }

                    MPI_Recv(toRecvFromRight, sizeOfRightData, MPI_BYTE, rank + 1 , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasBeenSentToLeft < iNeedToSendLeftCount);
        }
        // We have to receive from left and transfere to right
        else if(iDoNotHaveEnoughtToSendRight){

            do{
                arrayIdxLeft = 0;
                while(arrayIdxLeft < sizeOfLeftData && hasBeenSentToRight < iNeedToSendRightCount){
                    const int particlesInThisLeaf = *(int*)&toRecvFromLeft[arrayIdxLeft];
                    arrayIdxLeft += sizeof(int) + sizeof(ParticleClass) * particlesInThisLeaf;
                    --hasToBeReceivedFromLeft;
                    ++hasBeenSentToRight;
                }

                MPI_Send(toRecvFromLeft, arrayIdxLeft, MPI_BYTE , rank + 1, FMpi::TagSandSettling, MPI_COMM_WORLD);
                if(hasBeenSentToRight < iNeedToSendRightCount){
                    MPI_Status probStatus;
                    MPI_Probe(MPI_ANY_SOURCE, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfLeftData);
                    if(sizeOfLeftBuffer < sizeOfLeftData){
                        sizeOfLeftBuffer = sizeOfLeftData;
                        delete[] toRecvFromLeft;
                        toRecvFromLeft = new char[sizeOfLeftData];
                    }
                    MPI_Recv(toRecvFromLeft, sizeOfLeftData, MPI_BYTE, rank - 1 , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasBeenSentToRight < iNeedToSendRightCount);
        }

        if(iWillReceiveFromLeft ){ // I need to wait
            do{
                while(arrayIdxLeft < sizeOfLeftData){
                    const int particlesInThisLeaf = *(int*)&toRecvFromLeft[arrayIdxLeft];
                    arrayIdxLeft += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&toRecvFromLeft[arrayIdxLeft]);
                    arrayIdxLeft += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }

                    --hasToBeReceivedFromLeft;
                }
                arrayIdxLeft = 0;

                if(hasToBeReceivedFromLeft){
                    MPI_Status probStatus;
                    MPI_Probe( rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfLeftData);
                    if(sizeOfLeftBuffer < sizeOfLeftData){
                        sizeOfLeftBuffer = sizeOfLeftData;
                        delete[] toRecvFromLeft;
                        toRecvFromLeft = new char[sizeOfLeftData];
                    }
                    MPI_Recv(toRecvFromLeft, sizeOfLeftData, MPI_BYTE, rank - 1 , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasToBeReceivedFromLeft);
        }

        if(iWillReceiveFromRight){
            do{
                while(arrayIdxRight < sizeOfRightData){
                    const int particlesInThisLeaf = *(int*)&toRecvFromRight[arrayIdxRight];
                    arrayIdxRight += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&toRecvFromRight[arrayIdxRight]);
                    arrayIdxRight += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }

                    --hasToBeReceivedFromRight;
                }
                arrayIdxRight = 0;

                if(hasToBeReceivedFromRight){
                    MPI_Status probStatus;

                    MPI_Probe( rank + 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfRightData);

                    if(sizeOfRightBuffer < sizeOfRightData){
                        sizeOfRightBuffer = sizeOfRightData;
                        delete[] toRecvFromRight;
                        toRecvFromRight = new char[sizeOfRightData];
                    }
                    MPI_Recv(toRecvFromRight, sizeOfRightData, MPI_BYTE, rank + 1 , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasToBeReceivedFromRight);
        }

        // Delete reception buffers
        delete[] toRecvFromLeft;
        delete[] toRecvFromRight;

        return true;
    }

    /*
    template <class OctreeClass>
    bool intervalsToTree(OctreeClass& realTree){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader" , __FILE__ , __LINE__) );
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();

        //////////////////////////////////////////////////////////////////////////////////
        // We inform the master proc about the data we have
        //////////////////////////////////////////////////////////////////////////////////
        FTRACE( FTrace::FRegion preprocessTrace("Preprocess", __FUNCTION__ , __FILE__ , __LINE__) );

        FSize currentNbLeafs = nbLeavesInIntervals;
        FSize currentLeafsOnMyLeft = 0;
        FSize currentLeafsOnMyRight = 0;

        // receive from left and right
        if((rank == 0)){
            MPI_Request requests[2];
            MPI_Isend(&currentNbLeafs, sizeof(FSize), MPI_BYTE , 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(&currentLeafsOnMyRight, sizeof(FSize), MPI_BYTE, 1 , FMpi::TagExchangeNbLeafs , MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
        }
        else if(rank == nbProcs - 1){
            MPI_Request requests[2];
            MPI_Isend(&currentNbLeafs, sizeof(FSize), MPI_BYTE , rank - 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(&currentLeafsOnMyLeft, sizeof(FSize), MPI_BYTE, rank - 1 , FMpi::TagExchangeNbLeafs , MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
        }
        else { //rank != 0) && rank != nbProcs - 1
            for(int idxToReceive = 0 ; idxToReceive < 2 ; ++idxToReceive){
                int source(0);
                FSize temp = 0;
                receiveDataFromTag(sizeof(FSize), FMpi::TagExchangeNbLeafs, &temp, &source);
                if(source < rank){ // come from left
                    currentLeafsOnMyLeft = temp;
                    temp += currentNbLeafs;
                    MPI_Send(&temp, sizeof(FSize), MPI_BYTE , rank + 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD);
                }
                else { // come from right
                    currentLeafsOnMyRight = temp;
                    temp += currentNbLeafs;
                    MPI_Send(&temp, sizeof(FSize), MPI_BYTE , rank - 1, FMpi::TagExchangeNbLeafs, MPI_COMM_WORLD);
                }
            }
        }
        FTRACE( preprocessTrace.end() );
        //////////////////////////////////////////////////////////////////////////////////
        // We balance the data
        //////////////////////////////////////////////////////////////////////////////////

        const FSize totalNbLeafs = (currentLeafsOnMyLeft + currentNbLeafs + currentLeafsOnMyRight);
        const FSize nbLeafsOnMyLeft = GetLeft(totalNbLeafs);
        const FSize nbLeafsOnMyRight = GetRight(totalNbLeafs);

        const bool iNeedToSendToLeft = currentLeafsOnMyLeft < nbLeafsOnMyLeft;
        const bool iNeedToSendToRight = nbLeafsOnMyRight < currentLeafsOnMyLeft + currentNbLeafs;

        const bool iWillReceiveFromRight = currentLeafsOnMyLeft + currentNbLeafs < nbLeafsOnMyRight;
        const bool iWillReceiveFromLeft = currentLeafsOnMyLeft > nbLeafsOnMyLeft;

        const bool iDoNotHaveEnoughtToSendRight = nbLeafsOnMyRight < currentLeafsOnMyLeft;
        const bool iDoNotHaveEnoughtToSendLeft = currentLeafsOnMyLeft + currentNbLeafs < nbLeafsOnMyLeft;


        const FSize iNeedToSendLeftCount = nbLeafsOnMyLeft - currentLeafsOnMyLeft;
        const FSize iCanSendToLeft = currentNbLeafs;

        const FSize iNeedToSendRightCount = currentLeafsOnMyLeft + currentNbLeafs - nbLeafsOnMyRight;
        const FSize iCanSendToRight = currentNbLeafs;

        int hasToBeReceivedFromLeft  = int(currentLeafsOnMyLeft - nbLeafsOnMyLeft);
        int hasToBeReceivedFromRight = int(nbLeafsOnMyRight - (currentLeafsOnMyLeft + currentNbLeafs));

        ///////////////////////////////
        // Manage data we already have
        ///////////////////////////////
        FTRACE( FTrace::FRegion step1Trace("Step1", __FUNCTION__ , __FILE__ , __LINE__) );
        char* particlesToSend = intervals;

        // We have enought to send right and left it means we will not receive anything
        if( !iWillReceiveFromRight && !iWillReceiveFromLeft ){
            printf("I may send %d to left and %d to right\n", iNeedToSendLeftCount, iNeedToSendRightCount);
            MPI_Request requests[2];
            int iterRequest = 0;

            FSize currentLeafPosition = 0;

            //Send to Left (the first leaves)
            if(iNeedToSendToLeft){
                for(FSize idxLeaf = 0 ; idxLeaf < iNeedToSendLeftCount && idxLeaf < iCanSendToLeft ; ++idxLeaf){
                    currentLeafPosition += ((*(int*)&particlesToSend[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                MPI_Isend(particlesToSend, int(currentLeafPosition), MPI_BYTE , rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);
            }

            // Insert the particles I host and that belong to me
            FSize endForMe = currentNbLeafs;
            if(iNeedToSendToRight) endForMe -= iNeedToSendRightCount;

            for(FSize idxLeaf = iNeedToSendLeftCount ; idxLeaf < endForMe ; ++idxLeaf){
                const int nbPartInLeaf = (*(int*)&particlesToSend[currentLeafPosition]);
                ParticleClass* const particles = reinterpret_cast<ParticleClass*>(&particlesToSend[currentLeafPosition] + sizeof(int));

                for(int idxPart = 0 ; idxPart < nbPartInLeaf ; ++idxPart){
                    realTree.insert(particles[idxPart]);
                }
                currentLeafPosition += (nbPartInLeaf * sizeof(ParticleClass)) + sizeof(int);
            }

            //Send to Right (the right-est leaves)
            if(iNeedToSendToRight){
                const FSize beginWriteIndex = currentLeafPosition;
                for(int idxLeaf = 0 ; idxLeaf < iNeedToSendRightCount && idxLeaf < iCanSendToRight ; ++idxLeaf){
                    currentLeafPosition += (*(int*)&particlesToSend[currentLeafPosition]* sizeof(ParticleClass)) + sizeof(int);
                }
                MPI_Isend( &particlesToSend[beginWriteIndex], int(currentLeafPosition - beginWriteIndex), MPI_BYTE , rank + 1, FMpi::TagSandSettling,
                          MPI_COMM_WORLD, &requests[iterRequest++]);
            }
            printf("Wait all\n");
            MPI_Waitall( iterRequest, requests, MPI_STATUSES_IGNORE);
        }
        // I will receive from left and right means I will not send anything
        else if(iWillReceiveFromLeft && iWillReceiveFromRight){
            // I receive from both part and insert
            printf("I may receive %d from left and %d from right\n", hasToBeReceivedFromLeft, hasToBeReceivedFromRight);

            char* buffer    = 0;
            int bufferSize  = 0;

            while(hasToBeReceivedFromLeft || hasToBeReceivedFromRight){
                MPI_Status probStatus;
                MPI_Probe( MPI_ANY_SOURCE, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);

                int sizeOfLeftData    = 0;
                MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfLeftData);
                if(bufferSize < sizeOfLeftData){
                    bufferSize = sizeOfLeftData;
                    delete[] buffer;
                    buffer = new char[bufferSize];
                }
                MPI_Recv(buffer, sizeOfLeftData, MPI_BYTE, probStatus.MPI_SOURCE , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int nbLeafsInThisBuffer = 0;
                int currentIdxArray = 0;
                while(currentIdxArray < sizeOfLeftData){
                    const int particlesInThisLeaf = *(int*)&buffer[currentIdxArray];
                    currentIdxArray += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&buffer[currentIdxArray]);
                    currentIdxArray += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }

                    ++nbLeafsInThisBuffer;
                }

                if( probStatus.MPI_SOURCE == rank - 1) {
                    hasToBeReceivedFromLeft -= nbLeafsInThisBuffer;
                }
                else {
                    hasToBeReceivedFromRight -= nbLeafsInThisBuffer;
                }
            }

            // I insert all my particles because I did not send anything
            FSize currentLeafPosition = 0;
            for(FSize idxLeaf = 0 ; idxLeaf < currentNbLeafs ; ++idxLeaf){
                const int nbPartInLeaf = (*(int*)&intervals[currentLeafPosition]);
                ParticleClass* const particles = reinterpret_cast<ParticleClass*>(&intervals[currentLeafPosition] + sizeof(int));

                for(int idxPart = 0 ; idxPart < nbPartInLeaf ; ++idxPart){
                    realTree.insert(particles[idxPart]);
                }
                currentLeafPosition += (nbPartInLeaf * sizeof(ParticleClass)) + sizeof(int);
            }

        }
        // if I may send to right and may receive from left
        else if( iWillReceiveFromLeft || iNeedToSendToRight ){
            printf("I may receive %d from left and send %d to right\n", hasToBeReceivedFromLeft, iNeedToSendRightCount);
            int hasBeenSentToRight = 0;
            MPI_Request requests[nbProcs];
            int iterRequest = 0;
            // If we do not have enought it means we can send all the buffer
            if(iDoNotHaveEnoughtToSendRight){
                FSize currentLeafPosition = 0;
                // Send all the data!
                for(int idxLeaf = 0 ; idxLeaf < currentNbLeafs ; ++idxLeaf){
                    currentLeafPosition += (*(int*)&particlesToSend[currentLeafPosition]* sizeof(ParticleClass)) + sizeof(int);
                }
                if(currentLeafPosition){
                    MPI_Isend( particlesToSend, currentLeafPosition, MPI_BYTE , rank + 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);
                }
                hasBeenSentToRight = currentNbLeafs;
            }
            else {
                // Else we have to send a part only
                FSize currentLeafPosition = 0;
                for(FSize idxLeaf = 0 ; idxLeaf < iNeedToSendRightCount && idxLeaf < iCanSendToRight ; ++idxLeaf){
                    currentLeafPosition += ((*(int*)&particlesToSend[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                hasBeenSentToRight = FMath::Min(iNeedToSendRightCount, iCanSendToRight);
                if(currentLeafPosition){
                    MPI_Isend(particlesToSend, int(currentLeafPosition), MPI_BYTE , rank + 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);
                }
                // and insert the other part into the tree
                for(FSize idxLeaf = hasBeenSentToRight ; idxLeaf < currentNbLeafs ; ++idxLeaf){
                    const int particlesInThisLeaf = *(int*)&particlesToSend[currentLeafPosition];
                    currentLeafPosition += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&particlesToSend[currentLeafPosition]);
                    currentLeafPosition += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }
                }
            }

            // Now we have to receive the data
            struct BufferDescriptor {
                char* buffer;
                int bufferSize;
            };
            int currentBufferIdx = 0;
            BufferDescriptor buffers[nbProcs];

            // While we have to receive some data
            while( hasToBeReceivedFromLeft ){
                // We must wait the message size
                MPI_Status probStatus;
                MPI_Probe( rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);

                BufferDescriptor currentMessage;
                MPI_Get_count( &probStatus,  MPI_BYTE, &currentMessage.bufferSize);
                currentMessage.buffer = new char[currentMessage.bufferSize];
                // Receive data
                MPI_Recv(currentMessage.buffer, currentMessage.bufferSize, MPI_BYTE, probStatus.MPI_SOURCE , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int nbLeafsInThisBuffer = 0;
                int currentIdxArray = 0;

                // If needed we transfer some data to the next proc
                if(hasBeenSentToRight < iNeedToSendRightCount){
                    FSize currentLeafPosition = 0;
                    while(currentLeafPosition < currentMessage.bufferSize && hasBeenSentToRight + nbLeafsInThisBuffer < iNeedToSendRightCount){
                        ++nbLeafsInThisBuffer;
                        currentLeafPosition += ((*(int*)&currentMessage.buffer[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                    }
                    hasBeenSentToRight += nbLeafsInThisBuffer;
                    MPI_Isend(currentMessage.buffer, int(currentLeafPosition), MPI_BYTE , rank + 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);

                    currentIdxArray = currentLeafPosition;
                }

                // Then we insert the other data
                while(currentIdxArray < currentMessage.bufferSize){
                    const int particlesInThisLeaf = *(int*)&currentMessage.buffer[currentIdxArray];
                    currentIdxArray += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&currentMessage.buffer[currentIdxArray]);
                    currentIdxArray += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }

                    ++nbLeafsInThisBuffer;
                }
                // We keep the buffer alive
                buffers[currentBufferIdx++] = currentMessage;
                hasToBeReceivedFromLeft -= nbLeafsInThisBuffer;
            }
            // Wait message to be sent
            MPI_Waitall( iterRequest, requests, MPI_STATUSES_IGNORE);
            for(int idxBuffer = 0 ; idxBuffer < currentBufferIdx ; ++idxBuffer){
                delete[] buffers[idxBuffer].buffer;
            }
        }
        // if I may send to left and may receive from right
        else if( iWillReceiveFromRight || iNeedToSendToLeft ){
            printf("I may send %d to left and receive %d from right\n", iNeedToSendLeftCount, hasToBeReceivedFromRight);
            int hasBeenSentToLeft = 0;
            MPI_Request requests[nbProcs];
            int iterRequest = 0;
            // If we do not have enought it means we can send all the buffer
            if(iDoNotHaveEnoughtToSendLeft){
                printf("iDoNotHaveEnoughtToSendLeft, first send %d\n", currentNbLeafs);
                FSize currentLeafPosition = 0;
                // Send all the data!
                for(int idxLeaf = 0 ; idxLeaf < currentNbLeafs ; ++idxLeaf){
                    currentLeafPosition += (*(int*)&particlesToSend[currentLeafPosition]* sizeof(ParticleClass)) + sizeof(int);
                }
                if(currentLeafPosition){
                    MPI_Isend( particlesToSend, currentLeafPosition, MPI_BYTE , rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);
                }
                hasBeenSentToLeft = currentNbLeafs;
            }
            else {
                // Else we have to send a part only
                FSize currentLeafPosition = 0;
                for(FSize idxLeaf = 0 ; idxLeaf < iNeedToSendLeftCount && idxLeaf < iCanSendToLeft ; ++idxLeaf){
                    currentLeafPosition += ((*(int*)&particlesToSend[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                hasBeenSentToLeft = FMath::Min(iNeedToSendLeftCount, iCanSendToLeft);
                if(currentLeafPosition){
                    MPI_Isend(particlesToSend, int(currentLeafPosition), MPI_BYTE , rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);
                }

                printf("I have enought to send left, first send %d\n", hasBeenSentToLeft);
                // and insert the other part into the tree
                for(FSize idxLeaf = hasBeenSentToLeft ; idxLeaf < currentNbLeafs ; ++idxLeaf){
                    const int particlesInThisLeaf = *(int*)&particlesToSend[currentLeafPosition];
                    currentLeafPosition += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&particlesToSend[currentLeafPosition]);
                    currentLeafPosition += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }
                }
            }

            // Now we have to receive the data
            struct BufferDescriptor {
                char* buffer;
                int bufferSize;
            };
            int currentBufferIdx = 0;
            BufferDescriptor buffers[nbProcs];

            // While we have to receive some data
            while( hasToBeReceivedFromRight ){
                printf("Go into loop...\n");
                // We must wait the message size
                MPI_Status probStatus;
                MPI_Probe( rank + 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &probStatus);

                BufferDescriptor currentMessage;
                MPI_Get_count( &probStatus,  MPI_BYTE, &currentMessage.bufferSize);
                currentMessage.buffer = new char[currentMessage.bufferSize];
                // Receive data
                MPI_Recv(currentMessage.buffer, currentMessage.bufferSize, MPI_BYTE, probStatus.MPI_SOURCE , FMpi::TagSandSettling , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("Received a message of size %d...\n", currentMessage.bufferSize);
                int nbLeafsInThisBuffer = 0;
                int currentIdxArray = 0;

                // If needed we transfer some data to the next proc
                if(hasBeenSentToLeft < iNeedToSendLeftCount){
                    FSize currentLeafPosition = 0;
                    while(currentLeafPosition < currentMessage.bufferSize && hasBeenSentToLeft + nbLeafsInThisBuffer < iNeedToSendLeftCount){
                        ++nbLeafsInThisBuffer;
                        currentLeafPosition += ((*(int*)&currentMessage.buffer[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                    }
                    hasBeenSentToLeft += nbLeafsInThisBuffer;
                    MPI_Isend(currentMessage.buffer, int(currentLeafPosition), MPI_BYTE , rank - 1, FMpi::TagSandSettling, MPI_COMM_WORLD, &requests[iterRequest++]);

                    currentIdxArray = currentLeafPosition;
                }

                // Then we insert the other data
                while(currentIdxArray < currentMessage.bufferSize){
                    const int particlesInThisLeaf = *(int*)&currentMessage.buffer[currentIdxArray];
                    currentIdxArray += sizeof(int);
                    ParticleClass*const particles = reinterpret_cast<ParticleClass*>(&currentMessage.buffer[currentIdxArray]);
                    currentIdxArray += sizeof(ParticleClass) * particlesInThisLeaf;

                    for( int idxPart = 0 ; idxPart < particlesInThisLeaf ; ++idxPart){
                        realTree.insert( particles[ idxPart ] );
                    }

                    ++nbLeafsInThisBuffer;
                }
                // We keep the buffer alive
                buffers[currentBufferIdx++] = currentMessage;
                hasToBeReceivedFromRight -= nbLeafsInThisBuffer;
                printf("nbLeafsInThisBuffer %d, hasToBeReceivedFromRight %d\n", nbLeafsInThisBuffer, hasToBeReceivedFromRight);
            }
            printf("Wait all\n");
            // Wait message to be sent
            MPI_Waitall( iterRequest, requests, MPI_STATUSES_IGNORE);
            printf("delete all\n");
            for(int idxBuffer = 0 ; idxBuffer < currentBufferIdx ; ++idxBuffer){
                delete[] buffers[idxBuffer].buffer;
            }
        }

        return true;
    }
      */

};

#endif // FMPITREEBUILDER_H
