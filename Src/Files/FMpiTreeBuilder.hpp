#ifndef FMPITREEBUILDER_H
#define FMPITREEBUILDER_H


#include "../Utils/FQuickSort.hpp"


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
        const float step = (float(inSize) / float(MpiGetNbProcs()));
        return T(FMath::Ceil(step * float(MpiGetRank())));
    }

    template< class T >
    static T GetRight(const T inSize) {
        const float step = (float(inSize) / float(MpiGetNbProcs()));
        const T res = T(FMath::Ceil(step * float(MpiGetRank()+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    template< class T >
    static T GetOtherRight(const T inSize, const int other) {
        const float step = (float(inSize) / MpiGetNbProcs());
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
    int nbLeavesInIntervals;

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
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();

        // First we create the particles that belong to us (our proc)
        // we compute their morton index to be able to sort them
        //

        IndexedParticle* outputArray = 0;
        long outputSize = 0;
        {
            // create particles
            IndexedParticle*const realParticlesIndexed = new IndexedParticle[loader.getNumberOfParticles()];
            memset(realParticlesIndexed, 0, sizeof(IndexedParticle)* loader.getNumberOfParticles());
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
            MortonIndex otherFirstIndex = -1;
            {
                FMpi::Request req[2];
                MPI_Status status[2];
                int reqiter = 0;
                // can I send my first index? == I am not rank 0 & I have data
                if( 0 < rank && outputSize){
                    MPI_Isend( &outputArray[0].index, 1, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, &req[reqiter++]);
                }
                if( rank != nbProcs - 1){
                    MPI_Irecv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD, &req[reqiter++]);
                }

                MPI_Waitall(reqiter,req,status);

                // I could not send because I do not have data, so I transmit the data coming
                // from my right neigbors
                if( 0 < rank && !outputSize){
                    MPI_Send( &otherFirstIndex, 1, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD);
                }
            }


            MPI_Request req[2];
            MPI_Status status[2];
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
                MPI_Probe(rank - 1, 0, MPI_COMM_WORLD, &probStatus);
                MPI_Get_count( &probStatus,  MPI_BYTE, &sendByOther);

                if(sendByOther){
                    sendByOther /= sizeof(IndexedParticle);
                    const IndexedParticle* const reallocOutputArray = outputArray;
                    const long reallocOutputSize = outputSize;

                    outputSize += sendByOther;
                    outputArray = new IndexedParticle[outputSize];
                    memcpy(&outputArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                    delete[] reallocOutputArray;

                    MPI_Recv(outputArray, sizeof(IndexedParticle) * sendByOther, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, &probStatus);
                }
                else{
                    MPI_Irecv(0, 0, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, &req[reqiter++]);
                }
            }

            if(rank != nbProcs - 1){

                long idxPart = outputSize - 1 ;
                while(idxPart >= 0 && outputArray[idxPart].index == otherFirstIndex){
                    --idxPart;
                }
                const long toSend = outputSize - 1 - idxPart;
                MPI_Isend( &outputArray[idxPart + 1], toSend * sizeof(IndexedParticle), MPI_BYTE, rank + 1, 0, MPI_COMM_WORLD, &req[reqiter++]);

                if( rank != 0 && !needToRecvBeforeSend && (rank != nbProcs - 1)){
                    int sendByOther = 0;

                    MPI_Status probStatus;
                    MPI_Probe(rank - 1, 0, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sendByOther);

                    if(sendByOther){
                        sendByOther /= sizeof(IndexedParticle);
                        char* const tempBuffer = new char[sizeof(IndexedParticle) * sendByOther];

                        MPI_Irecv(tempBuffer, sizeof(IndexedParticle) * sendByOther, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, &req[reqiter++]);

                        MPI_Waitall(reqiter,req, status);
                        reqiter = 0;

                        const IndexedParticle* const reallocOutputArray = outputArray;
                        const long reallocOutputSize = outputSize;

                        outputSize += sendByOther;
                        outputArray = new IndexedParticle[outputSize];
                        memcpy(&outputArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                        delete[] reallocOutputArray;
                        memcpy(outputArray, tempBuffer, sendByOther * sizeof(IndexedParticle));
                        delete[] tempBuffer;
                    }
                    else{
                        MPI_Irecv( 0, 0, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, &req[reqiter++]);
                    }
                }
            }
            MPI_Waitall(reqiter,req,status);
        }

        // We now copy the data from a sorted type into real particles array + counter

        nbLeavesInIntervals = 0;
        if(outputSize){
            intervals = new char[outputSize * (sizeof(ParticleClass) + sizeof(int))];

            MortonIndex previousIndex = -1;
            char* writeIndex = intervals;
            int* writeCounter = 0;

            for( int idxPart = 0; idxPart < outputSize ; ++idxPart){
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
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();

        FTic counter;
        counter.tic();
        //////////////////////////////////////////////////////////////////////////////////
        // We inform the master proc about the data we have
        //////////////////////////////////////////////////////////////////////////////////

        int nbLeafs = nbLeavesInIntervals;

        // receive from left and right
        if((rank == 0)){
            MPI_Send(&nbLeafs, sizeof(int), MPI_BYTE , 1, 0, MPI_COMM_WORLD);
        }
        else if(rank == nbProcs - 1){
            MPI_Send(&nbLeafs, sizeof(int), MPI_BYTE , rank - 1, 0, MPI_COMM_WORLD);
        }
        // receive
        int leftLeafs = 0;
        int rightLeafs = 0;
        if(!(rank == 0) && rank != nbProcs - 1){
            for(int idxToReceive = 0 ; idxToReceive < 2 ; ++idxToReceive){
                int source(0);
                int temp = 0;
                receiveDataFromTag(sizeof(int), 0, &temp, &source);
                if(source < rank){ // come from left
                    leftLeafs = temp;
                    temp += nbLeafs;
                    MPI_Send(&temp, sizeof(int), MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);
                }
                else { // come from right
                    rightLeafs = temp;
                    temp += nbLeafs;
                    MPI_Send(&temp, sizeof(int), MPI_BYTE , rank - 1, 0, MPI_COMM_WORLD);
                }
            }
        }
        else {
            if((rank == 0)){ // come from right
                receiveDataFromTag(sizeof(int), 0, &rightLeafs);
            }
            else { // come from left
                receiveDataFromTag(sizeof(int), 0, &leftLeafs);
            }
        }


        //////////////////////////////////////////////////////////////////////////////////
        // We balance the data
        //////////////////////////////////////////////////////////////////////////////////

        const int totalNbLeafs = (leftLeafs + nbLeafs + rightLeafs);
        const int myLeftLeaf = GetLeft(totalNbLeafs);
        const int myRightLeaf = GetRight(totalNbLeafs);

        const bool iNeedToSendToLeft = leftLeafs < myLeftLeaf;
        const bool iNeedToSendToRight = myRightLeaf < leftLeafs + nbLeafs;

        const bool iWillReceiveFromRight = leftLeafs + nbLeafs < myRightLeaf;
        const bool iWillReceiveFromLeft = leftLeafs > myLeftLeaf;

        const bool iDoNotHaveEnoughtToSendRight = myRightLeaf < leftLeafs;
        const bool iDoNotHaveEnoughtToSendLeft = leftLeafs + nbLeafs < myLeftLeaf;


        const int iNeedToSendLeftCount = myLeftLeaf - leftLeafs;
        const int iCanSendToLeft = nbLeafs;

        const int iNeedToSendRightCount = leftLeafs + nbLeafs - myRightLeaf;
        const int iCanSendToRight = nbLeafs;

        MPI_Request requests[2];
        MPI_Status status[2];
        int iterRequest = 0;

        int hasBeenSentToLeft = 0;
        int hasBeenSentToRight = 0;

        char* particlesToSend = 0;


        printf("on my left %d on my right %d\n",leftLeafs, rightLeafs);
        printf("iNeedToSendLeftCount %d iCanSendToLeft %d\n",iNeedToSendLeftCount, iCanSendToLeft);
        printf("iNeedToSendRightCount %d iCanSendToRight %d \n",iNeedToSendRightCount, iCanSendToRight);
        printf("Elapsed %lf\n", counter.tacAndElapsed());

        ///////////////////////////////
        // Manage data we already have
        ///////////////////////////////

        if(nbLeafs){
            particlesToSend = intervals;

            int currentLeafPosition = 0;

            //Send to Left (the first leaves
            if(iNeedToSendToLeft){
                for(int idxLeaf = 0 ; idxLeaf < iNeedToSendLeftCount && idxLeaf < iCanSendToLeft ; ++idxLeaf){
                    currentLeafPosition += ((*(int*)&particlesToSend[currentLeafPosition]) * sizeof(ParticleClass)) + sizeof(int);
                }
                hasBeenSentToLeft = FMath::Min(iNeedToSendLeftCount, iCanSendToLeft);
                MPI_Isend(particlesToSend, currentLeafPosition, MPI_BYTE , rank - 1, 0, MPI_COMM_WORLD, &requests[iterRequest++]);
                printf("I send to left %d bytes %d leaves\n", currentLeafPosition, hasBeenSentToLeft);
            }
            printf("Elapsed %lf\n", counter.tacAndElapsed());

            // Insert the particles I host and that belong to me
            const int beginForMe = (iNeedToSendToLeft ? FMath::Min(iNeedToSendLeftCount,iCanSendToLeft) : 0);
            const int endForMe = nbLeafs - (iNeedToSendToRight ? FMath::Min(iNeedToSendRightCount,iCanSendToRight) : 0);
            printf("I insert my particles from %d to %d \n", beginForMe, endForMe);
            for(int idxLeaf = beginForMe ; idxLeaf < endForMe ; ++idxLeaf){

                const int nbPartInLeaf = (*(int*)&particlesToSend[currentLeafPosition]);
                ParticleClass* const particles = reinterpret_cast<ParticleClass*>(&particlesToSend[currentLeafPosition] + sizeof(int));

                for(int idxPart = 0 ; idxPart < nbPartInLeaf ; ++idxPart){
                    realTree.insert(particles[idxPart]);
                }

                currentLeafPosition += (nbPartInLeaf * sizeof(ParticleClass)) + sizeof(int);
            }
            printf("Done\n");
            printf("Elapsed %lf\n", counter.tacAndElapsed());

            //Send to Right (the right-est leaves
            if(iNeedToSendToRight){
                const int beginWriteIndex = currentLeafPosition;

                for(int idxLeaf = 0 ; idxLeaf < iNeedToSendRightCount && idxLeaf < iCanSendToRight ; ++idxLeaf){
                    currentLeafPosition += (*(int*)&particlesToSend[currentLeafPosition]* sizeof(ParticleClass)) + sizeof(int);
                }

                hasBeenSentToRight = FMath::Min(iNeedToSendRightCount, iCanSendToRight);
                MPI_Isend( &particlesToSend[beginWriteIndex], currentLeafPosition - beginWriteIndex, MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD, &requests[iterRequest++]);
                printf("I send to right %d bytes %d leaves\n", currentLeafPosition - beginWriteIndex, hasBeenSentToRight);
            }
            printf("Elapsed %lf\n", counter.tacAndElapsed());
        }

        char* toRecvFromLeft = 0;
        char* toRecvFromRight = 0;
        int countReceive = int(iWillReceiveFromLeft) + int(iWillReceiveFromRight);
        int sizeOfLeftBuffer = 0;
        int sizeOfRightBuffer = 0;
        int sizeOfRightData = 0;
        int sizeOfLeftData = 0;

        int sourceToWhileRecv = MPI_ANY_SOURCE;

        printf("countReceive %d\n", countReceive);
        printf("Elapsed %lf\n", counter.tacAndElapsed());

        // Now prepare to receive data
        while(countReceive--){
            MPI_Status recvStatus;
            MPI_Probe(sourceToWhileRecv, 0, MPI_COMM_WORLD, &recvStatus);
            // receive from left
            if(recvStatus.MPI_SOURCE == rank - 1){
                MPI_Get_count( &recvStatus,  MPI_BYTE, &sizeOfLeftBuffer);
                toRecvFromLeft = new char[sizeOfLeftBuffer];
                sizeOfLeftData = sizeOfLeftBuffer;
                MPI_Irecv(toRecvFromLeft, sizeOfLeftBuffer, MPI_BYTE, rank - 1 , 0 , MPI_COMM_WORLD, &requests[iterRequest++]);
                sourceToWhileRecv = rank + 1;
                printf("I will receive %d bytes from left\n", sizeOfLeftBuffer);
            }
            // receive from right
            else{
                MPI_Get_count( &recvStatus,  MPI_BYTE, &sizeOfRightBuffer);
                toRecvFromRight = new char[sizeOfRightBuffer];
                sizeOfRightData = sizeOfRightBuffer;
                MPI_Irecv(toRecvFromRight, sizeOfRightBuffer, MPI_BYTE, rank + 1 , 0 , MPI_COMM_WORLD, &requests[iterRequest++]);
                sourceToWhileRecv = rank - 1;
                printf("I will receive %d bytes from right\n", sizeOfRightBuffer);
            }
        }

        ///////////////////////////////
        // Wait send receive
        ///////////////////////////////
        MPI_Waitall(iterRequest, requests, status);
        // We can delete the buffer use to send our particles only

        printf("Wait passed...\n");
        printf("Elapsed %lf\n", counter.tacAndElapsed());

        ///////////////////////////////
        // Process received data
        // and transfer if needed
        ///////////////////////////////
        // We have to receive from right and transfere to left
        int hasToBeReceivedFromLeft = leftLeafs - myLeftLeaf;
        int hasToBeReceivedFromRight = myRightLeaf - (leftLeafs + nbLeafs);
        int arrayIdxRight = 0;
        int arrayIdxLeft = 0;

        if(iDoNotHaveEnoughtToSendLeft){
            printf("iDoNotHaveEnoughtToSendLeft\n");
            do{
                arrayIdxRight = 0;
                while(arrayIdxRight < sizeOfRightData && hasBeenSentToLeft < iNeedToSendLeftCount){
                    const int particlesInThisLeaf = *(int*)&toRecvFromRight[arrayIdxRight];
                    arrayIdxRight += sizeof(int) + sizeof(ParticleClass) * particlesInThisLeaf;
                    --hasToBeReceivedFromRight;
                    ++hasBeenSentToLeft;
                }
                printf("Send %d to left, total leaves sent %d / %d\n", arrayIdxRight, hasBeenSentToLeft, iNeedToSendLeftCount);
                printf("Elapsed %lf\n", counter.tacAndElapsed());
                MPI_Send(toRecvFromRight, arrayIdxRight, MPI_BYTE , rank - 1, 0, MPI_COMM_WORLD);
                if(hasBeenSentToLeft < iNeedToSendLeftCount){
                    MPI_Status probStatus;
                    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfRightData);
                    if(sizeOfRightBuffer < sizeOfRightData){
                        sizeOfRightBuffer = sizeOfRightData;
                        delete[] toRecvFromRight;
                        toRecvFromRight = new char[sizeOfRightData];
                    }
                    printf("Receive %d bytes from right\n", sizeOfRightData);
                    MPI_Recv(toRecvFromRight, sizeOfRightData, MPI_BYTE, rank + 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasBeenSentToLeft < iNeedToSendLeftCount);
        }
        // We have to receive from left and transfere to right
        else if(iDoNotHaveEnoughtToSendRight){
            printf("iDoNotHaveEnoughtToSendRight\n");
            do{
                arrayIdxLeft = 0;
                while(arrayIdxLeft < sizeOfLeftData && hasBeenSentToRight < iNeedToSendRightCount){
                    const int particlesInThisLeaf = *(int*)&toRecvFromLeft[arrayIdxLeft];
                    arrayIdxLeft += sizeof(int) + sizeof(ParticleClass) * particlesInThisLeaf;
                    --hasToBeReceivedFromLeft;
                    ++hasBeenSentToRight;
                }
                printf("Send %d to right, total leaves sent %d / %d\n", arrayIdxLeft, hasBeenSentToRight, iNeedToSendRightCount);
                printf("Elapsed %lf\n", counter.tacAndElapsed());
                MPI_Send(toRecvFromLeft, arrayIdxLeft, MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);
                if(hasBeenSentToRight < iNeedToSendRightCount){
                    MPI_Status probStatus;
                    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfLeftData);
                    if(sizeOfLeftBuffer < sizeOfLeftData){
                        sizeOfLeftBuffer = sizeOfLeftData;
                        delete[] toRecvFromLeft;
                        toRecvFromLeft = new char[sizeOfLeftData];
                    }
                    printf("Receive %d bytes from left", sizeOfLeftData);
                    MPI_Recv(toRecvFromLeft, sizeOfLeftData, MPI_BYTE, rank - 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasBeenSentToRight < iNeedToSendRightCount);
        }

        printf("Finished to send\n");

        if(iWillReceiveFromLeft ){ // I need to wait
            printf("iWillReceiveFromLeft (hasToBeReceivedFromLeft %d leaves)\n", hasToBeReceivedFromLeft);
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

                printf("Remains hasToBeReceivedFromLeft %d leaves to receive\n",hasToBeReceivedFromLeft);
                printf("Elapsed %lf\n", counter.tacAndElapsed());

                if(hasToBeReceivedFromLeft){
                    MPI_Status probStatus;
                    MPI_Probe( rank - 1, 0, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfLeftData);
                    if(sizeOfLeftBuffer < sizeOfLeftData){
                        sizeOfLeftBuffer = sizeOfLeftData;
                        delete[] toRecvFromLeft;
                        toRecvFromLeft = new char[sizeOfLeftData];
                    }
                    printf("Received %d bytes from left\n", sizeOfLeftData);
                    MPI_Recv(toRecvFromLeft, sizeOfLeftData, MPI_BYTE, rank - 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasToBeReceivedFromLeft);
        }
        if(iWillReceiveFromRight){
            printf("iWillReceiveFromRight (hasToBeReceivedFromRight %d leaves)\n", hasToBeReceivedFromRight);
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

                printf("Remains hasToBeReceivedFromRight %d leaves to receive\n",hasToBeReceivedFromRight);
                printf("Elapsed %lf\n", counter.tacAndElapsed());

                if(hasToBeReceivedFromRight){
                    MPI_Status probStatus;
                    printf("Probe\n");
                    MPI_Probe( rank + 1, 0, MPI_COMM_WORLD, &probStatus);
                    MPI_Get_count( &probStatus,  MPI_BYTE, &sizeOfRightData);
                    printf("Receive %d bytes from right\n", sizeOfRightData);
                    if(sizeOfRightBuffer < sizeOfRightData){
                        sizeOfRightBuffer = sizeOfRightData;
                        delete[] toRecvFromRight;
                        toRecvFromRight = new char[sizeOfRightData];
                    }
                    MPI_Recv(toRecvFromRight, sizeOfRightData, MPI_BYTE, rank + 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasToBeReceivedFromRight);
        }

        // Delete reception buffers
        delete[] toRecvFromLeft;
        delete[] toRecvFromRight;

        return true;
    }
};

#endif // FMPITREEBUILDER_H
