#ifndef FMPITREEBUILDER_H
#define FMPITREEBUILDER_H


#include "../Utils/FQuickSort.hpp"

template<class ContainerClass, class ParticleClass>
class FMpiTreeBuilder{
    static long getTreeCoordinate(const FReal inRelativePosition, const FReal boxWidthAtLeafLevel) {
            const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
            const long index = FMath::dfloor(indexFReal);
            if( index && FMath::LookEqual(inRelativePosition, boxWidthAtLeafLevel * index ) ){
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

    static void receiveDataFromTag(const int inSize, const int inTag, void* const inData, int* const inSource = 0, int* const inFilledSize = 0){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, inTag, MPI_COMM_WORLD, &status);
        if(inSource) *inSource = status.MPI_SOURCE;
        if(inFilledSize) MPI_Get_count(&status,MPI_CHAR,inFilledSize);
    }

    template< class T >
    static T GetLeft(const T inSize) {
        const float step = (float(inSize) / MpiGetNbProcs());
        return T(FMath::Ceil(step * MpiGetRank()));
    }

    template< class T >
    static T GetRight(const T inSize) {
        const float step = (float(inSize) / MpiGetNbProcs());
        const T res = T(FMath::Ceil(step * (MpiGetRank()+1)));
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

    struct ParticlesGroup {
        int number;
        int positionInArray;
        MortonIndex index;
        ParticlesGroup(const int inNumber = 0 , const int inPositionInArray = 0, const MortonIndex inIndex = 0)
            : number(inNumber), positionInArray(inPositionInArray), index(inIndex) {
        }
    };


    struct IndexedParticle{
        MortonIndex index;
        ParticleClass particle;

        operator MortonIndex(){
            return this->index;
        }
    };

public:
    template <class OctreeClass, class LoaderClass>
    static bool SplitAndSortFile(OctreeClass& treeInterval, int& myNbParticlesCounter, LoaderClass& loader){
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();
        const int NbLevels = treeInterval.getHeight();

        IndexedParticle* outputArray = 0;
        long outputSize = 0;
        {
            // create particles
            IndexedParticle*const realParticlesIndexed = new IndexedParticle[loader.getNumberOfParticles()];
            F3DPosition boxCorner(loader.getCenterOfBox() - (loader.getBoxWidth()/2));
            FTreeCoordinate host;
            const FReal boxWidthAtLeafLevel = loader.getBoxWidth() / (1 << (NbLevels - 1) );
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
        {
            MortonIndex otherFirstIndex = -1;
            {
                FMpi::Request req[2];
                int reqiter = 0;
                if( 0 < rank && outputSize){
                    MPI_Isend( &outputArray[0].index, 1, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, &req[reqiter++]);
                }
                if( rank != nbProcs - 1){
                    MPI_Irecv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD, &req[reqiter++]);
                }

                MPI_Waitall(reqiter,req,0);

                if( 0 < rank && !outputSize){
                    MPI_Send( &otherFirstIndex, 1, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD);
                }
            }
            // at this point every one know the first index of his right neighbors
            const bool needToRecvBeforeSend = (rank != 0 && ((outputSize && outputArray[0].index == otherFirstIndex) || !outputSize));
            if( needToRecvBeforeSend || (rank == nbProcs - 1) ){
                long sendByOther = 0;
                MPI_Recv(&sendByOther, 1, MPI_LONG, rank - 1, 0, MPI_COMM_WORLD,0);
                if(sendByOther){
                    const IndexedParticle* const reallocOutputArray = outputArray;
                    const long reallocOutputSize = outputSize;

                    outputSize += sendByOther;
                    outputArray = new IndexedParticle[outputSize];
                    memcpy(&outputArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                    delete[] reallocOutputArray;

                    MPI_Recv(outputArray, sizeof(IndexedParticle) * sendByOther, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD, 0);
                }
            }
            if(rank != nbProcs - 1){
                long idxPart = outputSize - 1 ;
                while(idxPart >= 0 && outputArray[idxPart].index == otherFirstIndex){
                    --idxPart;
                }
                long toSend = outputSize - 1 - idxPart;
                MPI_Send( &toSend, 1, MPI_LONG, rank + 1, 0, MPI_COMM_WORLD);
                if(toSend){
                    MPI_Send( &outputArray[idxPart + 1], toSend * sizeof(IndexedParticle), MPI_BYTE, rank + 1, 0, MPI_COMM_WORLD);
                }
                if( rank != 0 && !needToRecvBeforeSend ){
                    long sendByOther = 0;
                    MPI_Recv(&sendByOther, 1, MPI_LONG, rank - 1, 0, MPI_COMM_WORLD, 0);
                    if(sendByOther){
                        const IndexedParticle* const reallocOutputArray = outputArray;
                        const long reallocOutputSize = outputSize;

                        outputSize += sendByOther;
                        outputArray = new IndexedParticle[outputSize];
                        memcpy(&outputArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                        delete[] reallocOutputArray;

                        MPI_Recv(outputArray, sizeof(IndexedParticle) * sendByOther, MPI_BYTE, rank - 1, 0, MPI_COMM_WORLD,0);
                    }
                }
            }
        }

        // we can insert into the tree
        for(int idxPart = 0 ; idxPart < outputSize ; ++idxPart){
            treeInterval.insert(outputArray[idxPart].particle);
        }
        myNbParticlesCounter = outputSize;

        return true;
    }

    template <class OctreeClass, class LoaderClass>
    static bool SplitAndSortFileWithoutQS(OctreeClass& treeInterval, int& myNbParticlesCounter, LoaderClass& loader){
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();

        //////////////////////////////////////////////////////////////////////////////////
        // My method
        //////////////////////////////////////////////////////////////////////////////////
        FVector<ParticlesGroup> groups;
        ParticleClass*const realParticles = reinterpret_cast<ParticleClass*>(new char[loader.getNumberOfParticles() * sizeof(ParticleClass)]);

        {
            OctreeClass sortingTree(treeInterval.getHeight(), treeInterval.getSubHeight() ,loader.getBoxWidth(),loader.getCenterOfBox());

            ParticleClass particle;
            for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(particle);
                sortingTree.insert(particle);
            }

            //////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////
            int indexPart = 0;

            typename OctreeClass::Iterator octreeIterator(&sortingTree);
            octreeIterator.gotoBottomLeft();
            do{
                typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
                const MortonIndex indexAtThisLeaf = octreeIterator.getCurrentGlobalIndex();

                groups.push(ParticlesGroup(octreeIterator.getCurrentListTargets()->getSize(),indexPart, indexAtThisLeaf));

                while( iter.hasNotFinished() ){
                    realParticles[indexPart] = iter.data();
                    ++indexPart;
                    iter.gotoNext();
                }
            } while(octreeIterator.moveRight());

        }

        //////////////////////////////////////////////////////////////////////////////////
        // We send the particle that do not belong to us
        //////////////////////////////////////////////////////////////////////////////////

        MortonIndex min = 0;
        MortonIndex max = 0;

        MPI_Reduce( &groups[0].index, &min, 1, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD );
        MPI_Reduce( &groups[groups.getSize() - 1].index, &max, 1, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD );

        MPI_Bcast ( &min, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
        MPI_Bcast ( &max, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );

        // Used for print information
        //const MortonIndex startIndex = GetLeft(max - min + 1) + min;
        //const MortonIndex endIndex = GetRight(max - min + 1) + min;

        int*const needToReceive = new int[nbProcs * nbProcs];
        memset(needToReceive,0,nbProcs * nbProcs * sizeof(int));

        FMpi::Request requests[nbProcs];
        {
            int needToSend[nbProcs];
            memset(needToSend, 0, sizeof(int) * nbProcs);

            MortonIndex rightMortonIndex = min;
            int groudIndex = 0;
            for(int idxProc = 0 ; idxProc < nbProcs && groudIndex < groups.getSize() ; ++idxProc){
                rightMortonIndex = GetOtherRight(max - min + 1, idxProc) + min;

                if(idxProc != rank){
                    int size = 0;
                    int currentGroupIndex = groudIndex;
                    while(groudIndex < groups.getSize() && groups[groudIndex].index < rightMortonIndex){
                        size += groups[groudIndex].number;
                        ++groudIndex;
                    }
                    needToSend[idxProc] = size;

                    MPI_Isend(&realParticles[groups[currentGroupIndex].positionInArray], sizeof(ParticleClass) * size, MPI_BYTE , idxProc, 1, MPI_COMM_WORLD, &requests[idxProc]);
                }
                else{
                    needToSend[idxProc] = 0;
                    while(groudIndex < groups.getSize() && groups[groudIndex].index < rightMortonIndex){
                        const int end = groups[groudIndex].positionInArray + groups[groudIndex].number;
                        for(int idxPart = groups[groudIndex].positionInArray ; idxPart < end ; ++idxPart){
                            //std::cout << "\t I keep (" << realParticles[idxPart].getPosition().getX() << ";" << realParticles[idxPart].getPosition().getY() << ";" << realParticles[idxPart].getPosition().getZ() << ")" << std::endl;
                            treeInterval.insert(realParticles[idxPart]);
                            ++myNbParticlesCounter;
                        }
                        ++groudIndex;
                    }
                }
            }

            MPI_Allgather( needToSend, nbProcs, MPI_INT, needToReceive, nbProcs, MPI_INT, MPI_COMM_WORLD);
        }


        //////////////////////////////////////////////////////////////////////////////////
        // We receive others particles and insert them in the tree
        //////////////////////////////////////////////////////////////////////////////////
        int CounterProcToReceive(0);
        int maxPartToReceive(0);
        for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
            if(idxProc != rank && needToReceive[nbProcs * idxProc + rank]){
                ++CounterProcToReceive;
                if(maxPartToReceive < needToReceive[nbProcs * idxProc + rank]){
                    maxPartToReceive = needToReceive[nbProcs * idxProc + rank];
                }
            }
        }


        ParticleClass*const iterParticles = reinterpret_cast<ParticleClass*>(new char[maxPartToReceive * sizeof(ParticleClass)]);
        // we receive message from nb proc - 1 (from every other procs
        for(int idxProc = 0 ; idxProc < CounterProcToReceive ; ++idxProc){
            MPI_Status status;
            MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            const int source = status.MPI_SOURCE;

            const int nbPartFromProc = needToReceive[nbProcs * source + rank];
            int received(0);

            MPI_Recv(iterParticles, sizeof(ParticleClass) * nbPartFromProc, MPI_BYTE, source, 1, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status,MPI_BYTE, &received);

            for(int idxPart = 0 ; idxPart < nbPartFromProc ; ++idxPart){
                //std::cout << "\t We receive a new particle (" << (*iterParticles).getPosition().getX() << ";" << (*iterParticles).getPosition().getY() << ";" << (*iterParticles).getPosition().getZ() << ")" << std::endl;
                treeInterval.insert(iterParticles[idxPart]);
                ++myNbParticlesCounter;
            }
        }

        for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
            if(idxProc != rank && needToReceive[nbProcs * rank + idxProc ]){
                MPI_Wait(&requests[idxProc], 0);
            }
        }

        delete [] reinterpret_cast<char*>(realParticles);
        delete [] needToReceive;

        return true;
    }


    template <class OctreeClass>
    static bool IntervalsToTree(OctreeClass& realTree, OctreeClass& treeInterval, const int myNbParticlesCounter){
        const int rank = MpiGetRank();
        const int nbProcs = MpiGetNbProcs();

        //////////////////////////////////////////////////////////////////////////////////
        // We inform the master proc about the data we have
        //////////////////////////////////////////////////////////////////////////////////

        int nbLeafs = 0;

        // we might now have any particles in our interval
        if(myNbParticlesCounter){
            typename OctreeClass::Iterator octreeIterator(&treeInterval);
            octreeIterator.gotoBottomLeft();
            do{
                ++nbLeafs;
            } while(octreeIterator.moveRight());
        }

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

        MPI_Request requests[10];
        int iterRequest = 0;

        int hasBeenSentToLeft = 0;
        int hasBeenSentToRight = 0;

        char* particlesToSend = 0;

        ///////////////////////////////
        // Manage data we already have
        ///////////////////////////////

        if(myNbParticlesCounter){
            particlesToSend = new char[sizeof(ParticleClass) * myNbParticlesCounter + sizeof(int) * nbLeafs];
            int idxWriteParticles = 0;

            typename OctreeClass::Iterator octreeIterator(&treeInterval);
            octreeIterator.gotoBottomLeft();

            //Send to Left (the first leaves
            if(iNeedToSendToLeft){
                for(int idxLeaf = 0 ; idxLeaf < iNeedToSendLeftCount && idxLeaf < iCanSendToLeft ; ++idxLeaf){
                    *(int*)&particlesToSend[idxWriteParticles] = octreeIterator.getCurrentListTargets()->getSize();
                    idxWriteParticles += sizeof(int);

                    memcpy(&particlesToSend[idxWriteParticles], octreeIterator.getCurrentListTargets()->data(), sizeof(ParticleClass) * octreeIterator.getCurrentListTargets()->getSize());
                    idxWriteParticles += sizeof(ParticleClass) * octreeIterator.getCurrentListTargets()->getSize();

                    octreeIterator.moveRight();
                }

                hasBeenSentToLeft = FMath::Min(iNeedToSendLeftCount, iCanSendToLeft);
                MPI_Send(particlesToSend, idxWriteParticles, MPI_BYTE , rank - 1, 0, MPI_COMM_WORLD);
            }

            // Insert the particles I host and that belong to me
            const int beginForMe = (iNeedToSendToLeft ? FMath::Min(iNeedToSendLeftCount,iCanSendToLeft) : 0);
            const int endForMe = nbLeafs - (iNeedToSendToRight ? FMath::Min(iNeedToSendRightCount,iCanSendToRight) : 0);
            for(int idxLeaf = beginForMe ; idxLeaf < endForMe ; ++idxLeaf){
                typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());

                while( iter.hasNotFinished() ){
                    realTree.insert(iter.data());
                    iter.gotoNext();
                }

                octreeIterator.moveRight();
            }

            //Send to Right (the right-est leaves
            if(iNeedToSendToRight){
                const int beginWriteIndex = idxWriteParticles;

                for(int idxLeaf = 0 ; idxLeaf < iNeedToSendRightCount && idxLeaf < iCanSendToRight ; ++idxLeaf){
                    *(int*)&particlesToSend[idxWriteParticles] = octreeIterator.getCurrentListTargets()->getSize();
                    idxWriteParticles += sizeof(int);

                    memcpy(&particlesToSend[idxWriteParticles], octreeIterator.getCurrentListTargets()->data(), sizeof(ParticleClass) * octreeIterator.getCurrentListTargets()->getSize());
                    idxWriteParticles += sizeof(ParticleClass) * octreeIterator.getCurrentListTargets()->getSize();

                    octreeIterator.moveRight();
                }

                hasBeenSentToRight = FMath::Min(iNeedToSendRightCount, iCanSendToRight);
                MPI_Isend( &particlesToSend[beginWriteIndex], idxWriteParticles - beginWriteIndex, MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD, &requests[iterRequest++]);
            }
        }

        char* toRecvFromLeft = 0;
        char* toRecvFromRight = 0;
        int countReceive = int(iWillReceiveFromLeft) + int(iWillReceiveFromRight);
        int bytesRecvFromLeft = 0;
        int bytesRecvFromRight = 0;

        // Now prepare to receive data
        while(countReceive--){
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            // receive from left
            if(status.MPI_SOURCE == rank - 1){
                MPI_Get_count( &status,  MPI_BYTE, &bytesRecvFromLeft);
                toRecvFromLeft = new char[bytesRecvFromLeft];
                MPI_Irecv(toRecvFromLeft, bytesRecvFromLeft, MPI_BYTE, rank - 1 , 0 , MPI_COMM_WORLD, &requests[iterRequest++]);
            }
            // receive from right
            else{
                MPI_Get_count( &status,  MPI_BYTE, &bytesRecvFromRight);
                toRecvFromRight = new char[bytesRecvFromRight];
                MPI_Irecv(toRecvFromRight, bytesRecvFromRight, MPI_BYTE, rank + 1 , 0 , MPI_COMM_WORLD, &requests[iterRequest++]);
            }
        }

        ///////////////////////////////
        // Wait send receive
        ///////////////////////////////
        MPI_Waitall(iterRequest, requests, 0);
        // We can delete the buffer use to send our particles only
        delete[] particlesToSend;
        particlesToSend = 0;

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
            do{
                arrayIdxRight = 0;
                while(arrayIdxRight < bytesRecvFromRight && hasBeenSentToLeft < iNeedToSendLeftCount){
                    const int particlesInThisLeaf = *(int*)&toRecvFromRight[arrayIdxRight];
                    arrayIdxRight += sizeof(int) + sizeof(ParticleClass) * particlesInThisLeaf;
                    --hasToBeReceivedFromRight;
                }
                MPI_Send(toRecvFromRight, arrayIdxRight, MPI_BYTE , rank - 1, 0, MPI_COMM_WORLD);
                if(hasBeenSentToLeft < iNeedToSendLeftCount){
                    MPI_Status status;
                    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                    int status_size = 0;
                    MPI_Get_count( &status,  MPI_BYTE, &status_size);
                    if(bytesRecvFromRight < status_size){
                        bytesRecvFromRight = status_size;
                        delete[] toRecvFromRight;
                        toRecvFromRight = new char[status_size];
                    }
                    MPI_Recv(toRecvFromRight, status_size, MPI_BYTE, rank + 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasBeenSentToLeft < iNeedToSendLeftCount);
        }
        // We have to receive from left and transfere to right
        else if(iDoNotHaveEnoughtToSendRight){
            do{
                arrayIdxLeft = 0;
                while(arrayIdxLeft < bytesRecvFromLeft && hasBeenSentToRight < iNeedToSendRightCount){
                    const int particlesInThisLeaf = *(int*)&toRecvFromLeft[arrayIdxLeft];
                    arrayIdxLeft += sizeof(int) + sizeof(ParticleClass) * particlesInThisLeaf;
                    --hasToBeReceivedFromLeft;
                }
                MPI_Send(toRecvFromLeft, arrayIdxLeft, MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);
                if(hasBeenSentToRight < iNeedToSendRightCount){
                    MPI_Status status;
                    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                    int status_size = 0;
                    MPI_Get_count( &status,  MPI_BYTE, &status_size);
                    if(bytesRecvFromLeft < status_size){
                        bytesRecvFromLeft = status_size;
                        delete[] toRecvFromLeft;
                        toRecvFromLeft = new char[status_size];
                    }
                    MPI_Recv(toRecvFromLeft, status_size, MPI_BYTE, rank - 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasBeenSentToRight < iNeedToSendRightCount);
        }

        if(iWillReceiveFromLeft ){ // I need to wait
            do{
                while(arrayIdxLeft < bytesRecvFromLeft){
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
                    MPI_Status status;
                    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                    int status_size = 0;
                    MPI_Get_count( &status,  MPI_BYTE, &status_size);
                    if(bytesRecvFromLeft < status_size){
                        bytesRecvFromLeft = status_size;
                        delete[] toRecvFromLeft;
                        toRecvFromLeft = new char[status_size];
                    }
                    MPI_Recv(toRecvFromLeft, status_size, MPI_BYTE, rank - 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            } while(hasToBeReceivedFromLeft);
        }
        if(iWillReceiveFromRight){
            do{
                while(arrayIdxRight < bytesRecvFromRight){
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
                    MPI_Status status;
                    MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                    int status_size = 0;
                    MPI_Get_count( &status,  MPI_BYTE, &status_size);
                    if(bytesRecvFromRight < status_size){
                        bytesRecvFromRight = status_size;
                        delete[] toRecvFromRight;
                        toRecvFromRight = new char[status_size];
                    }
                    MPI_Recv(toRecvFromRight, status_size, MPI_BYTE, rank + 1 , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
