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

        const MortonIndex startIndex = GetLeft(max - min + 1) + min;
        const MortonIndex endIndex = GetRight(max - min + 1) + min;

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

        FVector<ParticlesGroup> groups;
        ParticleClass*const realParticles = myNbParticlesCounter?new ParticleClass[myNbParticlesCounter]:0;

        int nbLeafs = 0;

        // we might now have any particles in our interval
        if(myNbParticlesCounter){
            typename OctreeClass::Iterator octreeIterator(&treeInterval);
            octreeIterator.gotoBottomLeft();
            int indexPart = 0;
            do{
                typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
                const MortonIndex indexAtThisLeaf = octreeIterator.getCurrentGlobalIndex();

                groups.push(ParticlesGroup(octreeIterator.getCurrentListTargets()->getSize(),indexPart, indexAtThisLeaf));

                while( iter.hasNotFinished() ){
                    realParticles[indexPart] = iter.data();
                    ++indexPart;
                    iter.gotoNext();
                }

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

        ParticleClass* rpart(0);
        int rpartSize(0);

        // Do I need to send to right?
        if(iNeedToSendToRight){
            int iNeedToSend = leftLeafs + nbLeafs - myRightLeaf;
            int iCanSend = nbLeafs;
            int idxSend (0);

            MPI_Send(&iNeedToSend, sizeof(int), MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);

            while(idxSend < iNeedToSend && idxSend < iCanSend){
                MPI_Send(&groups[nbLeafs - idxSend - 1].number, sizeof(int), MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);
                MPI_Send(&realParticles[groups[nbLeafs - idxSend - 1].positionInArray], sizeof(ParticleClass) * groups[nbLeafs - idxSend - 1].number,
                         MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);

                ++idxSend;
            }
            // I need to wait (idxSend == iCanSend && idxSend < iNeedToSend)
            if( iDoNotHaveEnoughtToSendRight ){ // I need to wait
                int nbLeafsToRead(0);
                receiveDataFromTag(sizeof(int), 0, &nbLeafsToRead);
                for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                    int nbPartToRead(0);
                    receiveDataFromTag(sizeof(int), 0, &nbPartToRead);

                    if(rpartSize < nbPartToRead){
                        rpartSize = nbPartToRead;
                        delete [] (reinterpret_cast<char*>(rpart));
                        rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                    }

                    receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 0, rpart);
                    if(idxSend < iNeedToSend){
                        MPI_Send(&nbPartToRead, sizeof(int), MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);
                        MPI_Send(rpart, sizeof(ParticleClass) * nbPartToRead, MPI_BYTE , rank + 1, 0, MPI_COMM_WORLD);

                        ++idxSend;
                    }
                    else{
                        //insert into tree
                        for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                            realTree.insert(rpart[idxPart]);
                        }
                    }
                }
            }
        }
        // will I receive from left
        if(iNeedToSendToLeft){
            int iNeedToSend = myLeftLeaf - leftLeafs;
            int iCanSend = nbLeafs;
            int idxSend (0);

            MPI_Send(&iNeedToSend, sizeof(int), MPI_BYTE , rank - 1, 1, MPI_COMM_WORLD);

            while(idxSend < iNeedToSend && idxSend < iCanSend){
                MPI_Send(&groups[idxSend].number, sizeof(int), MPI_BYTE , rank - 1, 1, MPI_COMM_WORLD);
                MPI_Send(&realParticles[groups[idxSend].positionInArray], sizeof(ParticleClass) * groups[idxSend].number,
                         MPI_BYTE , rank - 1, 1, MPI_COMM_WORLD);

                ++idxSend;
            }
            // Can I do it now?
            if( iDoNotHaveEnoughtToSendLeft ){
                int nbLeafsToRead(0);
                receiveDataFromTag(sizeof(int), 1, &nbLeafsToRead);
                for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                    int nbPartToRead(0);
                    receiveDataFromTag(sizeof(int), 1, &nbPartToRead);

                    if(rpartSize < nbPartToRead){
                        rpartSize = nbPartToRead;
                        delete [] (reinterpret_cast<char*>(rpart));
                        rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                    }

                    receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 1, rpart);
                    if(idxSend < iNeedToSend){
                        MPI_Send(&nbPartToRead, sizeof(int), MPI_BYTE , rank - 1, 1, MPI_COMM_WORLD);
                        MPI_Send(rpart, sizeof(ParticleClass) * nbPartToRead, MPI_BYTE , rank - 1, 1, MPI_COMM_WORLD);

                        ++idxSend;
                    }
                    else{
                        for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                            realTree.insert(rpart[idxPart]);
                        }
                    }
                }
            }
        }

        // If i will receive from left and I did no already have
        if(!(iNeedToSendToRight && iDoNotHaveEnoughtToSendRight) && iWillReceiveFromLeft){
            int nbLeafsToRead(0);
            receiveDataFromTag(sizeof(int), 0, &nbLeafsToRead);
            for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                int nbPartToRead(0);
                receiveDataFromTag(sizeof(int), 0, &nbPartToRead);
                //printf("%d I will receive %d particles\n",rank, nbPartToRead);
                if(rpartSize < nbPartToRead){
                    rpartSize = nbPartToRead;
                    delete [] (reinterpret_cast<char*>(rpart));
                    rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                }

                receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 0, rpart);
                for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                    realTree.insert(rpart[idxPart]);
                }
            }
        }
        // If i will receive from right and I did no already have
        if(!(iNeedToSendToLeft && iDoNotHaveEnoughtToSendLeft) && iWillReceiveFromRight){
            int nbLeafsToRead(0);
            receiveDataFromTag(sizeof(int), 1, &nbLeafsToRead);
            for(int idxToRead = 0 ; idxToRead < nbLeafsToRead ; ++idxToRead){
                int nbPartToRead(0);
                receiveDataFromTag(sizeof(int), 1, &nbPartToRead);
                //printf("%d I will receive %d particles\n",rank, nbPartToRead);
                if(rpartSize < nbPartToRead){
                    rpartSize = nbPartToRead;
                    delete [] (reinterpret_cast<char*>(rpart));
                    rpart = reinterpret_cast<ParticleClass*>(new char[nbPartToRead*sizeof(ParticleClass)]);
                }

                receiveDataFromTag(nbPartToRead*sizeof(ParticleClass), 1, rpart);
                for(int idxPart = 0 ; idxPart < nbPartToRead ; ++idxPart){
                    realTree.insert(rpart[idxPart]);
                }
            }
        }

        // insert the particles we already have
        if(leftLeafs != totalNbLeafs){
            for(int idxLeafInsert = FMath::Max(myLeftLeaf-leftLeafs,0) ; idxLeafInsert <  totalNbLeafs - rightLeafs - leftLeafs - FMath::Max(0,leftLeafs + nbLeafs - myRightLeaf) ; ++idxLeafInsert){
                for(int idxPart = 0 ; idxPart < groups[idxLeafInsert].number ; ++idxPart){
                    realTree.insert(realParticles[groups[idxLeafInsert].positionInArray + idxPart]);
                }
            }
        }

        delete [] reinterpret_cast<char*>(rpart);

        return true;
    }
};

#endif // FMPITREEBUILDER_H
