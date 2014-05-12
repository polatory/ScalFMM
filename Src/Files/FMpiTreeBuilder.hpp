// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FMPITREEBUILDER_H
#define FMPITREEBUILDER_H

#include "../Utils/FMpi.hpp"
#include "../Utils/FQuickSortMpi.hpp"
#include "../Utils/FBitonicSort.hpp"

#include "../Utils/FMemUtils.hpp"
#include "../Utils/FTrace.hpp"

#include "../BalanceTree/FLeafBalance.hpp"

/** This class manage the loading of particles for the mpi version.
 * It use a BinLoader and then sort the data with a parallel quick sort
 * or bitonic sort.
 * After it needs to merge the leaves to finally equalize the data.
 */
template<class ParticleClass>
class FMpiTreeBuilder{
public:
    /** What sorting algorithm to use */
    enum SortingType{
        QuickSort,
        BitonicSort,
    };

private:
    /** This method has been tacken from the octree
   * it computes a tree coordinate (x or y or z) from real position
   */
    static int getTreeCoordinate(const FReal inRelativePosition, const FReal boxWidthAtLeafLevel) {
        const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
        const int index = int(FMath::dfloor(indexFReal));
        if( index && FMath::LookEqual(inRelativePosition, boxWidthAtLeafLevel * FReal(index) ) ){
            return index - 1;
        }
        return index;
    }


    /** A particle may not have a MortonIndex Method
   * but they are sorted on their morton index
   * so this struct store a particle + its index
   */
    struct IndexedParticle{
        MortonIndex index;
        ParticleClass particle;

        operator MortonIndex(){
            return this->index;
        }
    };

    //////////////////////////////////////////////////////////////////////////
    // To sort tha particles we hold
    //////////////////////////////////////////////////////////////////////////


    template <class LoaderClass>
    static void SortParticles( const FMpi::FComm& communicator, LoaderClass& loader, const SortingType type,
                               const int TreeHeight, IndexedParticle**const outputArray, FSize* const outputSize){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__ , "Loader to Tree" , __FILE__ , __LINE__) );

        // create particles
        IndexedParticle*const realParticlesIndexed = new IndexedParticle[loader.getNumberOfParticles()];
        FMemUtils::memset(realParticlesIndexed, 0, sizeof(IndexedParticle) * loader.getNumberOfParticles());
        FPoint boxCorner(loader.getCenterOfBox() - (loader.getBoxWidth()/2));
        FTreeCoordinate host;

        const FReal boxWidthAtLeafLevel = loader.getBoxWidth() / FReal(1 << (TreeHeight - 1) );
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(realParticlesIndexed[idxPart].particle);
            host.setX( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
            host.setY( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
            host.setZ( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));

            realParticlesIndexed[idxPart].index = host.getMortonIndex(TreeHeight - 1);
        }

        // sort particles
        if(type == QuickSort){
            FQuickSortMpi<IndexedParticle,MortonIndex, FSize>::QsMpi(realParticlesIndexed, loader.getNumberOfParticles(), *outputArray, *outputSize,communicator);
            delete [] (realParticlesIndexed);
        }
        else {
            FBitonicSort<IndexedParticle,MortonIndex, FSize>::Sort( realParticlesIndexed, loader.getNumberOfParticles(), communicator );
            *outputArray = realParticlesIndexed;
            *outputSize = loader.getNumberOfParticles();
        }
    }

    static void SortParticlesFromArray( const FMpi::FComm& communicator, const ParticleClass array[], const FSize size, const SortingType type,
                                        const FPoint& centerOfBox, const FReal boxWidth,
                                        const int TreeHeight, IndexedParticle**const outputArray, FSize* const outputSize){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__ , "Loader to Tree" , __FILE__ , __LINE__) );

        // create particles
        IndexedParticle*const realParticlesIndexed = new IndexedParticle[size];
        FMemUtils::memset(realParticlesIndexed, 0, sizeof(IndexedParticle) * size);

        FPoint boxCorner(centerOfBox - (boxWidth/2));
        FTreeCoordinate host;

        const FReal boxWidthAtLeafLevel = boxWidth / FReal(1 << (TreeHeight - 1) );

        for(int idxPart = 0 ; idxPart < size ; ++idxPart){
            realParticlesIndexed[idxPart].particle = array[idxPart];
            host.setX( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
            host.setY( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
            host.setZ( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));

            realParticlesIndexed[idxPart].index = host.getMortonIndex(TreeHeight - 1);
        }

        // sort particles
        if(type == QuickSort){
            FQuickSortMpi<IndexedParticle,MortonIndex, FSize>::QsMpi(realParticlesIndexed, size, *outputArray, *outputSize,communicator);
            delete [] (realParticlesIndexed);
        }
        else {
            FBitonicSort<IndexedParticle,MortonIndex, FSize>::Sort( realParticlesIndexed, size, communicator );
            *outputArray = realParticlesIndexed;
            *outputSize = size;
        }
    }


    //////////////////////////////////////////////////////////////////////////
    // To merge the leaves
    //////////////////////////////////////////////////////////////////////////

    static void MergeLeaves(const FMpi::FComm& communicator, IndexedParticle*& workingArray, FSize& workingSize,
                            char** leavesArray, FSize* const leavesSize){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader to Tree" , __FILE__ , __LINE__) );
        const int rank = communicator.processId();
        const int nbProcs = communicator.processCount();

        // be sure there is no splited leaves
        // to do that we exchange the first index with the left proc
        {
            FTRACE( FTrace::FRegion regionTrace("Remove Splited leaves", __FUNCTION__ , __FILE__ , __LINE__) );

            MortonIndex otherFirstIndex = -1;
            if(workingSize != 0 && rank != 0 && rank != nbProcs - 1){
                MPI_Sendrecv(&workingArray[0].index, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs,
                        &otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs,
                        communicator.getComm(), MPI_STATUS_IGNORE);
            }
            else if( rank == 0){
                MPI_Recv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs, communicator.getComm(), MPI_STATUS_IGNORE);
            }
            else if( rank == nbProcs - 1){
                MPI_Send( &workingArray[0].index, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs, communicator.getComm());
            }
            else {
                MPI_Recv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs, communicator.getComm(), MPI_STATUS_IGNORE);
                MPI_Send(&otherFirstIndex, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs, communicator.getComm());
            }

            // at this point every one know the first index of his right neighbors
            const bool needToRecvBeforeSend = (rank != 0 && ((workingSize && otherFirstIndex == workingArray[0].index ) || !workingSize));
            MPI_Request requestSendLeaf;

            IndexedParticle* sendBuffer = 0;
            if(rank != nbProcs - 1 && needToRecvBeforeSend == false){
                FSize idxPart = workingSize - 1 ;
                while(idxPart >= 0 && workingArray[idxPart].index == otherFirstIndex){
                    --idxPart;
                }
                const int particlesToSend = int(workingSize - 1 - idxPart);
                if(particlesToSend){
                    workingSize -= particlesToSend;
                    sendBuffer = new IndexedParticle[particlesToSend];
                    memcpy(sendBuffer, &workingArray[idxPart + 1], particlesToSend * sizeof(IndexedParticle));

                    MPI_Isend( sendBuffer, particlesToSend * int(sizeof(IndexedParticle)), MPI_BYTE,
                               rank + 1, FMpi::TagSplittedLeaf, communicator.getComm(), &requestSendLeaf);
                }
                else{
                    MPI_Isend( 0, 0, MPI_BYTE, rank + 1, FMpi::TagSplittedLeaf, communicator.getComm(), &requestSendLeaf);
                }
            }

            if( rank != 0 ){
                int sendByOther = 0;

                MPI_Status probStatus;
                MPI_Probe(rank - 1, FMpi::TagSplittedLeaf, communicator.getComm(), &probStatus);
                MPI_Get_count( &probStatus,  MPI_BYTE, &sendByOther);

                if(sendByOther){
                    sendByOther /= int(sizeof(IndexedParticle));

                    const IndexedParticle* const reallocOutputArray = workingArray;
                    const FSize reallocOutputSize = workingSize;

                    workingSize += sendByOther;
                    workingArray = new IndexedParticle[workingSize];
                    FMemUtils::memcpy(&workingArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
                    delete[] reallocOutputArray;

                    MPI_Recv(workingArray, int(sizeof(IndexedParticle)) * sendByOther, MPI_BYTE,
                             rank - 1, FMpi::TagSplittedLeaf, communicator.getComm(), MPI_STATUS_IGNORE);
                }
                else{
                    MPI_Recv( 0, 0, MPI_BYTE, rank - 1, FMpi::TagSplittedLeaf, communicator.getComm(), MPI_STATUS_IGNORE);
                }
            }

            if(rank != nbProcs - 1 && needToRecvBeforeSend == true){
                MPI_Send( workingArray, int(workingSize * sizeof(IndexedParticle)), MPI_BYTE,
                          rank + 1, FMpi::TagSplittedLeaf, communicator.getComm());
                delete[] workingArray;
                workingArray = 0;
                workingSize  = 0;
            }
            else if(rank != nbProcs - 1){
                MPI_Wait( &requestSendLeaf, MPI_STATUS_IGNORE);
                delete[] sendBuffer;
                sendBuffer = 0;
            }
        }

        {
            FTRACE( FTrace::FRegion regionTrace("Remove Splited leaves", __FUNCTION__ , __FILE__ , __LINE__) );
            // We now copy the data from a sorted type into real particles array + counter

            (*leavesSize)  = 0;
            (*leavesArray) = 0;

            if(workingSize){
                (*leavesArray) = new char[workingSize * (sizeof(ParticleClass) + sizeof(int))];

                MortonIndex previousIndex = -1;
                char* writeIndex = (*leavesArray);
                int* writeCounter = 0;

                for( FSize idxPart = 0; idxPart < workingSize ; ++idxPart){
                    if( workingArray[idxPart].index != previousIndex ){
                        previousIndex = workingArray[idxPart].index;
                        ++(*leavesSize);

                        writeCounter = reinterpret_cast<int*>( writeIndex );
                        writeIndex += sizeof(int);

                        (*writeCounter) = 0;
                    }

                    memcpy(writeIndex, &workingArray[idxPart].particle, sizeof(ParticleClass));

                    writeIndex += sizeof(ParticleClass);
                    ++(*writeCounter);
                }
            }

            delete [] workingArray;

            workingArray = 0;
            workingSize = 0;
        }
    }

    //////////////////////////////////////////////////////////////////////////
    // To equalize (same number of leaves among the procs)
    //////////////////////////////////////////////////////////////////////////

    /** Put the interval into a tree */
    template <class ContainerClass>
    static void EqualizeAndFillTree(const FMpi::FComm& communicator,  ContainerClass* particlesSaver,
                                    char*& leavesArray, FSize& nbLeavesInIntervals, FAbstractBalanceAlgorithm * balancer){


        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader to Tree" , __FILE__ , __LINE__) );
        const int myRank = communicator.processId();
        const int nbProcs = communicator.processCount();
        const FSize myNumberOfLeaves = nbLeavesInIntervals;

        // We have to know the number of leaves each procs holds
        int*const numberOfLeavesPerProc = new int[nbProcs];
        memset(numberOfLeavesPerProc, 0, sizeof(int) * nbProcs);
        int intNbLeavesInIntervals = int(nbLeavesInIntervals);
        MPI_Allgather(&intNbLeavesInIntervals, 1, MPI_INT, numberOfLeavesPerProc, 1, MPI_INT, communicator.getComm());
        printf("Proc : %d : Currently have %lld leaves \n", myRank, myNumberOfLeaves);

        //Start working THERE

        //We need the max number of leafs over th procs
        int*const leavesOffsetPerProc = new int[nbProcs + 1];
        FSize totalNumberOfLeaves  = 0;

        leavesOffsetPerProc[0] = 0;
        totalNumberOfLeaves   += numberOfLeavesPerProc[0];
        for(int idxProc = 1 ; idxProc < nbProcs ; ++idxProc ){
            leavesOffsetPerProc[idxProc] = leavesOffsetPerProc[idxProc-1] + numberOfLeavesPerProc[idxProc-1];
            totalNumberOfLeaves += numberOfLeavesPerProc[idxProc];
        }
        leavesOffsetPerProc[nbProcs] = int(totalNumberOfLeaves);

        const FSize currentLeafsOnMyLeft  = leavesOffsetPerProc[myRank];
        const FSize currentRightLeafIdx   = leavesOffsetPerProc[myRank+1];

        //Creation of an array to store how many parts are in each leaf
        int*const numberOfParticlesPerLeaf = new int[totalNumberOfLeaves];
        int*const myParticlesCounterArray  = &numberOfParticlesPerLeaf[leavesOffsetPerProc[myRank]];

        memset(numberOfParticlesPerLeaf, 0, sizeof(int)*totalNumberOfLeaves);

        //Loop over leafArray to fill myParts
        size_t idxOfParticlesNumber = 0;
        for(int idxLeaf = 0 ; idxLeaf < nbLeavesInIntervals ; ++idxLeaf){
            const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&leavesArray[idxOfParticlesNumber]));
            myParticlesCounterArray[idxLeaf] += numberOfParticlesInThisLeaf;
            idxOfParticlesNumber += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
        }

        MPI_Allgatherv(myParticlesCounterArray,numberOfLeavesPerProc[myRank], MPI_INT,
                       numberOfParticlesPerLeaf, numberOfLeavesPerProc, leavesOffsetPerProc, MPI_INT, communicator.getComm());

        FSize totalNumberOfParticles = 0;
        for(int idxLeaf = 0 ; idxLeaf < totalNumberOfLeaves ; ++idxLeaf){
            totalNumberOfParticles += numberOfParticlesPerLeaf[idxLeaf];
        }

        const FSize correctLeftLeavesNumber     = balancer->getLeft( totalNumberOfLeaves,numberOfParticlesPerLeaf,totalNumberOfParticles,
                                                                     NULL,nbProcs,myRank);
        const FSize correctRightLeavesIndex     = balancer->getRight(totalNumberOfLeaves,numberOfParticlesPerLeaf,totalNumberOfParticles,
                                                                     NULL,nbProcs,myRank);

        //// TODO REMOVE WHEN DEBUG printf("Proc [%d] :: will work from leaf %lld \t to leaf %lld \n",myRank,correctLeftLeavesNumber,correctRightLeavesIndex);

        MPI_Request* requests = new MPI_Request[nbProcs * 2];
        int counterRequest = 0;

        if(currentLeafsOnMyLeft < correctLeftLeavesNumber || correctRightLeavesIndex < currentRightLeafIdx){
            size_t offsetLeafToSend = 0;
            int counterLeafToSend = 0;
            int idxProcToProceed  = 0;

            while( idxProcToProceed < nbProcs && (balancer->getLeft(totalNumberOfLeaves,numberOfParticlesPerLeaf,totalNumberOfParticles,NULL,nbProcs,idxProcToProceed) < currentRightLeafIdx)){
                const FSize procToProceedRightIdx = balancer->getRight(totalNumberOfLeaves,numberOfParticlesPerLeaf,totalNumberOfParticles,NULL,nbProcs,idxProcToProceed);
                const FSize procToProceedLeftIdx  = balancer->getLeft(totalNumberOfLeaves,numberOfParticlesPerLeaf,totalNumberOfParticles,NULL,nbProcs,idxProcToProceed);
                const bool procToProceedHasLeftInMyInterval  = (currentLeafsOnMyLeft <= procToProceedLeftIdx && procToProceedLeftIdx < currentRightLeafIdx);
                const bool procToProceedHasRightInMyInterval = (currentLeafsOnMyLeft <= procToProceedRightIdx && procToProceedRightIdx < currentRightLeafIdx);
                const bool procIncludeMyInterval = (procToProceedLeftIdx <= currentLeafsOnMyLeft  && currentRightLeafIdx <= procToProceedRightIdx);
                //// TODO REMOVE WHEN DEBUG printf("%d] idxProcToProceed %d procToProceedRightIdx %llu procToProceedLeftIdx %llu procToProceedHasLeftInMyInterval %d procToProceedHasRightInMyInterval %d\n",
                //// TODO REMOVE WHEN DEBUG        myRank, idxProcToProceed, procToProceedRightIdx, procToProceedLeftIdx, procToProceedHasLeftInMyInterval, procToProceedHasRightInMyInterval);

                if(idxProcToProceed != myRank && (procToProceedHasLeftInMyInterval || procToProceedHasRightInMyInterval || procIncludeMyInterval) ){
                    const int firstLeafToSend = FMath::Max(int(procToProceedLeftIdx - currentLeafsOnMyLeft), 0);
                    const int lastLeafToSend  = int(FMath::Min(procToProceedRightIdx - currentLeafsOnMyLeft, myNumberOfLeaves ));

                    //// TODO REMOVE WHEN DEBUG printf("Proc :: %d (from leaf %d to %d)\n", myRank, firstLeafToSend, lastLeafToSend);

                    while(counterLeafToSend != firstLeafToSend){
                        const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&leavesArray[offsetLeafToSend]));
                        offsetLeafToSend  += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
                        counterLeafToSend += 1;
                    }
                    const size_t offetSetToSend = offsetLeafToSend;
                    while(counterLeafToSend != lastLeafToSend){
                        const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&leavesArray[offsetLeafToSend]));
                        offsetLeafToSend  += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
                        counterLeafToSend += 1;
                    }

                    //// TODO REMOVE WHEN DEBUG printf("Proc :: %d send %d bytes to %d (from leaf %d to %d)\n",
                    //// TODO REMOVE WHEN DEBUG        myRank, int(offsetLeafToSend - offetSetToSend), idxProcToProceed, firstLeafToSend, lastLeafToSend);
                    MPI_Isend(&leavesArray[offetSetToSend], int(offsetLeafToSend - offetSetToSend), MPI_BYTE,
                              idxProcToProceed, firstLeafToSend + int(currentLeafsOnMyLeft), communicator.getComm(), &requests[counterRequest++]);
                }
                idxProcToProceed += 1;
            }
        }

        struct RecvBlockInfo{
            char* buffer;
            int nbLeaves;
        };
        RecvBlockInfo* recvBlockInfo = new RecvBlockInfo[nbProcs];
        int nbBlocksToRecv = 0;

        if(correctLeftLeavesNumber < currentLeafsOnMyLeft || currentRightLeafIdx < correctRightLeavesIndex){
            FSize iterCorrectLeafIdx = correctLeftLeavesNumber;
            int idxProcToProceed = 0;
            while(iterCorrectLeafIdx < correctRightLeavesIndex){
                if(currentLeafsOnMyLeft <= iterCorrectLeafIdx && iterCorrectLeafIdx < currentRightLeafIdx){
                    //// TODO REMOVE WHEN DEBUG printf("%d] currentLeafsOnMyLeft %llu iterCorrectLeafIdx %llu iterCorrectLeafIdx %llu currentRightLeafIdx %llu\n",
                    //// TODO REMOVE WHEN DEBUG        myRank, currentLeafsOnMyLeft, iterCorrectLeafIdx, iterCorrectLeafIdx, currentRightLeafIdx);
                    iterCorrectLeafIdx = currentRightLeafIdx;
                    idxProcToProceed   = myRank + 1;
                }
                else{
                    //// TODO REMOVE WHEN DEBUG printf("%d] currentLeafsOnMyLeft %llu iterCorrectLeafIdx %llu iterCorrectLeafIdx %llu currentRightLeafIdx %llu correctRightLeavesIndex %llu\n",
                    //// TODO REMOVE WHEN DEBUG        myRank, currentLeafsOnMyLeft, iterCorrectLeafIdx, iterCorrectLeafIdx, currentRightLeafIdx, correctRightLeavesIndex);
                    while(leavesOffsetPerProc[idxProcToProceed+1] <= iterCorrectLeafIdx){
                    //// TODO REMOVE WHEN DEBUG     printf("%d] leavesOffsetPerProc[%d+1] %llu iterCorrectLeafIdx %lld\n",
                    //// TODO REMOVE WHEN DEBUG            myRank, idxProcToProceed, leavesOffsetPerProc[idxProcToProceed+1], iterCorrectLeafIdx);
                        idxProcToProceed += 1;
                    }
                    const int nbLeafToReceive  = FMath::Min(leavesOffsetPerProc[idxProcToProceed+1], int(correctRightLeavesIndex)) - int(iterCorrectLeafIdx);
                    FSize nbParticlesToReceive = 0;
                    for(int idxLeaf = 0 ; idxLeaf < nbLeafToReceive ; ++idxLeaf){
                        nbParticlesToReceive += numberOfParticlesPerLeaf[idxLeaf + iterCorrectLeafIdx];
                    }

                    FSize bytesToRecv     = (sizeof(ParticleClass)*nbParticlesToReceive) + sizeof(int)*nbLeafToReceive;
                    char* bufferToReceive = new char[bytesToRecv];

                    //// TODO REMOVE WHEN DEBUG printf("Proc :: %d recv %d bytes to %d (from leaf %d to %d)\n",
                    //// TODO REMOVE WHEN DEBUG        myRank, bytesToRecv, idxProcToProceed, iterCorrectLeafIdx, iterCorrectLeafIdx + nbLeafToReceive);

                    MPI_Irecv(bufferToReceive, int(bytesToRecv), MPI_BYTE, idxProcToProceed, int(iterCorrectLeafIdx),
                              communicator.getComm(), &requests[counterRequest++]);

                    recvBlockInfo[nbBlocksToRecv].buffer   = bufferToReceive;
                    recvBlockInfo[nbBlocksToRecv].nbLeaves = nbLeafToReceive;
                    nbBlocksToRecv += 1;

                    iterCorrectLeafIdx += nbLeafToReceive;
                }
            }
        }

        //// TODO REMOVE WHEN DEBUG printf("%d Wait!\n", myRank);
        MPI_Waitall(counterRequest, requests, MPI_STATUSES_IGNORE);
        //// TODO REMOVE WHEN DEBUG printf("%d Done!\n", myRank);

        int idxBlockRecvInLeft = 0;
        if(correctLeftLeavesNumber < currentLeafsOnMyLeft){
            const int nbLeavesRecv = int(FMath::Min(currentLeafsOnMyLeft, correctRightLeavesIndex) - correctLeftLeavesNumber);
            //// TODO REMOVE WHEN DEBUG printf("%d] has receive %d from left\n", myRank, nbLeavesRecv);
            int idxLeaf  = 0;
            while(idxLeaf < nbLeavesRecv){
                //// TODO REMOVE WHEN DEBUG printf("%d] block %d has %d leaves\n", myRank, idxBlockRecvInLeft, recvBlockInfo[idxBlockRecvInLeft].nbLeaves);
                size_t offsetBuffer = 0;
                for(int idxLeafInBlock = 0 ; idxLeafInBlock < recvBlockInfo[idxBlockRecvInLeft].nbLeaves ; ++idxLeafInBlock){
                    const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&recvBlockInfo[idxBlockRecvInLeft].buffer[offsetBuffer]));
                    const ParticleClass*const particles   = reinterpret_cast<ParticleClass*>(&recvBlockInfo[idxBlockRecvInLeft].buffer[offsetBuffer] + sizeof(int));
                    //// TODO REMOVE WHEN DEBUG printf("%d] block %d leaf %d has %d part\n", myRank, idxBlockRecvInLeft, idxLeafInBlock, numberOfParticlesInThisLeaf);
                    for(int idxParticle = 0 ; idxParticle < numberOfParticlesInThisLeaf ; ++idxParticle){
                        particlesSaver->push(particles[idxParticle]);
                    }
                    offsetBuffer  += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
                }
                idxLeaf += recvBlockInfo[idxBlockRecvInLeft].nbLeaves;
                delete[] recvBlockInfo[idxBlockRecvInLeft].buffer;
                idxBlockRecvInLeft += 1;
            }
        }
        //// TODO REMOVE WHEN DEBUG printf("currentLeafsOnMyLeft %lld correctLeftLeavesNumber %lld currentRightLeafIdx %lld correctRightLeavesIndex %lld \n",
        //// TODO REMOVE WHEN DEBUG        currentLeafsOnMyLeft, correctLeftLeavesNumber, currentRightLeafIdx, correctRightLeavesIndex);
        if((currentLeafsOnMyLeft <= correctLeftLeavesNumber && correctLeftLeavesNumber < currentRightLeafIdx)
                || (currentLeafsOnMyLeft < correctRightLeavesIndex && correctRightLeavesIndex <= currentRightLeafIdx)){
	  const int nbLeavesToSkip = int(correctLeftLeavesNumber-currentLeafsOnMyLeft);
            size_t offsetBuffer = 0;
            //// TODO REMOVE WHEN DEBUG printf("%d] skip %d leaves\n", myRank, nbLeavesToSkip);
            for(int idxToSkip = 0 ; idxToSkip < nbLeavesToSkip ; ++idxToSkip){
                const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&leavesArray[offsetBuffer]));
                offsetBuffer += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
                //// TODO REMOVE WHEN DEBUG printf("%d] leaf %d had %d part\n", myRank, idxToSkip, numberOfParticlesInThisLeaf);
            }
            const int nbLeafToCopy = int(FMath::Min(currentRightLeafIdx, correctRightLeavesIndex) - FMath::Max(currentLeafsOnMyLeft, correctLeftLeavesNumber));
            //// TODO REMOVE WHEN DEBUG printf("%d] Need to copy %d leaves\n", myRank, nbLeafToCopy);
            for(int idxToProcess = 0 ; idxToProcess < nbLeafToCopy ; ++idxToProcess){
                const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&leavesArray[offsetBuffer]));
                //// TODO REMOVE WHEN DEBUG printf("%d] leaf %d had %d part\n", myRank, idxToProcess, numberOfParticlesInThisLeaf);
                const ParticleClass*const particles   = reinterpret_cast<ParticleClass*>(&leavesArray[offsetBuffer] + sizeof(int));
                for(int idxParticle = 0 ; idxParticle < numberOfParticlesInThisLeaf ; ++idxParticle){
                    particlesSaver->push(particles[idxParticle]);
                }

                offsetBuffer += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
            }
        }
        if(currentRightLeafIdx < correctRightLeavesIndex){
            const int nbLeavesRecv = int(correctRightLeavesIndex - FMath::Max(currentRightLeafIdx, correctLeftLeavesNumber));
            //// TODO REMOVE WHEN DEBUG printf("%d] has receive %d from right\n", myRank, nbLeavesRecv);
            int idxLeaf  = 0;
            while(idxLeaf < nbLeavesRecv){
                //// TODO REMOVE WHEN DEBUG printf("%d] block %d has %d leaves\n", myRank, idxBlockRecvInLeft, recvBlockInfo[idxBlockRecvInLeft].nbLeaves);
                size_t offsetBuffer = 0;
                for(int idxLeafInBlock = 0 ; idxLeafInBlock < recvBlockInfo[idxBlockRecvInLeft].nbLeaves ; ++idxLeafInBlock){
                    const int numberOfParticlesInThisLeaf = (*reinterpret_cast<int*>(&recvBlockInfo[idxBlockRecvInLeft].buffer[offsetBuffer]));
                    const ParticleClass*const particles   = reinterpret_cast<ParticleClass*>(&recvBlockInfo[idxBlockRecvInLeft].buffer[offsetBuffer] + sizeof(int));
                    //// TODO REMOVE WHEN DEBUG printf("%d] block %d leaf %d has %d part\n", myRank, idxBlockRecvInLeft, idxLeafInBlock, numberOfParticlesInThisLeaf);
                    for(int idxParticle = 0 ; idxParticle < numberOfParticlesInThisLeaf ; ++idxParticle){
                        particlesSaver->push(particles[idxParticle]);
                    }
                    offsetBuffer  += (sizeof(ParticleClass)*numberOfParticlesInThisLeaf+sizeof(int));
                }
                idxLeaf += recvBlockInfo[idxBlockRecvInLeft].nbLeaves;
                delete[] recvBlockInfo[idxBlockRecvInLeft].buffer;
                idxBlockRecvInLeft += 1;
            }
        }

        delete[] leavesArray;
        delete[] numberOfLeavesPerProc;
        delete[] leavesOffsetPerProc;
        delete[] numberOfParticlesPerLeaf;
        delete[] requests;
        delete[] recvBlockInfo;
    }

public:

    //////////////////////////////////////////////////////////////////////////
    // The builder function
    //////////////////////////////////////////////////////////////////////////

    template <class ContainerClass>
    static void ArrayToTree(const FMpi::FComm& communicator, const ParticleClass array[], const FSize size,
                            const FPoint& boxCenter, const FReal boxWidth, const int treeHeight,
                            ContainerClass* particleSaver, FAbstractBalanceAlgorithm* balancer,const SortingType type = QuickSort){

        IndexedParticle* particlesArray = 0;
        FSize particlesSize = 0;
        SortParticlesFromArray(communicator, array, size, type, boxCenter, boxWidth, treeHeight,
                               &particlesArray, &particlesSize);

        char* leavesArray = 0;
        FSize leavesSize = 0;
        MergeLeaves(communicator, particlesArray, particlesSize, &leavesArray, &leavesSize);

        EqualizeAndFillTree(communicator, particleSaver, leavesArray, leavesSize, balancer);
    }


};

#endif // FMPITREEBUILDER_H
