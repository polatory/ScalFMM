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
#ifndef FQUICKSORTMPI_HPP
#define FQUICKSORTMPI_HPP

#include "FQuickSort.hpp"
#include "FMpi.hpp"


template <class SortType, class CompareType, class IndexType>
class FQuickSortMpi : public FQuickSort< SortType, CompareType, IndexType> {
public:
    /* the mpi qs */
    static void QsMpi(const SortType originalArray[], IndexType size, SortType* & outputArray, IndexType& outputSize, const FMpi::FComm& originalComm){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Quicksort" , __FILE__ , __LINE__) );
        // We need a structure see the algorithm detail to know more
        struct Fix{
            IndexType pre;
            IndexType suf;
        };

        // first we copy data into our working buffer : outputArray
        outputArray = new SortType[size];
        FMemUtils::memcpy(outputArray, originalArray, sizeof(SortType) * size);
        outputSize = size;

        // alloc outputArray to store pre/sufixe, maximum needed is nb procs[comm world] + 1
        Fix fixes[originalComm.processCount() + 1];
        Fix fixesSum[originalComm.processCount() + 1];
        memset(fixes,0,sizeof(Fix) * originalComm.processCount());
        memset(fixesSum,0,sizeof(Fix) * (originalComm.processCount() + 1) );

        // receiving buffer
        IndexType bufferSize = 0;
        SortType* buffer = 0;

        // Create the first com
        FMpi::FComm currentComm(originalComm.getComm());

        // While I am not working alone on my own data
        while( currentComm.processCount() != 1 ){
            const int currentRank = currentComm.processId();
            const int currentNbProcs = currentComm.processCount();

            MPI_Request requests[currentNbProcs * 2];
            int iterRequest = 0;

            /////////////////////////////////////////////////
            // Local sort
            /////////////////////////////////////////////////

            // sort QsLocal part of the outputArray
            const CompareType pivot = currentComm.reduceAverageAll( CompareType(outputArray[size/2]) );
            Fix myFix;
            FQuickSort< SortType, CompareType, IndexType>::QsLocal(outputArray, pivot, 0, size - 1, myFix.pre, myFix.suf);

            // exchange fixes
            FMpi::Assert( MPI_Allgather( &myFix, sizeof(Fix), MPI_BYTE, fixes, sizeof(Fix), MPI_BYTE, currentComm.getComm()),  __LINE__ );

            // each procs compute the summation
            fixesSum[0].pre = 0;
            fixesSum[0].suf = 0;
            for(int idxProc = 0 ; idxProc < currentNbProcs ; ++idxProc){
                fixesSum[idxProc + 1].pre = fixesSum[idxProc].pre + fixes[idxProc].pre;
                fixesSum[idxProc + 1].suf = fixesSum[idxProc].suf + fixes[idxProc].suf;
            }

            // then I need to know which procs will be in the middle
            int splitProc = FMpi::GetProc(fixesSum[currentNbProcs].pre - 1, fixesSum[currentNbProcs].pre + fixesSum[currentNbProcs].suf, currentNbProcs);
            if(splitProc == currentNbProcs - 1){
                --splitProc;
            }

            /////////////////////////////////////////////////
            // Send my data
            /////////////////////////////////////////////////

            // above pivot (right part)
            if( fixes[currentRank].suf ){
                const int procsInSuf = currentNbProcs - 1 - splitProc;
                const int firstProcInSuf = splitProc + 1;
                const IndexType elementsInSuf = fixesSum[currentNbProcs].suf;

                const int firstProcToSend = FMpi::GetProc(fixesSum[currentRank].suf, elementsInSuf, procsInSuf) + firstProcInSuf;
                const int lastProcToSend = FMpi::GetProc(fixesSum[currentRank + 1].suf - 1, elementsInSuf, procsInSuf) + firstProcInSuf;

                IndexType sent = 0;
                for(int idxProc = firstProcToSend ; idxProc <= lastProcToSend ; ++idxProc){
                    const IndexType thisProcRight = FMpi::GetRight(elementsInSuf, idxProc - firstProcInSuf, procsInSuf);
                    IndexType sendToProc = thisProcRight - fixesSum[currentRank].suf - sent;

                    if(sendToProc + sent > fixes[currentRank].suf){
                        sendToProc = fixes[currentRank].suf - sent;
                    }
                    if( sendToProc ){
                        FMpi::Assert( MPI_Isend(&outputArray[sent + fixes[currentRank].pre], int(sendToProc * sizeof(SortType)), MPI_BYTE , idxProc, FMpi::TagQuickSort, currentComm.getComm(), &requests[iterRequest++]),  __LINE__ );
                    }
                    sent += sendToProc;
                }
            }

            // under pivot (left part)
            if( fixes[currentRank].pre ){
                const int procsInPre = splitProc + 1;
                const IndexType elementsInPre = fixesSum[currentNbProcs].pre;

                const int firstProcToSend = FMpi::GetProc(fixesSum[currentRank].pre, elementsInPre, procsInPre);
                const int lastProcToSend = FMpi::GetProc(fixesSum[currentRank + 1].pre - 1, elementsInPre, procsInPre);

                IndexType sent = 0;
                for(int idxProc = firstProcToSend ; idxProc <= lastProcToSend ; ++idxProc){
                    const IndexType thisProcRight = FMpi::GetRight(elementsInPre, idxProc, procsInPre);
                    IndexType sendToProc = thisProcRight - fixesSum[currentRank].pre - sent;

                    if(sendToProc + sent > fixes[currentRank].pre){
                        sendToProc = fixes[currentRank].pre - sent;
                    }
                    if(sendToProc){
                        FMpi::Assert( MPI_Isend(&outputArray[sent], int(sendToProc * sizeof(SortType)), MPI_BYTE , idxProc, FMpi::TagQuickSort, currentComm.getComm(), &requests[iterRequest++]),  __LINE__ );
                    }
                    sent += sendToProc;
                }
            }

            /////////////////////////////////////////////////
            // Receive data that belong to me
            /////////////////////////////////////////////////

            if( currentRank <= splitProc ){
                // I am in S-Part (smaller than pivot)
                const int procsInPre = splitProc + 1;
                const IndexType elementsInPre = fixesSum[currentNbProcs].pre;

                IndexType myLeft = FMpi::GetLeft(elementsInPre, currentRank, procsInPre);
                IndexType myRightLimit = FMpi::GetRight(elementsInPre, currentRank, procsInPre);

                size = myRightLimit - myLeft;
                if(bufferSize < size){
                    bufferSize = size;
                    delete[] buffer;
                    buffer = new SortType[bufferSize];
                }

                int idxProc = 0;
                while( idxProc < currentNbProcs && fixesSum[idxProc + 1].pre <= myLeft ){
                    ++idxProc;
                }

                IndexType indexArray = 0;

                while( idxProc < currentNbProcs && indexArray < myRightLimit - myLeft){
                    const IndexType firstIndex = FMath::Max(myLeft , fixesSum[idxProc].pre );
                    const IndexType endIndex = FMath::Min(fixesSum[idxProc + 1].pre,  myRightLimit);
                    if( (endIndex - firstIndex) ){
                        FMpi::Assert( MPI_Irecv(&buffer[indexArray], int((endIndex - firstIndex) * sizeof(SortType)), MPI_BYTE, idxProc, FMpi::TagQuickSort, currentComm.getComm(), &requests[iterRequest++]),  __LINE__ );
                    }
                    indexArray += endIndex - firstIndex;
                    ++idxProc;
                }
                // Proceed all send/receive
                FMpi::Assert( MPI_Waitall(iterRequest, requests, MPI_STATUSES_IGNORE),  __LINE__ );

                currentComm.groupReduce( 0, splitProc);
            }
            else{
                // I am in L-Part (larger than pivot)
                const int procsInSuf = currentNbProcs - 1 - splitProc;
                const IndexType elementsInSuf = fixesSum[currentNbProcs].suf;

                const int rankInL = currentRank - splitProc - 1;
                IndexType myLeft = FMpi::GetLeft(elementsInSuf, rankInL, procsInSuf);
                IndexType myRightLimit = FMpi::GetRight(elementsInSuf, rankInL, procsInSuf);

                size = myRightLimit - myLeft;
                if(bufferSize < size){
                    bufferSize = size;
                    delete[] buffer;
                    buffer = new SortType[bufferSize];
                }

                int idxProc = 0;
                while( idxProc < currentNbProcs && fixesSum[idxProc + 1].suf <= myLeft ){
                    ++idxProc;
                }

                IndexType indexArray = 0;

                while( idxProc < currentNbProcs && indexArray < myRightLimit - myLeft){
                    const IndexType firstIndex = FMath::Max(myLeft , fixesSum[idxProc].suf );
                    const IndexType endIndex = FMath::Min(fixesSum[idxProc + 1].suf,  myRightLimit);
                    if( (endIndex - firstIndex) ){
                        FMpi::Assert( MPI_Irecv(&buffer[indexArray], int((endIndex - firstIndex) * sizeof(SortType)), MPI_BYTE, idxProc, FMpi::TagQuickSort, currentComm.getComm(), &requests[iterRequest++]),  __LINE__ );
                    }
                    indexArray += endIndex - firstIndex;
                    ++idxProc;
                }
                // Proceed all send/receive
                FMpi::Assert( MPI_Waitall(iterRequest, requests, MPI_STATUSES_IGNORE),  __LINE__ );

                currentComm.groupReduce( splitProc + 1, currentNbProcs - 1);
            }



            // Copy res into outputArray
            if(outputSize < size){
                delete[] outputArray;
                outputArray = new SortType[size];
                outputSize = size;
            }

            FMemUtils::memcpy(outputArray, buffer, sizeof(SortType) * size);
        }

        /////////////////////////////////////////////////
        // End QsMpi sort
        /////////////////////////////////////////////////

        // Clean
        delete[] buffer;

        // Normal Quick sort
        FQuickSort< SortType, CompareType, IndexType>::QsOmp(outputArray, size);
        outputSize = size;
    }

};

#endif // FQUICKSORTMPI_HPP
