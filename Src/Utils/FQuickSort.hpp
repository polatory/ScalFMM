#ifndef FQUICKSORT_HPP
#define FQUICKSORT_HPP

#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>

#include <mpi.h>

#include "FGlobal.hpp"
#include "FMemUtils.hpp"
#include "FTrace.hpp"
#include "FMpi.hpp"

#include "FOmpBarrier.hpp"

template <class SortType, class CompareType, class IndexType>
class FQuickSort {
    ////////////////////////////////////////////////////////////
    // Miscialenous functions
    ////////////////////////////////////////////////////////////


    /** To get max between 2 values */
    template <class NumType>
    static NumType Max(const NumType inV1, const NumType inV2){
        return (inV1 > inV2 ? inV1 : inV2);
    }

    /** To get min between 2 values */
    template <class NumType>
    static NumType Min(const NumType inV1, const NumType inV2){
        return (inV1 < inV2 ? inV1 : inV2);
    }

    /** swap to value */
    template <class NumType>
    static inline void Swap(NumType& value, NumType& other){
        NumType temp = value;
        value = other;
        other = temp;
    }


    ////////////////////////////////////////////////////////////
    // Split information
    ////////////////////////////////////////////////////////////

    static IndexType getLeft(const IndexType inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        return IndexType(ceil(step * inIdProc));
    }

    static IndexType getRight(const IndexType inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        const IndexType res = IndexType(ceil(step * (inIdProc+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    static int getProc(const IndexType position, const IndexType inSize, const int inNbProc) {
        const double step = double(inSize) / double(inNbProc);
        return int(double(position)/step);
    }


    ////////////////////////////////////////////////////////////
    // MPI Function
    ////////////////////////////////////////////////////////////

    /** generic mpi assert function */
    static void mpiassert(const int test, const unsigned line, const char* const message = 0){
        if(test != MPI_SUCCESS){
            printf("[ERROR-QS] Test failled at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

    /** get current rank */
    static int MpiGetRank(MPI_Comm comm){
        int rank(0);
        mpiassert( MPI_Comm_rank(comm, &rank),  __LINE__ );
        return rank;
    }

    /** get current nb procs */
    static int MpiGetNbProcs(MPI_Comm comm){
        int nb(0);
        mpiassert( MPI_Comm_size(comm, &nb), __LINE__);
        return nb;
    }

    /** get the pivot from several procs in Comm */
    static CompareType MpiGetPivot(CompareType myData, MPI_Comm& comm, const int nbProcs){
        CompareType result[nbProcs];
        mpiassert( MPI_Allgather( &myData, sizeof(CompareType), MPI_BYTE, result, sizeof(CompareType), MPI_BYTE, comm),  __LINE__ );
        // We do an average of the first element of each proc array
        CompareType sum = 0;
        for(int idxProc = 0 ; idxProc < nbProcs ;++idxProc){
            sum += result[idxProc] / nbProcs;
        }
        return sum ;
    }

    /** change the group and communicator */
    static void MpiChangeGroup(MPI_Group& currentGroup, MPI_Comm& currentComm, const int from , const int to){
        int procsIdArray[to - from + 1];
        for(int idxProc = from ;idxProc <= to ; ++idxProc){
            procsIdArray[idxProc - from] = idxProc;
        }

        MPI_Group previousGroup = currentGroup;
        mpiassert( MPI_Group_incl(previousGroup, to - from + 1 , procsIdArray, &currentGroup),  __LINE__ );

        MPI_Comm previousComm = currentComm;
        mpiassert( MPI_Comm_create(previousComm, currentGroup, &currentComm),  __LINE__ );

        MPI_Comm_free(&previousComm);
        MPI_Group_free(&previousGroup);
    }

    ////////////////////////////////////////////////////////////
    // Quick sort
    ////////////////////////////////////////////////////////////

    /* use in the sequential qs */
    static IndexType QsPartition(SortType array[], IndexType left, IndexType right){
        const IndexType part = right;
        Swap(array[part],array[(right + left ) / 2]);
        --right;

        while(true){
            while(CompareType(array[left]) < CompareType(array[part])){
                ++left;
            }
            while(right >= left && CompareType(array[part]) <= CompareType(array[right])){
                --right;
            }
            if(right < left) break;

            Swap(array[left],array[right]);
            ++left;
            --right;
        }

        Swap(array[part],array[left]);

        return left;
    }

    /* a local iteration of qs */
    static void QsLocal(SortType array[], const CompareType& pivot,
               IndexType myLeft, IndexType myRight,
               IndexType& prefix, IndexType& sufix){

        IndexType leftIter = myLeft;
        IndexType rightIter = myRight;

        while(true){
            while(CompareType(array[leftIter]) <= pivot && leftIter < rightIter){
                ++leftIter;
            }
            while(leftIter <= rightIter && pivot < CompareType(array[rightIter])){
                --rightIter;
            }
            if(rightIter < leftIter) break;

            Swap(array[leftIter],array[rightIter]);
            ++leftIter;
            --rightIter;
        }

        prefix = leftIter - myLeft;
        sufix = myRight - myLeft - prefix + 1;
    }

    /* a sequential qs */
    static void QsSequentialStep(SortType array[], const IndexType left, const IndexType right){
        if(left < right){
            const IndexType part = QsPartition(array, left, right);
            QsSequentialStep(array,part + 1,right);
            QsSequentialStep(array,left,part - 1);
        }
    }

    /** A task dispatcher */
    static void QsOmpTask(SortType array[], const IndexType left, const IndexType right, const int deep){
        if(left < right){
            if( deep ){
                const IndexType part = QsPartition(array, left, right);
                #pragma omp task
                QsOmpTask(array,part + 1,right, deep - 1);
                #pragma omp task
                QsOmpTask(array,left,part - 1, deep - 1);
            }
            else {
                const IndexType part = QsPartition(array, left, right);
                QsSequentialStep(array,part + 1,right);
                QsSequentialStep(array,left,part - 1);
            }
        }
    }

    /* the openmp qs */
    static void QsOmpNoTask(SortType array[], const IndexType size){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Quicksort" , __FILE__ , __LINE__) );
        struct Fix{
            IndexType pre;
            IndexType suf;
        };

        const int NbOfThreads = omp_get_max_threads();

        Fix fixes[NbOfThreads + 1];
        Fix allFixesSum[NbOfThreads + 1][NbOfThreads];

        memset(fixes,0,sizeof(Fix) * NbOfThreads);
        memset(allFixesSum,0,sizeof(Fix) * (NbOfThreads + 1) * NbOfThreads);

        SortType*const temporaryArray = reinterpret_cast<SortType*>(new char[sizeof(SortType) * size]);

        FOmpBarrier barriers[NbOfThreads];

        #pragma omp parallel
        {
            const int myThreadId = omp_get_thread_num();

            IndexType myLeft = getLeft(size, myThreadId, omp_get_num_threads());
            IndexType myRight = getRight(size, myThreadId, omp_get_num_threads()) - 1;

            IndexType startIndex = 0;
            IndexType endIndex = size - 1;

            int firstProc = 0;
            int lastProc = omp_get_num_threads() - 1;

            while( firstProc != lastProc && (endIndex - startIndex + 1) != 0){
                Fix* const fixesSum = &allFixesSum[0][firstProc];
                const IndexType nbElements = endIndex - startIndex + 1;

                if(myThreadId == firstProc){
                    barriers[firstProc].setNbThreads( lastProc - firstProc + 1);
                }

                // sort QsLocal part of the array
                const CompareType pivot = (CompareType(array[startIndex]) + CompareType(array[endIndex]) )/2;
                barriers[firstProc].wait();

                QsLocal(array, pivot, myLeft, myRight, fixes[myThreadId].pre, fixes[myThreadId].suf);

                // wait others that work on this part
                #pragma omp flush(array)
                barriers[firstProc].wait();

                // merge result
                if(myThreadId == firstProc){
                    fixesSum[firstProc].pre = 0;
                    fixesSum[firstProc].suf = 0;
                    for(int idxProc = firstProc ; idxProc <= lastProc ; ++idxProc){
                        fixesSum[idxProc + 1].pre = fixesSum[idxProc].pre + fixes[idxProc].pre;
                        fixesSum[idxProc + 1].suf = fixesSum[idxProc].suf + fixes[idxProc].suf;
                    }
                    #pragma omp flush(fixesSum)
                }
                // prepare copy
                if(myThreadId == firstProc + 1){
                    FMemUtils::memcpy(&temporaryArray[startIndex], &array[startIndex], sizeof(SortType) * nbElements );
                    #pragma omp flush(temporaryArray)
                }

                barriers[firstProc].wait();

                // copy my result where it belong (< pivot)
                FMemUtils::memcpy(&array[startIndex + fixesSum[myThreadId].pre], &temporaryArray[myLeft], sizeof(SortType) * fixes[myThreadId].pre);

                // copy my result where it belong (> pivot)
                const IndexType sufoffset = fixesSum[lastProc + 1].pre + startIndex;
                FMemUtils::memcpy(&array[sufoffset + fixesSum[myThreadId].suf], &temporaryArray[myLeft + fixes[myThreadId].pre ], sizeof(SortType) * fixes[myThreadId].suf);

                barriers[firstProc].wait();

                // find my next QsLocal part
                int splitProc = getProc(sufoffset - startIndex, nbElements, lastProc - firstProc + 1) + firstProc;
                if(splitProc == lastProc){
                    --splitProc;
                }

                if( myThreadId <= splitProc ){
                    endIndex = sufoffset - 1;
                    lastProc = splitProc;
                }
                else{
                    startIndex = sufoffset;
                    firstProc = splitProc + 1;
                }

                myLeft = getLeft(endIndex - startIndex + 1, myThreadId - firstProc, lastProc - firstProc + 1) + startIndex;
                myRight = getRight(endIndex - startIndex + 1, myThreadId - firstProc, lastProc - firstProc + 1) + startIndex - 1;
            }

            QsSequentialStep(array,myLeft,myRight);
        }

        delete[] reinterpret_cast<char*>(temporaryArray);
    }

public:
    /* a sequential qs */
    static void QsSequential(SortType array[], const IndexType size){
        QsSequentialStep(array,0,size-1);
    }

    /** The openmp quick sort */
    static void QsOmp(SortType array[], const IndexType size){
        #if _OPENMP >= 200805
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                QsOmpTask(array, 0, size - 1 , 15);
            }
        }
        #else
        QsOmpNoTask(array, size);
        #endif
    }

    /* the mpi qs */
    static void QsMpi(const SortType originalArray[], IndexType size, SortType* & outputArray, IndexType& outputSize, MPI_Comm originalComm = MPI_COMM_WORLD){
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
        Fix fixes[MpiGetNbProcs(MPI_COMM_WORLD) + 1];
        Fix fixesSum[MpiGetNbProcs(MPI_COMM_WORLD) + 1];
        memset(fixes,0,sizeof(Fix) * MpiGetNbProcs(MPI_COMM_WORLD));
        memset(fixesSum,0,sizeof(Fix) * (MpiGetNbProcs(MPI_COMM_WORLD) + 1) );

        // receiving buffer
        IndexType bufferSize = 0;
        SortType* buffer = 0;

        // Create the first group
        MPI_Group currentGroup;
        // Create the first com
        MPI_Comm currentComm;

        mpiassert( MPI_Comm_dup(originalComm, &currentComm),  __LINE__ );
        mpiassert( MPI_Comm_group(currentComm, &currentGroup),  __LINE__ );

        // While I am not working alone on my own data
        while( MpiGetNbProcs(currentComm) != 1 ){
            const int currentRank = MpiGetRank(currentComm);
            const int currentNbProcs = MpiGetNbProcs(currentComm);

            MPI_Request requests[currentNbProcs * 2];
            int iterRequest = 0;

            /////////////////////////////////////////////////
            // Local sort
            /////////////////////////////////////////////////

            // sort QsLocal part of the outputArray
            const CompareType pivot = MpiGetPivot(outputArray[size/2], currentComm, currentNbProcs);
            Fix myFix;
            QsLocal(outputArray, pivot, 0, size - 1, myFix.pre, myFix.suf);

            // exchange fixes
            mpiassert( MPI_Allgather( &myFix, sizeof(Fix), MPI_BYTE, fixes, sizeof(Fix), MPI_BYTE, currentComm),  __LINE__ );

            // each procs compute the summation
            fixesSum[0].pre = 0;
            fixesSum[0].suf = 0;
            for(int idxProc = 0 ; idxProc < currentNbProcs ; ++idxProc){
                fixesSum[idxProc + 1].pre = fixesSum[idxProc].pre + fixes[idxProc].pre;
                fixesSum[idxProc + 1].suf = fixesSum[idxProc].suf + fixes[idxProc].suf;
            }

            // then I need to know which procs will be in the middle
            int splitProc = getProc(fixesSum[currentNbProcs].pre - 1, fixesSum[currentNbProcs].pre + fixesSum[currentNbProcs].suf, currentNbProcs);
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

                const int firstProcToSend = getProc(fixesSum[currentRank].suf, elementsInSuf, procsInSuf) + firstProcInSuf;
                const int lastProcToSend = getProc(fixesSum[currentRank + 1].suf - 1, elementsInSuf, procsInSuf) + firstProcInSuf;

                IndexType sent = 0;
                for(int idxProc = firstProcToSend ; idxProc <= lastProcToSend ; ++idxProc){
                    const IndexType thisProcRight = getRight(elementsInSuf, idxProc - firstProcInSuf, procsInSuf);
                    IndexType sendToProc = thisProcRight - fixesSum[currentRank].suf - sent;

                    if(sendToProc + sent > fixes[currentRank].suf){
                        sendToProc = fixes[currentRank].suf - sent;
                    }
                    if( sendToProc ){
                        mpiassert( MPI_Isend(&outputArray[sent + fixes[currentRank].pre], int(sendToProc * sizeof(SortType)), MPI_BYTE , idxProc, FMpi::TagQuickSort, currentComm, &requests[iterRequest++]),  __LINE__ );
                    }
                    sent += sendToProc;
                }
            }

            // under pivot (left part)
            if( fixes[currentRank].pre ){
                const int procsInPre = splitProc + 1;
                const IndexType elementsInPre = fixesSum[currentNbProcs].pre;

                const int firstProcToSend = getProc(fixesSum[currentRank].pre, elementsInPre, procsInPre);
                const int lastProcToSend = getProc(fixesSum[currentRank + 1].pre - 1, elementsInPre, procsInPre);

                IndexType sent = 0;
                for(int idxProc = firstProcToSend ; idxProc <= lastProcToSend ; ++idxProc){
                    const IndexType thisProcRight = getRight(elementsInPre, idxProc, procsInPre);
                    IndexType sendToProc = thisProcRight - fixesSum[currentRank].pre - sent;

                    if(sendToProc + sent > fixes[currentRank].pre){
                        sendToProc = fixes[currentRank].pre - sent;
                    }
                    if(sendToProc){
                        mpiassert( MPI_Isend(&outputArray[sent], int(sendToProc * sizeof(SortType)), MPI_BYTE , idxProc, FMpi::TagQuickSort, currentComm, &requests[iterRequest++]),  __LINE__ );
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

                IndexType myLeft = getLeft(elementsInPre, currentRank, procsInPre);
                IndexType myRightLimit = getRight(elementsInPre, currentRank, procsInPre);

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
                    const IndexType firstIndex = Max(myLeft , fixesSum[idxProc].pre );
                    const IndexType endIndex = Min(fixesSum[idxProc + 1].pre,  myRightLimit);
                    if( (endIndex - firstIndex) ){
                        mpiassert( MPI_Irecv(&buffer[indexArray], int((endIndex - firstIndex) * sizeof(SortType)), MPI_BYTE, idxProc, FMpi::TagQuickSort, currentComm, &requests[iterRequest++]),  __LINE__ );
                    }
                    indexArray += endIndex - firstIndex;
                    ++idxProc;
                }
                // Proceed all send/receive
                mpiassert( MPI_Waitall(iterRequest, requests, MPI_STATUSES_IGNORE),  __LINE__ );

                MpiChangeGroup(currentGroup, currentComm, 0, splitProc);
            }
            else{
                // I am in L-Part (larger than pivot)
                const int procsInSuf = currentNbProcs - 1 - splitProc;
                const IndexType elementsInSuf = fixesSum[currentNbProcs].suf;

                const int rankInL = currentRank - splitProc - 1;
                IndexType myLeft = getLeft(elementsInSuf, rankInL, procsInSuf);
                IndexType myRightLimit = getRight(elementsInSuf, rankInL, procsInSuf);

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
                    const IndexType firstIndex = Max(myLeft , fixesSum[idxProc].suf );
                    const IndexType endIndex = Min(fixesSum[idxProc + 1].suf,  myRightLimit);
                    if( (endIndex - firstIndex) ){
                        mpiassert( MPI_Irecv(&buffer[indexArray], int((endIndex - firstIndex) * sizeof(SortType)), MPI_BYTE, idxProc, FMpi::TagQuickSort, currentComm, &requests[iterRequest++]),  __LINE__ );
                    }
                    indexArray += endIndex - firstIndex;
                    ++idxProc;
                }
                // Proceed all send/receive
                mpiassert( MPI_Waitall(iterRequest, requests, MPI_STATUSES_IGNORE),  __LINE__ );

                MpiChangeGroup(currentGroup, currentComm, splitProc + 1, currentNbProcs - 1);
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
        MPI_Comm_free(&currentComm);
        MPI_Group_free(&currentGroup);

        // Normal Quick sort
        QsOmp(outputArray, size);
        outputSize = size;
    }

};

#endif // FQUICKSORT_HPP
