#ifndef FQUICKSORT_HPP
#define FQUICKSORT_HPP

#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>

#include <mpi.h>


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

    static long getLeft(const long inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        return long(ceil(step * inIdProc));
    }

    static long getRight(const long inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        const long res = long(ceil(step * (inIdProc+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    static long getOtherRight(const long inSize, const int other, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        const long res = long(ceil(step * (other+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    static int getProc(const int position, const long inSize, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        return int(position/step);
    }

    ////////////////////////////////////////////////////////////
    // OMP Function
    ////////////////////////////////////////////////////////////

    /* custom barrier to wait proc from first to last, not all threads! */
    static void OmpBarrier(int mutex[], const int firstProc, const int lastProc, const int myThreadId){
        if(lastProc != firstProc){
            const int idRelative = myThreadId - firstProc;

            while(mutex[firstProc] != idRelative ){
                #pragma omp flush(mutex)
            }

            ++mutex[firstProc];
            #pragma omp flush(mutex)

            if(myThreadId == lastProc){
                mutex[firstProc] = idRelative - 1;
            }
            else{
                while(mutex[firstProc] != idRelative ){
                    #pragma omp flush(mutex)
                }
                if(idRelative != 0){
                    --mutex[firstProc];
                    #pragma omp flush(mutex)
                }
            }
        }
    }


    ////////////////////////////////////////////////////////////
    // MPI Function
    ////////////////////////////////////////////////////////////

    /** generic mpi assert function */
    static void mpiassert(const int test, const unsigned line, const char* const message = 0){
        if(test != MPI_SUCCESS){
            printf("[ERROR] Test failled at line %d, result is %d", line, test);
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
    template <class PivotType>
    static PivotType MpiGetPivot(PivotType myData, MPI_Comm& comm, const int nbProcs){
        PivotType result[nbProcs];
        mpiassert( MPI_Allgather( &myData, sizeof(PivotType), MPI_BYTE, result, sizeof(PivotType), MPI_BYTE, comm),  __LINE__ );
        // We do an average of the first element of each proc array
        PivotType sum = 0;
        for(int idxProc = 0 ; idxProc < nbProcs ;++idxProc){
            sum += result[idxProc];
        }
        return sum / nbProcs;
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
    template <class SortType, class PivotType>
    static long QsPartition(SortType array[], long left, long right){
        const long part = right;
        Swap(array[part],array[(right + left ) / 2]);
        --right;

        while(true){
            while(PivotType(array[left]) < PivotType(array[part])){
                ++left;
            }
            while(right >= left && PivotType(array[part]) <= PivotType(array[right])){
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
    template <class SortType, class PivotType>
    static void QsLocal(SortType array[], const PivotType& pivot,
               long myLeft, long myRight,
               int& prefix, int& sufix){

        long leftIter = myLeft;
        long rightIter = myRight;

        while(true){
            while(PivotType(array[leftIter]) <= pivot && leftIter < rightIter){
                ++leftIter;
            }
            while(leftIter <= rightIter && pivot < PivotType(array[rightIter])){
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

public:
    /* a sequential qs */
    template <class SortType, class PivotType>
    static void QsSequential(SortType array[], const long left, const long right){
        if(left < right){
            const long part = QsPartition<SortType,PivotType>(array, left, right);
            QsSequential<SortType,PivotType>(array,part + 1,right);
            QsSequential<SortType,PivotType>(array,left,part - 1);
        }
    }

    /* the openmp qs */
    template <class SortType, class PivotType>
    static void QsOmp(SortType array[], const long size){
        struct Fix{
            int pre;
            int suf;
        };

        const int NbOfThreads = omp_get_max_threads();

        Fix fixes[NbOfThreads + 1];
        Fix allFixesSum[NbOfThreads + 1][NbOfThreads];

        memset(fixes,0,sizeof(Fix) * NbOfThreads);
        memset(allFixesSum,0,sizeof(Fix) * (NbOfThreads + 1) * NbOfThreads);

        SortType*const temporaryArray = reinterpret_cast<SortType*>(new char[sizeof(SortType) * size]);

        int mutex[NbOfThreads];
        memset(mutex, 0, sizeof(int) * NbOfThreads);

        #pragma omp parallel
        {
            const int myThreadId = omp_get_thread_num();

            long myLeft = getLeft(size, myThreadId, omp_get_num_threads());
            long myRight = getRight(size, myThreadId, omp_get_num_threads()) - 1;

            long startIndex = 0;
            long endIndex = size - 1;

            int firstProc = 0;
            int lastProc = omp_get_num_threads() - 1;

            while( firstProc != lastProc && (endIndex - startIndex + 1) != 0){
                Fix* const fixesSum = &allFixesSum[0][firstProc];
                const long nbElements = endIndex - startIndex + 1;

                // sort QsLocal part of the array
                const PivotType pivot = (PivotType(array[startIndex]) + PivotType(array[endIndex]) )/2;
                OmpBarrier( mutex, firstProc, lastProc, myThreadId);

                QsLocal(array, pivot, myLeft, myRight, fixes[myThreadId].pre, fixes[myThreadId].suf);

                // wait others that work on this part
                #pragma omp flush(array)
                OmpBarrier( mutex, firstProc, lastProc, myThreadId);

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
                    memcpy(&temporaryArray[startIndex], &array[startIndex], sizeof(SortType) * nbElements );
                    #pragma omp flush(temporaryArray)
                }

                OmpBarrier( mutex, firstProc, lastProc, myThreadId);

                // copy my result where it belong (< pivot)
                memcpy(&array[startIndex + fixesSum[myThreadId].pre], &temporaryArray[myLeft], sizeof(SortType) * fixes[myThreadId].pre);

                // copy my result where it belong (> pivot)
                const int sufoffset = fixesSum[lastProc + 1].pre + startIndex;
                memcpy(&array[sufoffset + fixesSum[myThreadId].suf], &temporaryArray[myLeft + fixes[myThreadId].pre ], sizeof(SortType) * fixes[myThreadId].suf);

                OmpBarrier( mutex, firstProc, lastProc, myThreadId);

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

            QsSequential<SortType,PivotType>(array,myLeft,myRight);
        }

        delete[] reinterpret_cast<char*>(temporaryArray);
    }

    /* the mpi qs */
    template <class SortType, class PivotType>
    static void QsMpi(const SortType originalArray[], long size, SortType* & outputArray, long& outputSize, MPI_Comm originalComm = MPI_COMM_WORLD){
        // We need a structure see the algorithm detail to know more
        struct Fix{
            int pre;
            int suf;
        };

        // first we copy data into our working buffer : outputArray
        outputArray = new SortType[size];
        memcpy(outputArray, originalArray, sizeof(SortType) * size);
        outputSize = size;

        // alloc outputArray to store pre/sufixe, maximum needed is nb procs[comm world] + 1
        Fix fixes[MpiGetNbProcs(MPI_COMM_WORLD) + 1];
        Fix fixesSum[MpiGetNbProcs(MPI_COMM_WORLD) + 1];
        memset(fixes,0,sizeof(Fix) * MpiGetNbProcs(MPI_COMM_WORLD));
        memset(fixesSum,0,sizeof(Fix) * (MpiGetNbProcs(MPI_COMM_WORLD) + 1) );

        // receiving buffer
        long bufferSize = 0;
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
            const PivotType pivot = MpiGetPivot<PivotType>(outputArray[0], currentComm, currentNbProcs);
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
                const int elementsInSuf = fixesSum[currentNbProcs].suf;

                const int firstProcToSend = getProc(fixesSum[currentRank].suf, elementsInSuf, procsInSuf) + firstProcInSuf;
                const int lastProcToSend = getProc(fixesSum[currentRank + 1].suf - 1, elementsInSuf, procsInSuf) + firstProcInSuf;

                int sent = 0;
                for(int idxProc = firstProcToSend ; idxProc <= lastProcToSend ; ++idxProc){
                    const int thisProcRight = getRight(elementsInSuf, idxProc - firstProcInSuf, procsInSuf);
                    int sendToProc = thisProcRight - fixesSum[currentRank].suf - sent;

                    if(sendToProc + sent > fixes[currentRank].suf){
                        sendToProc = fixes[currentRank].suf - sent;
                    }

                    mpiassert( MPI_Isend(&outputArray[sent + fixes[currentRank].pre], sendToProc * sizeof(SortType), MPI_BYTE , idxProc, 0, currentComm, &requests[iterRequest++]),  __LINE__ );
                    sent += sendToProc;
                }
            }

            // under pivot (left part)
            if( fixes[currentRank].pre ){
                const int procsInPre = splitProc + 1;
                const int elementsInPre = fixesSum[currentNbProcs].pre;

                const int firstProcToSend = getProc(fixesSum[currentRank].pre, elementsInPre, procsInPre);
                const int lastProcToSend = getProc(fixesSum[currentRank + 1].pre - 1, elementsInPre, procsInPre);

                int sent = 0;
                for(int idxProc = firstProcToSend ; idxProc <= lastProcToSend ; ++idxProc){
                    const int thisProcRight = getRight(elementsInPre, idxProc, procsInPre);
                    int sendToProc = thisProcRight - fixesSum[currentRank].pre - sent;

                    if(sendToProc + sent > fixes[currentRank].pre){
                        sendToProc = fixes[currentRank].pre - sent;
                    }

                    mpiassert( MPI_Isend(&outputArray[sent], sendToProc * sizeof(SortType), MPI_BYTE , idxProc, 0, currentComm, &requests[iterRequest++]),  __LINE__ );
                    sent += sendToProc;
                }
            }

            /////////////////////////////////////////////////
            // Receive data that belong to me
            /////////////////////////////////////////////////

            if( currentRank <= splitProc ){
                // I am in S-Part (smaller than pivot)
                const int procsInPre = splitProc + 1;
                const int elementsInPre = fixesSum[currentNbProcs].pre;

                long myLeft = getLeft(elementsInPre, currentRank, procsInPre);
                long myRightLimit = getRight(elementsInPre, currentRank, procsInPre);

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

                int indexArray = 0;

                while( idxProc < currentNbProcs && indexArray < myRightLimit - myLeft){
                    const int firstIndex = Max(myLeft , (long int) fixesSum[idxProc].pre );
                    const int endIndex = Min((long int)fixesSum[idxProc + 1].pre,  myRightLimit);

                    mpiassert( MPI_Irecv(&buffer[indexArray], (endIndex - firstIndex) * sizeof(SortType), MPI_BYTE, idxProc, 0, currentComm, &requests[iterRequest++]),  __LINE__ );
                    indexArray += endIndex - firstIndex;
                    ++idxProc;
                }

                MpiChangeGroup(currentGroup, currentComm, 0, splitProc);
            }
            else{
                // I am in L-Part (larger than pivot)
                const int procsInSuf = currentNbProcs - 1 - splitProc;
                const int elementsInSuf = fixesSum[currentNbProcs].suf;

                const int rankInL = currentRank - splitProc - 1;
                long myLeft = getLeft(elementsInSuf, rankInL, procsInSuf);
                long myRightLimit = getRight(elementsInSuf, rankInL, procsInSuf);

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

                int indexArray = 0;

                while( idxProc < currentNbProcs && indexArray < myRightLimit - myLeft){
                    const int firstIndex = Max(myLeft , (long int)fixesSum[idxProc].suf );
                    const int endIndex = Min((long int)fixesSum[idxProc + 1].suf,  myRightLimit);

                    mpiassert( MPI_Irecv(&buffer[indexArray], (endIndex - firstIndex) * sizeof(SortType), MPI_BYTE, idxProc, 0, currentComm, &requests[iterRequest++]),  __LINE__ );
                    indexArray += endIndex - firstIndex;
                    ++idxProc;
                }

                MpiChangeGroup(currentGroup, currentComm, splitProc + 1, currentNbProcs - 1);
            }

            // Proceed all send/receive
            mpiassert( MPI_Waitall(iterRequest, requests, NULL),  __LINE__ );


            // Copy res into outputArray
            if(outputSize < size){
                delete[] outputArray;
                outputArray = new SortType[size];
                outputSize = size;
            }

            memcpy(outputArray, buffer, sizeof(SortType) * size);
        }

        /////////////////////////////////////////////////
        // End QsMpi sort
        /////////////////////////////////////////////////

        // Clean
        delete[] buffer;
        MPI_Comm_free(&currentComm);
        MPI_Group_free(&currentGroup);

        // Normal Quick sort
        QsOmp<SortType,PivotType>(outputArray, size);
        outputSize = size;
    }

};

#endif // FQUICKSORT_HPP
