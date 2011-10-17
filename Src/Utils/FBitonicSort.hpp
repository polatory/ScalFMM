#ifndef BITONICSORT_HPP
#define BITONICSORT_HPP

#include <cstdlib>
#include <cmath>

#include <mpi.h>

#include "FTrace.hpp"
#include "FMpi.hpp"
#include "FQuickSort.hpp"

class FBitonicSort {
private:

    ////////////////////////////////////////////////////////////////
    // Bitonic parallel sort !
    ////////////////////////////////////////////////////////////////

    // Mpi flag
    static const int FlagMin = 5;
    static const int FlagMax = 6;
    static const int FlagMinMess = 4;
    static const int FlagMaxMess = 3;

    // This function exchange data with the other rank,
    // its send the max value and receive min value
    template <class SortType, class CompareType, class IndexType>
    static void sendMaxAndGetMin(SortType array[], const IndexType size, const int otherRank){
        IndexType left  = -1;
        IndexType right = size - 1;
        IndexType pivot = left + (right - left + 1)/2;
        CompareType otherValue = -1;
        CompareType tempCompareValue = CompareType(array[pivot]);
        MPI_Sendrecv(&tempCompareValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMin,&otherValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMax,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        while( pivot != left && pivot != right  && array[pivot] != otherValue) {

            if( array[pivot] < otherValue ){
                left = pivot;
            }
            else {
                right = pivot;
            }
            pivot = left + (right - left + 1)/2;
            tempCompareValue = CompareType(array[pivot]);

            MPI_Sendrecv(&tempCompareValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMin,&otherValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMax,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        if( otherValue <= array[pivot] ){
            MPI_Sendrecv_replace(&array[pivot], int((size - pivot) * sizeof(SortType)) , MPI_BYTE,
                                   otherRank, FlagMinMess, otherRank, FlagMaxMess,
                                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
        else if( array[pivot] < otherValue){
            if(pivot != size - 1){
                MPI_Sendrecv_replace(&array[pivot + 1], int((size - pivot - 1) * sizeof(SortType)) , MPI_BYTE,
                                       otherRank, FlagMinMess, otherRank, FlagMaxMess,
                                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

    }

    // This function exchange data with the other rank,
    // its send the min value and receive max value
    template <class SortType, class CompareType, class IndexType>
    static void sendMinAndGetMax(SortType array[], const IndexType size, const int otherRank){
        IndexType left  = 0;
        IndexType right = size ;
        IndexType pivot = left + (right - left)/2;
        CompareType otherValue = -1;
        CompareType tempCompareValue = CompareType(array[pivot]);
        MPI_Sendrecv(&tempCompareValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMax,&otherValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMin,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        while(  pivot != left  && array[pivot] != otherValue) {

            if( array[pivot] < otherValue ){
                left = pivot;
            }
            else {
                right = pivot;
            }
            pivot = left + (right - left)/2;
            tempCompareValue = CompareType(array[pivot]);
            MPI_Sendrecv(&tempCompareValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMax,&otherValue,sizeof(CompareType),MPI_BYTE,otherRank,FlagMin,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }


        if( array[pivot] <= otherValue ){
            MPI_Sendrecv_replace(&array[0], int((pivot + 1) * sizeof(SortType)) , MPI_BYTE,
                                   otherRank, FlagMaxMess, otherRank, FlagMinMess,
                                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if( otherValue < array[pivot]){
            if(pivot != 0){
                MPI_Sendrecv_replace(&array[0], int((pivot) * sizeof(SortType)) , MPI_BYTE,
                                       otherRank, FlagMaxMess, otherRank, FlagMinMess,
                                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }


public:

    /*
    From :
    http://web.mst.edu/~ercal/387/P3/pr-proj-3.pdf

    Parallel Bitonic Sort Algorithm for processor Pk (for k := 0 . . . P − 1)
    d:= log P
    // cube dimension
    sort(local − datak ) // sequential sort
    // Bitonic Sort follows
    for i:=1 to d do
        window-id = Most Significant (d-i) bits of Pk
        for j:=(i-1) down to 0 do
            if((window-id is even AND j th bit of Pk = 0)
            OR (window-id is odd AND j th bit of Pk = 1))
                then call CompareLow(j)
            else
                call CompareHigh(j)
            endif
        endfor
    endfor
      */
    template <class SortType, class CompareType, class IndexType>
    static void sort(SortType array[], const IndexType size, const int np, const int rank){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Bitonic" , __FILE__ , __LINE__) );

        FQuickSort<SortType,CompareType,IndexType>::QsOmp(array, size);

        const int logNp = int(log2(np));
        for(int bitIdx = 1 ; bitIdx <= logNp ; ++bitIdx){
            // window-id = Most Significant (d-i) bits of Pk
            const int diBit =  (rank >> bitIdx) & 0x1;

            for(int otherBit = bitIdx - 1 ; otherBit >= 0 ; --otherBit){
                // if((window-id is even AND j th bit of Pk = 0)
                // OR (window-id is odd AND j th bit of Pk = 1))

                const int myOtherBit = (rank >> otherBit) & 0x1;
                const int otherRank = rank ^ (1 << otherBit);

                if( diBit != myOtherBit ){
                    sendMinAndGetMax<SortType,CompareType>(array, size, otherRank);
                }
                else{
                    sendMaxAndGetMin<SortType,CompareType>(array, size, otherRank);
                }
                // A merge sort is possible since the array is composed
                // by two part already sorted, but we want to do this in space
                FQuickSort<SortType,CompareType,IndexType>::QsOmp(array, size);
            }
        }
    }

    template <class SortType, class CompareType, class IndexType>
    static void sort(SortType array[], const IndexType size){
        int rank = 0;
        int nprocs = 0;

        MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        sort<SortType, CompareType>(array, size, nprocs, rank);
    }
};

#endif // BITONICSORT_HPP
