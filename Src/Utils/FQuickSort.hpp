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
#ifndef FQUICKSORT_HPP
#define FQUICKSORT_HPP

#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>

#include "FGlobal.hpp"
#include "FMemUtils.hpp"
#include "FTrace.hpp"

#include "FOmpBarrier.hpp"

/** This class is parallel quick sort
  * It hold a mpi version
  * + 2 openmp versions (one on tasks and the other like mpi)
  * + a sequential version
  *
  * The task based algorithm is easy to undestand,
  * for the mpi/openmp2nd please see
  * Introduction to parallel computing (Grama Gupta Karypis Kumar)
  */

template <class SortType, class CompareType, class IndexType>
class FQuickSort {
protected:
    ////////////////////////////////////////////////////////////
    // Miscialenous functions
    ////////////////////////////////////////////////////////////

    template <class T>
    static T GetLeft(const T inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        return T(ceil(step * inIdProc));
    }

    /** Compute a right index from data */
    template <class T>
    static T GetRight(const T inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        const T res = T(ceil(step * (inIdProc+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    /** Compute a proc id from index & data */
    template <class T>
    static int GetProc(const T position, const T inSize, const int inNbProc) {
        const double step = double(inSize) / double(inNbProc);
        return int(double(position)/step);
    }

    /** swap to value */
    template <class NumType>
    static inline void Swap(NumType& value, NumType& other){
        NumType temp = value;
        value = other;
        other = temp;
    }

    ////////////////////////////////////////////////////////////
    // Quick sort
    ////////////////////////////////////////////////////////////

    /* Use in the sequential qs */
    static IndexType QsPartition(SortType array[], IndexType left, IndexType right){
        const IndexType part = right;
        Swap(array[part],array[((right - left ) / 2) + left]);
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

    /* A local iteration of qs */
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

    /* The sequential qs */
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

    /* The openmp qs */
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

            IndexType myLeft = GetLeft(size, myThreadId, omp_get_num_threads());
            IndexType myRight = GetRight(size, myThreadId, omp_get_num_threads()) - 1;

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
                int splitProc = GetProc(sufoffset - startIndex, nbElements, lastProc - firstProc + 1) + firstProc;
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

                myLeft = GetLeft(endIndex - startIndex + 1, myThreadId - firstProc, lastProc - firstProc + 1) + startIndex;
                myRight = GetRight(endIndex - startIndex + 1, myThreadId - firstProc, lastProc - firstProc + 1) + startIndex - 1;
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


};

#endif // FQUICKSORT_HPP
