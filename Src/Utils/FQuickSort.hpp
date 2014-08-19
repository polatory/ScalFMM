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
#include <vector> // For parallel without task

#include "FGlobal.hpp"
#include "FMemUtils.hpp"
#include "FTrace.hpp"

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
        Swap(array[right],array[((right - left ) / 2) + left]);

        IndexType idx = left;
        while( idx < right && CompareType(array[idx]) <= CompareType(array[right])){
            idx += 1;
        }
        left = idx;

        for( ; idx < right ; ++idx){
            if( CompareType(array[idx]) <= CompareType(array[right]) ){
                Swap(array[idx],array[left]);
                left += 1;
            }
        }

        Swap(array[left],array[right]);

        return left;
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
            const IndexType part = QsPartition(array, left, right);
            if( deep ){
                #pragma omp task
                QsOmpTask(array,part + 1,right, deep - 1);
                #pragma omp task
                QsOmpTask(array,left,part - 1, deep - 1);
            }
            else {
                QsSequentialStep(array,part + 1,right);
                QsSequentialStep(array,left,part - 1);
            }
        }
    }

public:
    /* a sequential qs */
    static void QsSequential(SortType array[], const IndexType size){
        QsSequentialStep(array,0,size-1);
    }

#if _OPENMP < 200805 || defined(__ICC) || defined(__INTEL_COMPILER)
    class TaskInterval{
        IndexType left;
        IndexType right;
        int deep;
    public:
        TaskInterval(const IndexType inLeft, const IndexType inRight, const int inDeep)
            : left(inLeft), right(inRight), deep(inDeep){
        }

        IndexType getLeft() const{
            return left;
        }
        IndexType getRight() const{
            return right;
        }
        int getDeep() const{
            return deep;
        }
    };

    static void QsOmp(SortType elements[], const int nbElements){
        const int nbTasksRequiere = (omp_get_max_threads() * 5);
        int deep = 0;
        while( (1 << deep) < nbTasksRequiere ) deep += 1;

        std::vector<TaskInterval> tasks;
        tasks.push_back(TaskInterval(0, nbElements-1, deep));
        int numberOfThreadProceeding = 0;
        omp_lock_t mutexShareVariable;
        omp_init_lock(&mutexShareVariable);

        #pragma omp parallel
        {
            bool hasWorkToDo = true;
            while(hasWorkToDo){
                // Ask for the mutex
                omp_set_lock(&mutexShareVariable);
                if(tasks.size()){
                    // There is tasks to proceed
                    const TaskInterval ts(tasks.back());
                    tasks.pop_back();

                    // Does this task should create some other?
                    if(ts.getDeep() == 0){
                        // No release the mutex and run in seq
                        omp_unset_lock(&mutexShareVariable);
                        QsSequentialStep(elements , ts.getLeft(), ts.getRight());
                    }
                    else{
                        // Yes so inform other and release the mutex
                        numberOfThreadProceeding += 1;
                        omp_unset_lock(&mutexShareVariable);

                        // Partition
                        const IndexType part = QsPartition(elements, ts.getLeft(), ts.getRight());

                        // Push the new task in the vector
                        omp_set_lock(&mutexShareVariable);
                        tasks.push_back(TaskInterval(part+1, ts.getRight(), ts.getDeep()-1));
                        tasks.push_back(TaskInterval(ts.getLeft(), part-1, ts.getDeep()-1));
                        // We create new task but we are not working so inform other
                        numberOfThreadProceeding -= 1;
                        omp_unset_lock(&mutexShareVariable);
                    }
                }
                else{
                    // There is not task in the vector
                    #pragma omp flush(numberOfThreadProceeding)
                    if(numberOfThreadProceeding == 0){
                        // And there is no thread that may create some tasks so stop here
                        hasWorkToDo = false;
                    }
                    // Release mutex
                    omp_unset_lock(&mutexShareVariable);
                }
            }
        }

        omp_destroy_lock(&mutexShareVariable);
    }
#else
    /** The openmp quick sort */
    static void QsOmp(SortType array[], const IndexType size){
        const int nbTasksRequiere = (omp_get_max_threads() * 5);
        int deep = 0;
        while( (1 << deep) < nbTasksRequiere ) deep += 1;

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                QsOmpTask(array, 0, size - 1 , deep);
            }
        }
    }
#endif
};

#endif // FQUICKSORT_HPP
