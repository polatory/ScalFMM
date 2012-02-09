// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#include "FUTester.hpp"
#include "../Src/Utils/FQuickSort.hpp"

#include <unistd.h>

/**
* This file is a unit test for the quick sort
*/


/** this class test the list container */
class TestQuickSort : public FUTester<TestQuickSort> {
    static bool IsSorted(long long array[], const long size){
        for(int idx = 1; idx < size ; ++idx){
            if(array[idx-1] > array[idx]){
                return false;
            }
        }
        return true;
    }

    void manyThreads(){
        const long Size = 100000;
        long long* const array = new long long[Size];
        srand(0);
        const int originalThreadsNumber = omp_get_num_threads();

        for(int idxThread = 1 ; idxThread <= omp_get_max_threads() ; idxThread *= 2){
            omp_set_num_threads(idxThread);

            for(long idx = 0 ; idx < Size ; ++idx){
                array[idx] = rand();
            }

            FQuickSort<long long, long long, long>::QsOmp(array, Size);

            assert(IsSorted(array,Size));
        }

        omp_set_num_threads(originalThreadsNumber);
        delete [] array;
    }

    void bigSize(){
        const long Size = 10000000;//100000000;
        long long* const array = new long long[Size];

        for(long idx = 0 ; idx < Size ; ++idx){
            array[idx] = rand();
        }

        FQuickSort<long long, long long, long>::QsOmp(array, Size);
        assert(IsSorted(array,Size));

        delete [] array;
    }

    // set test
    void SetTests(){
        AddTest(&TestQuickSort::manyThreads,"Many threads");
        AddTest(&TestQuickSort::bigSize,"Big sort");
    }
};

// You must do this
TestClass(TestQuickSort)



