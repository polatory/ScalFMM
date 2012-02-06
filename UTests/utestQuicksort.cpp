// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
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



