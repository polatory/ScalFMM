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
#ifndef FOMPBARRIER_HPP
#define FOMPBARRIER_HPP

#include <omp.h>
#include <climits>

/** This function is a custom omp barrier
  * Because openmp give only a global barrier we need
  * to be ablo to peform a barrier operation between a group
  * of thread only.
  */

class FOmpBarrier {
private:
    int nbThreads;          //<The number of threads for this barrier
    int currentNbThread;    //<The current number of threads waiting
    bool sense;             //<Direct barrier feedback protection
    omp_lock_t mutex;       //<To have an atomic int

    FOmpBarrier(FOmpBarrier&){}
    FOmpBarrier& operator=(FOmpBarrier&){return *this;}

public:
    /** Constructor with the number of threads */
    FOmpBarrier(const int inNbThreads = INT_MAX)
        : nbThreads(inNbThreads), currentNbThread(0), sense(false) {
        omp_init_lock( &mutex );
    }

    /** Destructor, release the omp lock */
    ~FOmpBarrier(){
        omp_destroy_lock( &mutex );
    }

    /** Perform a barrier */
    void wait(){
        const bool mySense = sense;
        omp_set_lock( &mutex );
        const int nbThreadsArrived = (++currentNbThread);
        omp_unset_lock( &mutex );

        if(nbThreadsArrived == nbThreads) {
            currentNbThread = 0;
            sense = !sense;
            #pragma omp flush
        }
        else {
            volatile const bool* const ptSense = &sense;
            while( (*ptSense) == mySense){
            }
        }
    }

    /** Change the number of threads */
    void setNbThreads(const int inNbThread){
        omp_set_lock( &mutex );
        nbThreads = inNbThread;
        omp_unset_lock( &mutex );
    }
};


#endif // FOMPBARRIER_HPP
