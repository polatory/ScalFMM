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
            #pragma omp flush(sense)
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
