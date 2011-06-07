#ifndef FTIC_HPP
#define FTIC_HPP

#include <time.h>

#include "FGlobal.hpp"

// We need an os specific function
#if defined(WINDOWS)
    #include <windows.h>
#else
    #ifndef POSIX
        #warning Posix used withoug being explicitly defined
    #endif
    #include <sys/time.h>
#endif


/** @author Berenger Bramas (berenger.bramas@inria.fr)
  * This class is a easy to use time counter
  * With it you can very easyli estimate elapsed
  * time between two moments
  * <code>
  * FTic counter;<br>
  * counter.tic();<br>
  * //...<br>
  * counter.tac();<br>
  * counter.elapsed(); //> time in s<br>
  * </code>
  */
class FTic {
private:
    double start;       //< start time (tic)
    double end;         //< stop time (tac)
    double cumulate;    //< the cumulate time

public:
    /** Constructor */
    FTic() : start(0.0), end(0.0), cumulate(0.0) {
        tic();
    }

    /** Tic : start <= current time */
    void tic(){
        this->start = FTic::GetTime();
    }

    /** Tac : end <= current time */
    void tac(){
        this->end = FTic::GetTime();
        cumulate += elapsed();
    }

    /** Return end - start
      * @return the time elapsed between tic & tac in second */
    double elapsed() const{
        return this->end - this->start;
    }

    /** Return cumulate
      * @return the time elapsed between ALL tic & tac in second */
    double cumulated() const{
        return cumulate;
    }

    /** Return end - start
      * @return the time elapsed between tic & tac in second */
    double tacAndElapsed() {
        tac();
        return elapsed();
    }

    /** Global get time
      * @return a global time
      * GetTickCount on windows
      * gettimeofday on linux
      */
    static double GetTime(){
#ifdef WINDOWS
        return static_cast<double>(GetTickCount())/1000.0;
#else
        timeval t;
        gettimeofday(&t, NULL);
        return (t.tv_sec) + (t.tv_usec/1000000.0);
#endif
    }
};


#endif

