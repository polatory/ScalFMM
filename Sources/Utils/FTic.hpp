#ifndef FTIC_HPP
#define FTIC_HPP

#if defined(_WIN32) || defined(ming)
    #define WINDOWS
    #include <windows.h>
#else
    #define POSIX
    #include <sys/time.h>
#endif

#include <time.h>

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
    double start;   //< start time (tic)
    double end;     //<stop time (tac)

public:
    /** Constructor */
    FTic() : start(0.0), end(0.0) {
    }

    /** Tic : start <= current time */
    void tic(){
        this->start = FTic::GetTime();
    }

    /** Tac : end <= current time */
    void tac(){
        this->end = FTic::GetTime();
    }

    /** Return end - start
      * @return the time elapsed between tic & tac in second */
    double elapsed() const{
        return this->end - this->start;
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

