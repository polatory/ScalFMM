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
#ifndef FTIC_HPP
#define FTIC_HPP

#include "FGlobal.hpp"
#include <time.h>

// We need an os specific function
#if defined(WINDOWS)
    #include <windows.h>
#else
    #ifndef POSIX
        #warning Posix used withoug being explicitly defined
    #endif
    #include <sys/time.h>
    #include <unistd.h>
    #include <stdint.h>
#endif


/** @author Berenger Bramas (berenger.bramas@inria.fr)
  * This class is a easy to use time counter
  * With it you can very easyli estimate elapsed time between two moments
  * <code>
  * FTic counter;<br>
  * counter.tic();<br>
  * //...<br>
  * counter.tac();<br>
  * counter.elapsed(); //> time in s<br>
  * </code>
  *
  * The special method that uses asm register is based on the code by : Alexandre DENIS
  * http://dept-info.labri.fr/~denis/Enseignement/2006-SSECPD/timing.h
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
      * gettimeofday on linux or a direct ASM method
      */
    static double GetTime(){
#ifdef WINDOWS
        return static_cast<double>(GetTickCount())/1000.0;
#else // We are in linux/posix
        // better than gettimeofday but need -lrt
        struct timespec tp;
        clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
        return double(tp.tv_sec) + (double(tp.tv_nsec)/1000000000.0);
#endif
    }
};


#endif

