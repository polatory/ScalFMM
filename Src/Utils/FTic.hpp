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
    #if defined(__i386__) || defined(__pentium__) || defined(__pentiumpro__) || defined(__i586__) || defined(__i686__) || defined(__k6__) || defined(__k7__) || defined(__x86_64__)
        #define FSPECIAL_TIMER
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
#ifndef FSPECIAL_TIMER
        timeval t;
        gettimeofday(&t, NULL);
        return double(t.tv_sec) + (double(t.tv_usec)/1000000.0);
#else //FSPECIAL_TIMER is defined
        typedef union u_tick {
          uint64_t tick;
          struct {
            uint32_t low;
            uint32_t high;
          } sub;
        } tick_t;

        tick_t counter;
        __asm__ volatile("rdtsc" : "=a" (counter.sub.low), "=d" (counter.sub.high));
        return double(counter.tick) * Scale;
#endif
#endif
    }

#ifdef FSPECIAL_TIMER
    static const double Scale;
#endif
};


#endif

