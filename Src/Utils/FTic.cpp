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
#include "FTic.hpp"

// The fallowing code is only needed when using the special timer
#ifdef FSPECIAL_TIMER
// This function is used to get the scale between tick and second
double GetTickScale() {
        typedef union u_tick {
          uint64_t tick;
          struct {
            uint32_t low;
            uint32_t high;
          } sub;
        } tick_t;

        tick_t tickCounterStart, tickCounterEnd;
        timeval todStart, todEnd;

        __asm__ volatile("rdtsc" : "=a" (tickCounterStart.sub.low), "=d" (tickCounterStart.sub.high));
        gettimeofday(&todStart,0);
        usleep(500000);
        __asm__ volatile("rdtsc" : "=a" (tickCounterEnd.sub.low), "=d" (tickCounterEnd.sub.high));
        gettimeofday(&todEnd,0);
        return ( (double(todEnd.tv_sec) + (double(todEnd.tv_usec)/1e6)) - (double(todStart.tv_sec) + (double(todStart.tv_usec)/1e6)) ) /
                          double(tickCounterEnd.tick - tickCounterStart.tick);
}
// For the special timer
const double FTic::Scale = GetTickScale();
#endif
