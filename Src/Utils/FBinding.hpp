#ifndef FBINDING_HPP
#define FBINDING_HPP

#include "FGlobal.hpp"
#include "FAssert.hpp"

#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sched.h>
#include <sys/resource.h>
#include <thread>

namespace FBinding {

#if defined(linux) || defined(__linux__) || defined(unix) || defined(__unix__)
#define FBINDING_ENABLE
#endif

inline int GetThreadBinding(){
    // Mask will contain the current affinity
#ifdef FBINDING_ENABLE
    unsigned long mask = 0;
    // We need the thread pid (even if we are in openmp)
    pid_t tid = (pid_t) syscall(SYS_gettid);
    // Get the affinity
    FAssertLF(sched_getaffinity(tid, sizeof(mask), (cpu_set_t*)&mask) != -1);
    return int(mask>>1);
#endif
    return -1;
}

inline void SetThreadBinding(const int procId){
#ifdef FBINDING_ENABLE
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(procId, &set);

    pid_t tid = (pid_t) syscall(SYS_gettid);
    FAssertLF(sched_setaffinity(tid, sizeof(set), &set) != -1);
#endif
}

inline void BindThreadToAnyProcs(){
#ifdef FBINDING_ENABLE
    cpu_set_t set;
    CPU_ZERO(&set);

    const int nbProcs = (int)std::thread::hardware_concurrency();
    for( int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
        CPU_SET(idxProc, &set);
    }

    pid_t tid = (pid_t) syscall(SYS_gettid);
    FAssertLF(sched_setaffinity(tid, sizeof(set), &set) != -1);
#endif
}

}

#endif // FBINDING_HPP

