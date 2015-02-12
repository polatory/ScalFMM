#ifndef FSTARPUUTILS_HPP
#define FSTARPUUTILS_HPP

/////////////////////////////////////////////////////


extern "C"{
#include <starpu.h>
}

/////////////////////////////////////////////////////

#if (STARPU_MAJOR_VERSION >= 1) && (STARPU_MINOR_VERSION >= 2)
#define STARPU_SUPPORT_COMMUTE
#endif

/////////////////////////////////////////////////////

#include <functional>

class FStarPUUtils{
protected:
    static void ExecOnWorkersBind(void* ptr){
        std::function<void(void)>* func = (std::function<void(void)>*) ptr;
        (*func)();
    }

public:
    static void ExecOnWorkers(const unsigned int inWorkersType, std::function<void(void)> func){
        starpu_execute_on_each_worker(ExecOnWorkersBind, &func, inWorkersType);
    }
};

/////////////////////////////////////////////////////

#ifndef STARPU_SUPPORT_COMMUTE
    #define STARPU_COMMUTE STARPU_NONE
#endif

#endif // FSTARPUUTILS_HPP

