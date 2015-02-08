#ifndef FSTARPUUTILS_HPP
#define FSTARPUUTILS_HPP

extern "C"{
#include <starpu.h>
}

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

#endif // FSTARPUUTILS_HPP

