#ifndef FMEMSTATS_H
#define FMEMSTATS_H

#include "FGlobal.hpp"

#include <stdlib.h>
#include <cstring>

#ifdef SCALFMM_USE_MEM_STATS
#include <new>
#include <stdexcept>
#warning You are using meme stats
    void* operator new(std::size_t n) throw(std::bad_alloc);
    void* operator new(size_t n, std::nothrow_t const&) throw();
    void* operator new[](size_t n) throw(std::bad_alloc);
    void* operator new[](size_t n, std::nothrow_t const&) throw();
    void operator delete(void* p) throw();
    void operator delete(void* p, std::nothrow_t const&) throw();
    void operator delete[](void* p) throw();
    void operator delete[](void* p, std::nothrow_t const&) throw();
#endif

class FMemStats {
private:
    unsigned long long maxAllocated;
    unsigned long long totalAllocated;
    std::size_t currentAllocated;

    FMemStats()
        : maxAllocated(0), totalAllocated(0), currentAllocated(0) {
    }

    void allocate(const std::size_t size){
        currentAllocated += size;
        totalAllocated += size;

        if(maxAllocated < currentAllocated){
            maxAllocated = currentAllocated;
        }
    }

    void deallocate(const std::size_t size){
        currentAllocated -= size;
    }

#ifdef SCALFMM_USE_MEM_STATS
    friend void* operator new(std::size_t n) throw(std::bad_alloc);
    friend void* operator new(size_t n, std::nothrow_t const&) throw();
    friend void* operator new[](size_t n) throw(std::bad_alloc);
    friend void* operator new[](size_t n, std::nothrow_t const&) throw();
    friend void operator delete(void* p) throw();
    friend void operator delete(void* p, std::nothrow_t const&) throw();
    friend void operator delete[](void* p) throw();
    friend void operator delete[](void* p, std::nothrow_t const&) throw();
#endif

public:
    static FMemStats controler;

    unsigned long long getMaxAllocated() const{
        return maxAllocated;
    }

    unsigned long long getTotalAllocated() const{
        return totalAllocated;
    }

    std::size_t getCurrentAllocated() const{
        return currentAllocated;
    }

    float getMaxAllocatedMB() const{
        return float(getMaxAllocated()) / 1024 / 1024;
    }

    float getTotalAllocatedMB() const{
        return float(getTotalAllocated()) / 1024 / 1024;
    }

    float getCurrentAllocatedMB() const{
        return float(getCurrentAllocated()) / 1024 / 1024;
    }

    bool isUsed() const {
#ifdef SCALFMM_USE_MEM_STATS
        return true;
#else
        return false;
#endif
    }
};


#endif // FMEMSTATS_H
