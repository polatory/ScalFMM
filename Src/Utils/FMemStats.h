#ifndef FMEMSTATS_H
#define FMEMSTATS_H


#include <new>
#include <stdexcept>
#include <stdlib.h>

#ifdef SCALFMM_USE_MEM_STATS
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

    friend void* operator new(std::size_t n) throw(std::bad_alloc);
    friend void* operator new(size_t n, std::nothrow_t const&) throw();
    friend void* operator new[](size_t n) throw(std::bad_alloc);
    friend void* operator new[](size_t n, std::nothrow_t const&) throw();
    friend void operator delete(void* p) throw();
    friend void operator delete(void* p, std::nothrow_t const&) throw();
    friend void operator delete[](void* p) throw();
    friend void operator delete[](void* p, std::nothrow_t const&) throw();

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
};


#endif // FMEMSTATS_H
