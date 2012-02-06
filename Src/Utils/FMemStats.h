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
#ifndef FMEMSTATS_H
#define FMEMSTATS_H

#include "FGlobal.hpp"

#include <cstdlib>
#include <cstring>


/** Memstat has to be enabled in the cmake,
  * then it will use the know each allocate and deallocate
  * and give simple stats like max, total used, current used
  */

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
    /** Singleton */
    static FMemStats controler;

    /** return the max that has been allocated */
    unsigned long long getMaxAllocated() const{
        return maxAllocated;
    }

    /** return the total memory allocated during the running */
    unsigned long long getTotalAllocated() const{
        return totalAllocated;
    }

    /** return the current size allcoated */
    std::size_t getCurrentAllocated() const{
        return currentAllocated;
    }

    /** get the max in MB */
    float getMaxAllocatedMB() const{
        return float(getMaxAllocated()) / 1024 / 1024;
    }

    /** get the total in MB */
    float getTotalAllocatedMB() const{
        return float(getTotalAllocated()) / 1024 / 1024;
    }

    /** get the current in MB */
    float getCurrentAllocatedMB() const{
        return float(getCurrentAllocated()) / 1024 / 1024;
    }

    /** To know if mem stat has been enabled */
    bool isUsed() const {
#ifdef SCALFMM_USE_MEM_STATS
        return true;
#else
        return false;
#endif
    }
};


#endif // FMEMSTATS_H
