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
#warning You are using mem stats
void* operator new(std::size_t n);
void operator delete(void* p) noexcept;
void operator delete[](void* p) noexcept;
void operator delete  ( void* ptr, const std::nothrow_t& tag);
void operator delete[]( void* ptr, const std::nothrow_t& tag);
#endif

/** Give the memory allocation details
  *
  */
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
    friend void* operator new(std::size_t n);
    friend void operator delete(void* p) noexcept;
    friend void operator delete[](void* p) noexcept;
    friend void operator delete  ( void* ptr, const std::nothrow_t& tag);
    friend void operator delete[]( void* ptr, const std::nothrow_t& tag);
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
