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
#include "FMemStats.h"

FMemStats FMemStats::controler;

#include <cstdio>

#ifdef ScalFMM_USE_MEM_STATS
    // Regular scalar new
    void* operator new(std::size_t n)
    {
        using namespace std;

        for (;;) {
            void* allocated_memory = ::operator new(n, nothrow);
            if (allocated_memory != 0){
                return allocated_memory;
            }

            // Store the global new handler
            new_handler global_handler = set_new_handler(0);
            set_new_handler(global_handler);

            if (global_handler) {
                global_handler();
            } else {
                throw bad_alloc();
            }
        }
    }

    // Nothrow scalar new
    void* operator new(size_t n,const std::nothrow_t& nothrow_value) noexcept
    {
        //if (n == 0) n = 1;
        void* const allocated = malloc(n + 8);
        if(allocated){
            *(static_cast<size_t*>(allocated)) = n;
            FMemStats::controler.allocate(n);
            return static_cast<unsigned char*>(allocated) + 8;
        }
        return allocated;
    }

    // Regular array new
    void* operator new[](size_t n)
    {
        return ::operator new(n);
    }

    // Nothrow array new
    void* operator new[](size_t n,const std::nothrow_t& nothrow_value) noexcept
    {
        return ::operator new(n, std::nothrow);
    }

    // Regular scalar delete
    void operator delete(void* p) noexcept {
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - 8)) );
            free(static_cast<unsigned char*>(p) - 8);
        }
    }

    // Nothrow scalar delete
    void operator delete(void* p,const std::nothrow_t& nothrow_value) noexcept {
        ::operator delete(p);
    }

    // Regular array delete
    void operator delete[](void* p) noexcept
    {
        ::operator delete(p);
    }

    // Nothrow array delete
    void operator delete[](void* p,const std::nothrow_t& nothrow_value) noexcept
    {
        ::operator delete(p);
    }

#endif
