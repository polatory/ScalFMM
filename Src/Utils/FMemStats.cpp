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
#include <cstdlib>

#ifdef ScalFMM_USE_MEM_STATS
    // Regular scalar new
    void* operator new(std::size_t n) {
        void* const allocated = std::malloc(n + 8);
        if(allocated){
            *(static_cast<size_t*>(allocated)) = n;
            FMemStats::controler.allocate(n);
            return static_cast<unsigned char*>(allocated) + 8;
        }
        return allocated;
    }

    // Regular scalar delete
    void operator delete(void* p) noexcept{
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - 8)) );
            std::free(static_cast<unsigned char*>(p) - 8);
        }
    }

    void operator delete[](void* p) noexcept{
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - 8)) );
            std::free(static_cast<unsigned char*>(p) - 8);
        }
    }

    void operator delete  ( void* p, const std::nothrow_t& /*tag*/) {
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - 8)) );
            std::free(static_cast<unsigned char*>(p) - 8);
        }
    }

    void operator delete[]( void* p, const std::nothrow_t& /*tag*/) {
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - 8)) );
            std::free(static_cast<unsigned char*>(p) - 8);
        }
    }

#endif
