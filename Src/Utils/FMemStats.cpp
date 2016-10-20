// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#include "FMemStats.h"

FMemStats FMemStats::controler;

#include <cstdio>
#include <cstdlib>

#ifdef SCALFMM_USE_MEM_STATS
    // Regular scalar new
    void* operator new(std::size_t n) {
        void* const allocated = std::malloc(n + sizeof(size_t));
        if(allocated){
            *(static_cast<size_t*>(allocated)) = n;
            FMemStats::controler.allocate(n);
            return static_cast<unsigned char*>(allocated) + sizeof(size_t);
        }
        throw std::bad_alloc();
        return allocated;
    }

    void* operator new[]( std::size_t n ) {
        void* const allocated = std::malloc(n + sizeof(size_t));
        if(allocated){
            *(static_cast<size_t*>(allocated)) = n;
            FMemStats::controler.allocate(n);
            return static_cast<unsigned char*>(allocated) + sizeof(size_t);
        }
        throw std::bad_alloc();
        return allocated;
    }

    void* operator new  ( std::size_t n, const std::nothrow_t& tag){
        void* const allocated = std::malloc(n + sizeof(size_t));
        if(allocated){
            *(static_cast<size_t*>(allocated)) = n;
            FMemStats::controler.allocate(n);
            return static_cast<unsigned char*>(allocated) + sizeof(size_t);
        }
        return allocated;
    }

    void* operator new[] ( std::size_t n, const std::nothrow_t& tag){
        void* const allocated = std::malloc(n + sizeof(size_t));
        if(allocated){
            *(static_cast<size_t*>(allocated)) = n;
            FMemStats::controler.allocate(n);
            return static_cast<unsigned char*>(allocated) + sizeof(size_t);
        }
        return allocated;
    }

    // Regular scalar delete
    void operator delete(void* p) noexcept{
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - sizeof(size_t))) );
            std::free(static_cast<unsigned char*>(p) - sizeof(size_t));
        }
    }

    void operator delete[](void* p) noexcept{
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - sizeof(size_t))) );
            std::free(static_cast<unsigned char*>(p) - sizeof(size_t));
        }
    }

    void operator delete  ( void* p, const std::nothrow_t& /*tag*/) {
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - sizeof(size_t))) );
            std::free(static_cast<unsigned char*>(p) - sizeof(size_t));
        }
    }

    void operator delete[]( void* p, const std::nothrow_t& /*tag*/) {
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - sizeof(size_t))) );
            std::free(static_cast<unsigned char*>(p) - sizeof(size_t));
        }
    }

#endif
