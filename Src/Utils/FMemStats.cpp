#include "FMemStats.h"

FMemStats FMemStats::controler;


#ifdef SCALFMM_USE_MEM_STATS
    // Regular scalar new
    void* operator new(std::size_t n) throw(std::bad_alloc)
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
    void* operator new(size_t n, std::nothrow_t const&) throw()
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
    void* operator new[](size_t n) throw(std::bad_alloc)
    {
        return ::operator new(n);
    }

    // Nothrow array new
    void* operator new[](size_t n, std::nothrow_t const&) throw()
    {
        return ::operator new(n, std::nothrow);
    }

    // Regular scalar delete
    void operator delete(void* p) throw(){
        if(p){
            FMemStats::controler.deallocate( *(reinterpret_cast<size_t*>(static_cast<unsigned char*>(p) - 8)) );
            free(reinterpret_cast<unsigned char*>(p) - 8);
        }
    }

    // Nothrow scalar delete
    void operator delete(void* p, std::nothrow_t const&) throw(){
        ::operator delete(p);
    }

    // Regular array delete
    void operator delete[](void* p) throw()
    {
        ::operator delete(p);
    }

    // Nothrow array delete
    void operator delete[](void* p, std::nothrow_t const&) throw()
    {
        ::operator delete(p);
    }

#endif