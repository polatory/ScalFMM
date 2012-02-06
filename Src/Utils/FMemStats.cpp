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
#include "FMemStats.h"

FMemStats FMemStats::controler;

#include <cstdio>

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
            free(static_cast<unsigned char*>(p) - 8);
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
