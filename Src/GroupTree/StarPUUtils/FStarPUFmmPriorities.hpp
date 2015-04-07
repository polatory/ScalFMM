// @SCALFMM_PRIVATE
#ifndef FSTARPUFMMPRIORITIES_HPP
#define FSTARPUFMMPRIORITIES_HPP

#include "../../Utils/FGlobal.hpp"

#include "FStarPUKernelCapacities.hpp"

/**
 * @brief The FStarPUFmmPriorities class
 * This class should have an static method to be called by hetero getPrio.
 */
#ifdef STARPU_SUPPORT_SCHEDULER

#include "FStarPUHeteoprio.hpp"

class FStarPUFmmPriorities{
    int prioP2M;
    int prioM2M;

    int prioP2MSend;
    int prioM2MSend;

    int prioM2L;
    int prioM2LExtern;
    int prioL2L;
    int prioL2P;
    int prioP2P;
    int prioP2PExtern;
    int prioM2LMpi;
    int prioP2PMpi;

    int treeHeight;

    FStarPUKernelCapacities* capacities;

    static FStarPUFmmPriorities controller;


    FStarPUFmmPriorities(){
    }

public:
    static FStarPUFmmPriorities& Controller(){
        return controller;
    }

    static void InitSchedulerCallback(unsigned sched_ctx_id,
                                      struct _starpu_heteroprio_center_policy_heteroprio *heteroprio){
        Controller().initSchedulerCallback(sched_ctx_id, heteroprio);
    }

    void init(struct starpu_conf* conf, const int inTreeHeight,
              FStarPUKernelCapacities* inCapacities){
        capacities = inCapacities;

        conf->sched_policy = &_starpu_sched_heteroprio_policy,
                initialize_heteroprio_center_policy_callback = &InitSchedulerCallback;

        treeHeight  = inTreeHeight;

        int incPrio = 0;

        prioP2MSend = incPrio++;
        prioP2M     = incPrio++;

        prioM2MSend = incPrio++;
        prioM2M     = incPrio++;

        prioM2L     = incPrio;
        prioM2LExtern = incPrio;
        prioM2LMpi  = incPrio++;

        prioL2L     = incPrio++;

        incPrio += (treeHeight-2)-1 // M2L is done treeHeight-2 times
                   +(treeHeight-3)-1; // L2L is done treeHeight-3 times

        prioP2P     = incPrio;
        prioP2PExtern = incPrio;
        prioP2PMpi  = incPrio++;

        prioL2P     = incPrio++;
        assert(incPrio == 6 + (treeHeight-2) + (treeHeight-3));
    }

    void initSchedulerCallback(unsigned /*sched_ctx_id*/,
                               struct _starpu_heteroprio_center_policy_heteroprio *heteroprio){
#ifdef STARPU_USE_CPU
        // CPU follows the real prio
        {
            int cpuCountPrio = 0;
            //prioP2MSend = 0;
            //prioP2M     = prioP2MSend+1;
            if(capacities->supportP2M(FSTARPU_CPU_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioP2MSend;
                heteroprio->buckets[prioP2MSend].valide_archs |= STARPU_CPU;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioP2M;
                heteroprio->buckets[prioP2M].valide_archs |= STARPU_CPU;
            }
            //prioM2MSend = prioP2M+1;
            //prioM2M     = prioM2MSend+1;
            assert(cpuCountPrio == prioM2MSend); // True if CPU support all TODO
            if(capacities->supportM2M(FSTARPU_CPU_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioM2MSend;
                heteroprio->buckets[prioM2MSend].valide_archs |= STARPU_CPU;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioM2M;
                heteroprio->buckets[prioM2M].valide_archs |= STARPU_CPU;
            }

            // prioM2L       = prioM2M+1;
            // prioM2LExtern = prioM2L;
            // prioM2LMpi    = prioM2L;
            // prioL2L     = prioM2L+1;
            assert(cpuCountPrio == prioM2L); // True if CPU support all TODO
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(capacities->supportM2L(FSTARPU_CPU_IDX)){
                    const int prioM2LAtLevel = getPrioM2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioM2LAtLevel;
                    heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_CPU;
                }
                if(idxLevel != treeHeight-1 && capacities->supportL2L(FSTARPU_CPU_IDX)){
                    const int prioL2LAtLevel = getPrioL2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioL2LAtLevel;
                    heteroprio->buckets[prioL2LAtLevel].valide_archs |= STARPU_CPU;
                }
            }
            assert(cpuCountPrio == prioP2P); // True if CPU support all TODO

            //prioP2P       = prioL2L + (treeHeight-3)*2+1 +1;
            //prioP2PExtern = prioP2P;
            //prioP2PMpi    = prioP2P;
            if(capacities->supportP2P(FSTARPU_CPU_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioP2P;
                heteroprio->buckets[prioP2P].valide_archs |= STARPU_CPU;
            }

            assert(cpuCountPrio == prioL2P); // True if CPU support all TODO
            //prioL2P     = prioP2PMpi+1;
            if(capacities->supportL2P(FSTARPU_CPU_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CPU_IDX][cpuCountPrio++] = prioL2P;
                heteroprio->buckets[prioL2P].valide_archs |= STARPU_CPU;
            }

            heteroprio->nb_prio_per_arch_index[FSTARPU_CPU_IDX] = unsigned(cpuCountPrio);
        }
#endif
#ifdef STARPU_USE_OPENCL
        {
            int openclCountPrio = 0;

            //prioP2P       = prioL2L + (treeHeight-3)*2+1 +1;
            //prioP2PExtern = prioP2P;
            //prioP2PMpi    = prioP2P;
            if(capacities->supportP2P(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioP2P;
                heteroprio->buckets[prioP2P].factor_base_arch_index = FSTARPU_OPENCL_IDX;
                heteroprio->buckets[prioP2P].valide_archs |= STARPU_OPENCL;
#ifdef STARPU_USE_CPU
                heteroprio->buckets[prioP2P].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
            }

            // prioM2L       = prioM2M+1;
            // prioM2LExtern = prioM2L;
            // prioM2LMpi    = prioM2L;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(capacities->supportM2L(FSTARPU_OPENCL_IDX)){
                    const int prioM2LAtLevel = getPrioM2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioM2LAtLevel;
                    heteroprio->buckets[prioM2LAtLevel].factor_base_arch_index = FSTARPU_OPENCL_IDX;
                    heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_OPENCL;
#ifdef STARPU_USE_CPU
                    heteroprio->buckets[prioM2LAtLevel].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
                }
            }

            //prioP2MSend = 0;
            //prioP2M     = prioP2MSend+1;
            if(capacities->supportP2M(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioP2MSend;
                heteroprio->buckets[prioP2MSend].valide_archs |= STARPU_OPENCL;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioP2M;
                heteroprio->buckets[prioP2M].valide_archs |= STARPU_OPENCL;
            }

            //prioM2MSend = prioP2M+1;
            //prioM2M     = prioM2MSend+1;
            if(capacities->supportM2M(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioM2MSend;
                heteroprio->buckets[prioM2MSend].valide_archs |= STARPU_OPENCL;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioM2M;
                heteroprio->buckets[prioM2M].valide_archs |= STARPU_OPENCL;
            }

            // prioL2L     = prioM2L+1;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(idxLevel != treeHeight-1 && capacities->supportL2L(FSTARPU_OPENCL_IDX)){
                    const int prioL2LAtLevel = getPrioL2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioL2LAtLevel;
                    heteroprio->buckets[prioL2LAtLevel].valide_archs |= STARPU_OPENCL;
                }
            }

            //prioL2P     = prioP2PMpi+1;
            if(capacities->supportL2P(FSTARPU_OPENCL_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_OPENCL_IDX][openclCountPrio++] = prioL2P;
                heteroprio->buckets[prioL2P].valide_archs |= STARPU_OPENCL;
            }

            heteroprio->nb_prio_per_arch_index[FSTARPU_OPENCL_IDX] = unsigned(openclCountPrio);
        }
#endif
#ifdef STARPU_USE_CUDA
        {
            int openclCountPrio = 0;

            //prioP2P       = prioL2L + (treeHeight-3)*2+1 +1;
            //prioP2PExtern = prioP2P;
            //prioP2PMpi    = prioP2P;
            if(capacities->supportP2P(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioP2P;
                heteroprio->buckets[prioP2P].valide_archs |= STARPU_CUDA;
                heteroprio->buckets[prioP2P].factor_base_arch_index = FSTARPU_CUDA_IDX;
#ifdef STARPU_USE_CPU
                heteroprio->buckets[prioP2P].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
            }

            // prioM2L       = prioM2M+1;
            // prioM2LExtern = prioM2L;
            // prioM2LMpi    = prioM2L;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(capacities->supportM2L(FSTARPU_CUDA_IDX)){
                    const int prioM2LAtLevel = getPrioM2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioM2LAtLevel;
                    heteroprio->buckets[prioM2LAtLevel].valide_archs |= STARPU_CUDA;
                    heteroprio->buckets[prioM2LAtLevel].factor_base_arch_index = FSTARPU_CUDA_IDX;
#ifdef STARPU_USE_CPU
                    heteroprio->buckets[prioM2LAtLevel].slow_factors_per_index[FSTARPU_CPU_IDX] = 40.0f;
#endif
                }
            }

            //prioP2MSend = 0;
            //prioP2M     = prioP2MSend+1;
            if(capacities->supportP2M(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioP2MSend;
                heteroprio->buckets[prioP2MSend].valide_archs |= STARPU_CUDA;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioP2M;
                heteroprio->buckets[prioP2M].valide_archs |= STARPU_CUDA;
            }

            //prioM2MSend = prioP2M+1;
            //prioM2M     = prioM2MSend+1;
            if(capacities->supportM2M(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioM2MSend;
                heteroprio->buckets[prioM2MSend].valide_archs |= STARPU_CUDA;

                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioM2M;
                heteroprio->buckets[prioM2M].valide_archs |= STARPU_CUDA;
            }

            // prioL2L     = prioM2L+1;
            for(int idxLevel = 2 ; idxLevel < treeHeight ; ++idxLevel){
                if(idxLevel != treeHeight-1 && capacities->supportL2L(FSTARPU_CUDA_IDX)){
                    const int prioL2LAtLevel = getPrioL2L(idxLevel);
                    heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioL2LAtLevel;
                    heteroprio->buckets[prioL2LAtLevel].valide_archs |= STARPU_CUDA;
                }
            }

            //prioL2P     = prioP2PMpi+1;
            if(capacities->supportL2P(FSTARPU_CUDA_IDX)){
                heteroprio->prio_mapping_per_arch_index[FSTARPU_CUDA_IDX][openclCountPrio++] = prioL2P;
                heteroprio->buckets[prioL2P].valide_archs |= STARPU_CUDA;
            }

            heteroprio->nb_prio_per_arch_index[FSTARPU_CUDA_IDX] = int(openclCountPrio);
        }
#endif
    }

    int getPrioP2M() const {
        return prioP2M;
    }
    int getPrioM2M(const int /*inLevel*/) const {
        return prioM2M;
    }
    int getPrioP2M(bool willBeSend) const {
        return willBeSend?prioP2MSend:prioP2M;
    }
    int getPrioM2M(const int /*inLevel*/, bool willBeSend) const {
        return willBeSend?prioM2MSend:prioM2M;
    }
    int getPrioM2L(const int inLevel) const {
        return prioM2L + (inLevel - 2)*2;
    }
    int getPrioM2LExtern(const int inLevel) const {
        return prioM2LExtern + (inLevel - 2)*2;
    }
    int getPrioL2L(const int inLevel) const {
        return prioL2L + (inLevel - 2)*2;
    }
    int getPrioL2P() const {
        return prioL2P;
    }
    int getPrioP2P() const {
        return prioP2P;
    }
    int getPrioP2PExtern() const {
        return prioP2PExtern;
    }
    int getPrioM2LMpi(const int inLevel) const {
        return prioM2LMpi + (inLevel - 2)*2;
    }
    int getPrioP2PMpi() const {
        return prioP2PMpi;
    }
};

#else // STARPU_SUPPORT_SCHEDULER

class FStarPUFmmPriorities{
    static FStarPUFmmPriorities controller;


    FStarPUFmmPriorities(){
    }

public:
    static FStarPUFmmPriorities& Controller(){
        return controller;
    }


    void init(struct starpu_conf* /*conf*/, const int /*inTreeHeight*/,
              FStarPUKernelCapacities* /*inCapacities*/){

    }

    int getPrioP2M() const {
        return 0;
    }
    int getPrioM2M(const int /*inLevel*/) const {
        return 0;
    }
    int getPrioP2M(bool willBeSend) const {
        return 0;
    }
    int getPrioM2M(const int /*inLevel*/, bool willBeSend) const {
        return 0;
    }
    int getPrioM2L(const int inLevel) const {
        return 0;
    }
    int getPrioM2LExtern(const int inLevel) const {
        return 0;
    }
    int getPrioL2L(const int inLevel) const {
        return 0;
    }
    int getPrioL2P() const {
        return 0;
    }
    int getPrioP2P() const {
        return 0;
    }
    int getPrioP2PExtern() const {
        return 0;
    }
    int getPrioM2LMpi(const int inLevel) const {
        return 0;
    }
    int getPrioP2PMpi() const {
        return 0;
    }
};

#endif // STARPU_SUPPORT_SCHEDULER

FStarPUFmmPriorities FStarPUFmmPriorities::controller;


#endif // FSTARPUFMMPRIORITIES_HPP

