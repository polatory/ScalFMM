#ifndef __STARPU_CODELET_PARAMETERS_H__
#define __STARPU_CODELET_PARAMETERS_H__

#include <starpu.h>
#include "../StarPUUtils/FStarPUUtils.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

static inline void p2p_cl_in_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  int i;
  starpu_codelet_unpack_args(task->cl_arg,
	  	  	  	  &wrapperptr,
				  &i,
				  &parameters[0]);
}

static inline void p2p_cl_inout_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  std::vector<OutOfBlockInteraction>* outsideInteractions;
  int i;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &outsideInteractions,
				  &i,
				  &parameters[0],
				  &parameters[1]);
}

static inline void p2m_cl_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  int idxLevel, i;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &idxLevel,
				  &i,
			     	  &parameters[0]);
}
  
static inline void m2m_cl_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  int idxLevel, i;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &idxLevel,
				  &i,
			     	  &parameters[0],
				  &parameters[1]);
}
  
static inline void m2l_cl_in_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  int idxLevel, i;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &idxLevel,
				  &i,
			     	  &parameters[0],
				  &parameters[1]);
}

static inline void m2l_cl_inout_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  const std::vector<OutOfBlockInteraction>* outsideInteractions;
  int idxLevel, i, m;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &idxLevel,
			          &outsideInteractions,
				  &i,
                                  &m,
				  &parameters[0]);
}

static inline void l2l_cl_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  int idxLevel, i;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &idxLevel,
				  &i,
			     	  &parameters[0],
		        	  &parameters[1],
				  &parameters[2]);
}

static inline void l2p_cl_perf_func(struct starpu_task *task, double *parameters){
  FStarPUPtrInterface* wrapperptr;
  int i;
  starpu_codelet_unpack_args(task->cl_arg,
      	  	  	  	  &wrapperptr,
				  &i,
			     	  &parameters[0]);
}
  
#ifdef __cplusplus
}
#endif

#endif /* __STARPU_CODELET_PARAMETERS_H__ */
