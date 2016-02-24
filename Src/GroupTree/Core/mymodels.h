#ifndef __MYMODELS_H__
#define __MYMODELS_H__

#include <starpu.h>

#ifdef __cplusplus
extern "C"
{
#endif

static inline double M2M_cost_function(struct starpu_task *task, unsigned nimpl){
        printf("\nInside M2M cost function");
    	return 0.664437*1000; //Time is in milliseconds
    }



#ifdef __cplusplus
}
#endif

#endif /* __MYMODELS_H__ */
