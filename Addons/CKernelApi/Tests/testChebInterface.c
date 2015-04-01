#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//For timing monitoring
#include <time.h>
#include <sys/time.h>

#include "../Src/CScalfmmApi.h"

#include "../../Src/Kernels/Chebyshev/FChebInterface.h"


/**
 * This file is an example of the user defined kernel API, with link
 * to our Chebyshev Kernel
 **/

/**
 * @brief Wrapper to init internal ChebCell
 */
void* cheb_init_cell(int level, long long morton_index, int* tree_position, double* spatial_position){
    return ChebCellStruct_create(morton_index,tree_position);
}

/**
 * @brief Wrapper to free internal ChebCell
 */
void cheb_free_cell(void * inCell){
    ChebCellStruct_free(inCell);
}

/**
 * @brief Wrapper to FMM operators (refer to CScalfmmApi.h to get the
 * detailed descriptions)
 */
void cheb_p2m(void* cellData, FSize nbParticlesInLeaf, const FSize* particleIndexes, void* userData){
    ChebKernel_P2M(cellData,nbParticlesInLeaf,particleIndexes,userData);
}
void cheb_m2m(int level, void* parentCell, int childPosition, void* childCell, void* userData){
    ChebKernel_M2M(level,parentCell,childPosition,childCell,userData);
}
void cheb_m2l_full(int level, void* targetCell, void* sourceCell[343], void* userData){
    ChebKernel_M2L(level, targetCell, sourceCell, userData);
}
void cheb_l2l(int level, void* parentCell, int childPosition, void* childCell, void* userData){
    ChebKernel_L2L( level, parentCell, childPosition, childCell,  userData);
}
void cheb_l2p(void* leafCell, FSize nbParticles, const FSize* particleIndexes, void* userData){
    ChebKernel_L2P( leafCell, nbParticles, particleIndexes, userData);
}
void cheb_p2pFull(FSize nbParticles, const FSize* particleIndexes,
                  const FSize * sourceParticleIndexes[27], FSize sourceNbPart[27],void* userData) {
    ChebKernel_P2P(nbParticles, particleIndexes, sourceParticleIndexes, sourceNbPart, userData);
}


/**
 * @brief Wrapper on timeval struct
 */
typedef struct timer{
    struct timeval start;
    struct timeval end;
}Timer;

/**
 * @brief monitoring function (init start)
 */
void tic(Timer* inTimer){
    gettimeofday(&inTimer->start, NULL);
}

/**
 * @brief monitoring function (init end)
 */
void tac(Timer* inTimer){
    gettimeofday(&inTimer->end, NULL);
}

/**
 * @brief monitoring function (return elapsed time in micro second)
 */
long int get_elapsed(Timer* inTimer){
    return ((inTimer->end.tv_sec * 1000000 + inTimer->end.tv_usec)
            - (inTimer->start.tv_sec * 1000000 + inTimer->start.tv_usec));
}

/**
 * @brief monitoring function (print elapsed)
 */
void print_elapsed(Timer* inTimer){
    long int elapsed = get_elapsed(inTimer);
    printf("Elapsed : %ld us (or %f seconds)\n",elapsed,elapsed/1000000.0);
}

/**
 * @brief monitoring function (print difference :: First-Second)
 */
void print_difference_elapsed(Timer* inTimer1,Timer*inTimer2){
    long int diff = get_elapsed(inTimer1)-get_elapsed(inTimer2);
    printf("Timer Difference : %ld us (or %f seconds)\n",diff,(double)diff/1000000.0);
}


/**
 * @brief Do everything
 * @param number of particle (no default value)
 */
int main(int argc, char ** av){
    //omp_set_num_threads(1);
    printf("Start\n");
    if(argc<2){
        printf("Use : %s nb_part (optionnal : TreeHeight) \nexiting\n",av[0]);
        exit(0);
    }
    int nbPart= atoi(av[1]);
    int treeHeight = 5 ;
    if(argc>2){
        int treeHeight = atoi(av[2]);
    }


    double boxWidth = 1.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;

    //Allocation of the positions and physical values
    double* particleXYZ = malloc(sizeof(double)*3*nbPart);
    double* physicalValues = malloc(sizeof(double)*nbPart);


    {
        printf("Creating Particles:\n");
        FSize idxPart;
        for(idxPart = 0 ; idxPart < nbPart ; ++idxPart){
            particleXYZ[idxPart*3]   = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
            particleXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
            particleXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
            physicalValues[idxPart] = 1.0;

        }

    }

    {//This part will write generated particles to a file at ScalFMM
     //format in order to verify numercal results
        /* FILE * fd = fopen("input.fma","w"); */
        /* fprintf(fd,"8\t 4\n %d\n %f\t %f\t %f\t %f\n",nbPart,boxWidth/2.0, boxCenter[0],boxCenter[1],boxCenter[2]); */
        /* FSize idxPart; */
        /* for(idxPart=0 ; idxPart<nbPart ; ++idxPart){ */
        /*     fprintf(fd,"%e\t %e\t %e\t %e \n", */
        /*             particleXYZ[idxPart*3], */
        /*             particleXYZ[idxPart*3+1], */
        /*             particleXYZ[idxPart*3+2], */
        /*             physicalValues[idxPart]); */
        /* } */
        /* fclose(fd); */
    }

    scalfmm_handle handle = scalfmm_init(user_defined_kernel);

    //For Reference
    scalfmm_handle handle_ref = scalfmm_init(chebyshev);

    //Struct for user defined kernel
    struct User_Scalfmm_Cell_Descriptor cellDescriptor;
    cellDescriptor.user_init_cell = cheb_init_cell;
    cellDescriptor.user_free_cell = cheb_free_cell;

    //Struct for ref cheb kernel
    struct User_Scalfmm_Cell_Descriptor user_descr;
    user_descr.user_init_cell = NULL;
    user_descr.user_free_cell = NULL;


    // Init tree and cell
    printf("Building the tree and Initizalizing the cells:\n");

    scalfmm_build_tree(handle,treeHeight, boxWidth, boxCenter, cellDescriptor);
    scalfmm_build_tree(handle_ref,treeHeight, boxWidth, boxCenter, user_descr);

    // Insert particles
    printf("Inserting particles...\n");
    scalfmm_tree_insert_particles_xyz(handle, nbPart, particleXYZ);
    scalfmm_tree_insert_particles_xyz(handle_ref, nbPart, particleXYZ);

    //Set physical values for Cheb_ref
    scalfmm_set_physical_values(handle_ref,nbPart,physicalValues);
    //Set our callbacks
    struct User_Scalfmm_Kernel_Descriptor kernel;
    kernel.p2m = cheb_p2m;
    kernel.m2m = cheb_m2m;
    //init the other to NULL
    kernel.m2l = NULL;
    kernel.m2l_full = cheb_m2l_full;
    kernel.l2l = cheb_l2l;
    kernel.l2p = cheb_l2p;
    kernel.p2p_full = cheb_p2pFull;
    kernel.p2pinner = NULL;
    kernel.p2p = NULL;

    //Set my datas
    UserData userDatas;

    int nb_threads = omp_get_max_threads();
    int idThreads= 0;

    //Create as many kernels as there are threads in order to void concurrent write
    userDatas.kernelStruct = ChebKernelStruct_create(treeHeight,boxWidth,boxCenter);
    /* malloc(sizeof(void *)*nb_threads); */
    /* for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){ */
    /*     userDatas.kernelStruct[idThreads] = ChebKernelStruct_create(treeHeight,boxWidth,boxCenter); // Set my kernel inside userDatas */
    /* } */

    //Only read, so no split needed
    userDatas.insertedPositions = particleXYZ;                                       // Set the position
    userDatas.myPhyValues = physicalValues;                                          // Set the physical values

    //Create as many array of forces as there are threads in order to void concurrent write
    double ** forcesToStore = malloc(sizeof(double*)*nb_threads);
    //For each array, initialisation
    for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        forcesToStore[idThreads] =  malloc(sizeof(double)*nbPart*3); //allocate memory
        memset(forcesToStore[idThreads],0,sizeof(double)*nbPart*3);  //set to zero (IMPORTANT, as operators usually "+=" on forces)
    }
    userDatas.forcesComputed = forcesToStore;

    //Give ScalFMM the datas before calling fmm (this will set as well the kernel)
    scalfmm_user_kernel_config(handle,kernel,&userDatas);

    //Set timers
    Timer interface_timer,ref_timer;

    //Execute FMM
    tic(&interface_timer);
    scalfmm_execute_fmm(handle/*, kernel, &my_data*/);
    tac(&interface_timer);

    //Reduction on forces array
    {
        FSize idxPart;
        for(idThreads=1 ; idThreads<nb_threads ; ++idThreads){
            for(idxPart=0 ; idxPart<nbPart*3 ; ++idxPart){
                //Everything is stored in first array
                forcesToStore[0][idxPart] += forcesToStore[idThreads][idxPart];
            }
        }
    }

    printf("User defined Chebyshev done\n");
    print_elapsed(&interface_timer);

    tic(&ref_timer);
    scalfmm_execute_fmm(handle_ref/*, kernel, &my_data*/);
    tac(&ref_timer);

    printf("Intern Chebyshev done\n");
    print_elapsed(&ref_timer);

    //Print time results
    print_difference_elapsed(&interface_timer,&ref_timer);

    //get back the forces for ref_cheb execution
    double * forcesRef = malloc(sizeof(double)*3*nbPart);
    memset(forcesRef,0,sizeof(double)*3*nbPart);
    scalfmm_get_forces_xyz(handle_ref,nbPart,forcesRef);
    {//Comparison part
        FSize idxPart;
        int nbPartOkay = 0;
        for(idxPart=0 ; idxPart<nbPart ; ++idxPart ){
            double diffX,diffY,diffZ;
            diffX = forcesToStore[0][idxPart*3+0]-forcesRef[idxPart*3+0];
            diffY = forcesToStore[0][idxPart*3+1]-forcesRef[idxPart*3+1];
            diffZ = forcesToStore[0][idxPart*3+2]-forcesRef[idxPart*3+2];
            if(diffX < 0.00000001 && diffY < 0.00000001 && diffZ < 0.00000001){
                nbPartOkay++;
            }
            else{
                printf("id : %lld : %e, %e, %e\n",idxPart,diffX,diffY,diffZ);
            }
            //That part is to verify with our usual exec' if everything is alright
            if(idxPart == 0 || idxPart == nbPart/2 || idxPart == nbPart-1){
                printf("User one's id : %lld : %e, %e, %e\n",idxPart,
                       forcesToStore[0][idxPart*3+0],
                       forcesToStore[0][idxPart*3+1],
                       forcesToStore[0][idxPart*3+2]);
                printf("Chebyshev one's id : %lld : %e, %e, %e\n",idxPart,
                       forcesRef[idxPart*3+0],
                       forcesRef[idxPart*3+1],
                       forcesRef[idxPart*3+2]);
            }
        }
        printf("End of simulation \n \t Percentage of good parts : %d/%d (%f %%) \n",
               nbPartOkay,nbPart,(((double) nbPartOkay)/(double)nbPart)*100);
    }
    printf("Free the kernels\n");

    printf("Free the Handles ...\n");
    scalfmm_dealloc_handle(handle,cheb_free_cell);
    scalfmm_dealloc_handle(handle_ref,NULL);

    free(particleXYZ);
    free(physicalValues);
    free(forcesRef);
    //free the thread' specific datas
    for(idThreads=0 ; idThreads<nb_threads ; ++idThreads){
        free(forcesToStore[idThreads]);
    }
    ChebKernelStruct_free(userDatas.kernelStruct);

    free(userDatas.forcesComputed);
}
