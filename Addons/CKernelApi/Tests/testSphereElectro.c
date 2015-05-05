#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>




//For timing monitoring
#include <time.h>
#include <sys/time.h>

#include "../Src/CScalfmmApi.h"

#include "../../Src/Kernels/Chebyshev/FChebInterface.h"

double getRandom(){
    return (random()/(double)(RAND_MAX));
}

void generateSurfacePointOnUnitSphere(int N, double * points){
    double u, v, theta, phi, sinPhi ;
    //
    int j = 0,i=0 ;
    for ( i = 0 ; i< N ; ++i, j+=3)  {
        //
        u = getRandom() ;  v = getRandom() ;
        theta  = 2*M_PI*u ;
        phi     = acos(2*v-1);
        sinPhi = sin(phi);
        //
        points[j]   =     cos(theta)*sinPhi ;
        points[j+1] =   sin(theta)*sinPhi ;
        points[j+2] =   2*v-1 ;
        //
    }
}

void generateSurfacePoints(double rayon, double* centre, int nbDePoints, double* points){

    generateSurfacePointOnUnitSphere(nbDePoints , points) ;
    int j =0,i=0 ;
    for ( i = 0 ; i< nbDePoints ; ++i, j+=3)  {
        points[j]    *= rayon + centre[0];
        points[j+1]  *= rayon + centre[1];
        points[j+2]  *= rayon + centre[2];
    }
}

void generateInsidePoints(double rayon, double*centre, int nbDePoints,  double* points){
    generateSurfacePointOnUnitSphere(nbDePoints, points);
    int j=0;
    double u;
    for(j=0 ; j<nbDePoints ; ++j){
        u = getRandom();
        points[j]    *= (rayon + centre[0])*u;
        points[j+1]  *= (rayon + centre[1])*u;
        points[j+2]  *= (rayon + centre[2])*u;
    }
}

void displayPoints(int nbPoints, double * points){
    int i = 0;
    for(i=0 ; i<nbPoints ; ++i){
        printf("%e %e %e \n",points[i*3],points[i*3+1],points[i*3+2]);
    }
}

void displayArray(int nbValue, double * array){
    int i = 0;
    for(i=0 ; i<nbValue ; ++i){
        printf("%e \n",array[i]);
    }
}

void getNormal(double * positions, double * normeToFill){
    int i;
    double norme = sqrt(positions[0]*positions[0] + positions[1]*positions[1] + positions[2]*positions[2]);
    for(i=0 ; i<3 ; ++i){
        normeToFill[i] = positions[i]/norme;
    }
    printf("Tgt Norme %e - %e - %e\n",
           normeToFill[0],
           normeToFill[1],
           normeToFill[2]);
}

void computeNormalXForces(int nbPoints, double * forcesToRead, double * positionsToRead, double * arrayToFill){
    double * currentNormal = malloc(sizeof(double)*3);
    int idxPart,i;
    for(idxPart = 0 ; idxPart<nbPoints ; ++idxPart){
        getNormal(&positionsToRead[idxPart],currentNormal); //get the norme
        for(i=0 ; i<3 ; ++i){
            arrayToFill[idxPart] += currentNormal[i]*forcesToRead[idxPart+i];
        }
    }
    free(currentNormal);
}

int main(int argc, char ** av){
    printf("Start\n");
    if(argc<2){
        printf("Use : %s nb_part(cible) (optionnal : TreeHeight) \nexiting\n",av[0]);
        exit(0);
    }
    int nbPartTarget= atoi(av[1]);
    int treeHeight = 5 ;
    if(argc>2){
        int treeHeight = atoi(av[2]);
    }

    double boxWidth = 2.0;
    double boxCenter[3];
    boxCenter[0] = boxCenter[1] = boxCenter[2] = 0.0;

    int i;
    //Allocation of the target points
    double * targetsXYZ = malloc(sizeof(double)* 3*nbPartTarget);
    double * targetsPhiValues = malloc(sizeof(double)* nbPartTarget);
    //Memset (au cas ou)
    memset(targetsXYZ,0,sizeof(double)*3*nbPartTarget);
    memset(targetsPhiValues,0,sizeof(double)*nbPartTarget);
    //Fill
    for(i=0 ; i<nbPartTarget ; ++i){
        targetsPhiValues[i] = -1.0;
    }
    generateSurfacePoints(1.0,boxCenter,nbPartTarget,targetsXYZ);
    printf("Surface points generated \n");

    //Allocation of the sources points
    int nbPartSource = 10;
    double * sourceXYZ = malloc(sizeof(double)* 3*nbPartSource);
    double * sourcePhiValues = malloc(sizeof(double)* nbPartSource);
    //Set to Zero
    memset(sourceXYZ,0,3*sizeof(double)*nbPartSource);
    memset(sourcePhiValues,0,sizeof(double)*nbPartSource);
    //Fill
    for(i=0 ; i<nbPartSource ; ++i){
        sourcePhiValues[i] = 1.0;
    }
    generateInsidePoints(1.0,boxCenter,nbPartSource,sourceXYZ);
    //displayPoints(nbPartTarget,targetsXYZ);

    printf("Inside points generated \n");
    //displayPoints(nbPartSource,sourceXYZ);
    //Creation of arrays to store forces
    double * arrayOfForces = malloc(sizeof(double )* 3 * (nbPartSource+nbPartTarget));
    memset(arrayOfForces,0,sizeof(double)* 3 * (nbPartTarget));

    {//Start of computation

        //For handling the library
        scalfmm_handle handle = scalfmm_init(chebyshev,source_target);

        //Struct for ref cheb kernel
        struct User_Scalfmm_Cell_Descriptor user_descr;
        user_descr.user_init_cell = NULL;
        user_descr.user_free_cell = NULL;
        //Set algorithm to source target
        //scalfmm_algorithm_config(handle,source_target);
        //Build the tree
        scalfmm_build_tree(handle,treeHeight, boxWidth, boxCenter, user_descr);

        //Insert Sources and targets
        scalfmm_tree_insert_particles_xyz(handle,nbPartSource,sourceXYZ,SOURCE);
        printf("Sources inserted \n");
        scalfmm_tree_insert_particles_xyz(handle,nbPartTarget,targetsXYZ,TARGET);
        printf("Targets inserted \n");
        //Since we inserted first the sources, then sources will get
        //indices from 0 to (nbPartSource-1), and targets from
        //(nbPartSource) to nbPartSource+nbPartTarget-1).

        int * arrayofIndicesSource = malloc(sizeof(int)*nbPartSource);
        int * arrayofIndicesTarget = malloc(sizeof(int)*nbPartTarget);
        {//Set physical values

            //SRC
            int idPart;
            for(idPart = 0 ; idPart<nbPartSource ; ++idPart){
                arrayofIndicesSource[idPart] = idPart;
            }
            scalfmm_set_physical_values_npart(handle,nbPartSource,arrayofIndicesSource,sourcePhiValues,SOURCE);
            //TGT
            for(idPart = 0 ; idPart<nbPartTarget ; ++idPart){
                arrayofIndicesTarget[idPart] = idPart; // here, we add the number of sources previously inserted
            }
            scalfmm_set_physical_values_npart(handle,nbPartTarget,arrayofIndicesTarget,targetsPhiValues,TARGET);



        }
        //Computation
        scalfmm_execute_fmm(handle/*, kernel, &my_data*/);

        //Get back the forces
        scalfmm_get_forces_xyz(handle,nbPartTarget,arrayOfForces,TARGET);
        scalfmm_get_forces_xyz(handle,nbPartSource,&arrayOfForces[nbPartTarget],SOURCE);
        printf("Forces computed : \n");
        displayPoints(nbPartTarget+nbPartSource,arrayOfForces);
        printf("As expected, Source forces are 0\n \n");
        //Release memory used :
        free(arrayofIndicesSource);
        free(arrayofIndicesTarget);

        scalfmm_dealloc_handle(handle,NULL);

    }

    {//Let's check the result, we computed fr each target part its forces
        //Storage of reference forces
        double * arrayRefForces = malloc(sizeof(double)*nbPartTarget*3);
        memset(arrayRefForces,0,sizeof(double)*nbPartTarget*3);

        int idTgt;
        for(idTgt = 0 ; idTgt<nbPartTarget ; ++idTgt){
            int idSrc;
            double dx,dy,dz;
            for(idSrc = 0 ; idSrc<nbPartTarget ; ++idSrc){
                //First compute dist.
                dx = sourceXYZ[idSrc+0] - targetsXYZ[idTgt+0];
                dy = sourceXYZ[idSrc+1] - targetsXYZ[idTgt+1];
                dz = sourceXYZ[idSrc+2] - targetsXYZ[idTgt+2];

                //Secondly, compute coeff
                double coeffs = targetsPhiValues[idTgt] * sourcePhiValues[idSrc];
                double one_over_r = 1.0/(sqrt(dx*dx+dy*dy+dz*dz));
                double one_over_r3 = one_over_r * one_over_r * one_over_r;

                arrayRefForces[idTgt*3+0] += dx*coeffs*one_over_r3;
                arrayRefForces[idTgt*3+1] += dy*coeffs*one_over_r3;
                arrayRefForces[idTgt*3+2] += dz*coeffs*one_over_r3;

            }
        }

        {//Then, we compare
            double errorCumul = 0;
            int idArr;
            for(idArr = 0 ; idArr<nbPartTarget ; ++idArr){
                errorCumul += fabs(arrayRefForces[idArr+0]-arrayOfForces[idArr+0]);
                errorCumul += fabs(arrayRefForces[idArr+1]-arrayOfForces[idArr+1]);
                errorCumul += fabs(arrayRefForces[idArr+2]-arrayOfForces[idArr+2]);
                printf("Directly Computed %e %e %e\n",
                       arrayRefForces[idArr+0],
                       arrayRefForces[idArr+1],
                       arrayRefForces[idArr+2]);
            }
            printf("Error cumul : %e\n",errorCumul);
        }
    }


    //Part where we apply normal on target's forces vector
    //Copying each target's parts forces,
    double * targetsForces = malloc(sizeof(double) * 3 * nbPartTarget);
    memcpy(targetsForces,arrayOfForces,sizeof(double)*3*nbPartTarget);

    double * normeXForces =  malloc(sizeof(double) * nbPartTarget);
    memset(normeXForces,0,sizeof(double) * nbPartTarget);

    computeNormalXForces(nbPartTarget,targetsForces,targetsXYZ,normeXForces);
    printf("For each target, we display [Normal To Sphere] . [Force product] \n");
    displayArray(nbPartTarget,normeXForces);


    //Free memory
    free(sourceXYZ);
    free(sourcePhiValues);
    free(targetsXYZ);
    free(targetsPhiValues);
    free(arrayOfForces);
    return EXIT_SUCCESS;
}
