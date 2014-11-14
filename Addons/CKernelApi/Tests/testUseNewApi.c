#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../Src/CScalfmmApi.h"



void compute_displacement_from_forces(double * fXYZ, double* dXYZ, int nb_xyz){
    //Do some things
}



int main(int argc, char ** av){

    scalfmm_kernel_type myChoice = chebyshev;

    //Octree configuration
    int TreeHeight = 5;
    double boxWidth = 1.0;
    double boxCenter[3] = {0.0,0.0,0.0};

    //Init our lib
    scalfmm_handle handle = scalfmm_init(TreeHeight,boxWidth,boxCenter,myChoice); //The tree is built

    //Creation of an array of particles
    int nb_of_parts = 200;
    int idxPart;
    double * positionsXYZ = malloc(sizeof(double)*3*nb_of_parts);
    for(idxPart = 0; idxPart<nb_of_parts ; ++idxPart){
        positionsXYZ[idxPart*3+0] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[0];
        positionsXYZ[idxPart*3+1] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[1];
        positionsXYZ[idxPart*3+2] = (random()/(double)(RAND_MAX))*boxWidth - boxWidth/2 + boxCenter[2];
    }

    //Creation of charge for each part
    double * array_of_charge = malloc(sizeof(double)*nb_of_parts);
    for(idxPart = 0; idxPart<nb_of_parts ; ++idxPart){
        array_of_charge[idxPart] = (random()/(double)(RAND_MAX)); //charge in [-1,1]
    }

    //Inserting the array in the tree
    scalfmm_tree_insert_particles_xyz(handle,nb_of_parts,positionsXYZ);
    //Set the charge
    scalfmm_set_physical_values(handle,nb_of_parts,array_of_charge);


    //Computation Part
    int nb_iteration = 100;
    int curr_iteration = 0;

    //Array to store the forces
    double * array_of_forces = malloc(sizeof(double)*3*nb_of_parts);

    //Array to store the displacement
    double * array_of_displacement = malloc(sizeof(double)*3*nb_of_parts);

    //Start of an iteration loop
    while(curr_iteration < nb_iteration){
        //Execute
        scalfmm_execute_fmm(handle);

        //Get the resulting forces
        scalfmm_get_forces_xyz(handle,nb_of_parts,array_of_forces);

        //User function to compute the movement of each part
        compute_displacement_from_forces(array_of_forces,array_of_displacement,nb_of_parts);

        //Apply displacement computed
        scalfmm_add_to_positions_xyz(handle,nb_of_parts,array_of_displacement);

        curr_iteration++;
    }

    //End of Computation, useer get the position after a hundreds of iterations

    //Free memory
    free(positionsXYZ);
    free(array_of_charge);
    free(array_of_forces);
    free(array_of_displacement);

    return 0;
}
