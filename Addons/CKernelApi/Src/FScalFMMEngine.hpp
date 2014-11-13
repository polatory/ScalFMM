// ===================================================================================
// Copyright ScalFmm 2014 I
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


/**
 * @file This file contain a class, gathering all the function that
 * can be called by the ScalFMM API. Documentation for each function
 * can be found in the C Header.
 *
 */

#include "../../../Src/Utils/FAssert.hpp"


struct ScalFmmCoreHandle {
    struct ScalFmmCoreConfig {
        // Read/Write parameter
        int treeHeight;     //  Number of level in the octree
        FReal boxWidth;    // Simulation box size (root level)
        FPoint boxCenter; // Center position of the box simulation(FReal[3])
    };

    ScalFmmCoreConfig config;
    FScalFMMEngine* engine;
};



/**
 * @class FScalFMMEngine
 */
class FScalFMMEngine{
private:
    scalfmm_kernel_type kernel;
    FVector<bool> progress;

public:
    scalfmm_kernel_type getKernelType(){
        return this->kernel;
    }

    //Function about the tree (to be redefined)

    virtual void tree_insert_particles( int NbPositions, double * arrayX, double * arrayY, double * arrayZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }

    virtual void tree_insert_particles_xyz( int NbPositions, double * XYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }

    virtual void set_physical_values( int nbPhysicalValues, double * physicalValues){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }

    virtual void get_physical_values( int nbPhysicalValues, double * physicalValues){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }

    virtual void set_physical_values_npart( int nbPhysicalValues,
                                            int* idxOfParticles, double * physicalValues){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_physical_values_npart( int nbPhysicalValues,
                                            int* idxOfParticles, double * physicalValues){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }

    //To get the result
    virtual void get_forces( int nbParts, double * forcesToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_forces_npart( int nbParts, int* idxOfParticles, double * forcesToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_forces_xyz( int nbParts, double * fX, double* fY, double* fZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_forces_xyz_npart( int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    //To set initial condition
    virtual void set_forces( int nbParts, double * forcesToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_forces_npart( int nbParts, int* idxOfParticles, double * forcesToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_forces_xyz( int nbParts, double * fX, double* fY, double* fZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_forces_xyz_npart( int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    //To deal with potential
    virtual void get_potentials( int nbParts, double * potentialsToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_potentials( int nbParts, double * potentialsToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_potentials_npart( int nbParts, int* idxOfParticles, double * potentialsToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_potentials_npart( int nbParts, int* idxOfParticles, double * potentialsToFill){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }

    //To deal with positions
    virtual void add_to_positions( int NbPositions, double * updatedXYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void add_to_positions_xyz( int NbPositions, double * X, double * Y , double * Z){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void add_to_positions_npart( int NbPositions, int* idxOfParticles, double * updatedXYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void add_to_positions_xyz_npart( int NbPositions, int* idxOfParticles,
                                             double * X, double * Y , double * Z){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions( int NbPositions, double * updatedXYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions_xyz( int NbPositions, double * X, double * Y , double * Z){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions_npart( int NbPositions, int* idxOfParticles, double * updatedXYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void set_positions_xyz_npart( int NbPositions, int* idxOfParticles,
                                          double * X, double * Y , double * Z){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_positions( int NbPositions, double * updatedXYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_positions_xyz( int NbPositions, double * X, double * Y , double * Z){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_positions_npart( int NbPositions, int* idxOfParticles, double * updatedXYZ){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }
    virtual void get_positions_xyz_npart( int NbPositions, int* idxOfParticles,
                                          double * X, double * Y , double * Z){
        FAssertLF(1,"No tree instancied, exiting ...\n");
    }


    //User define Kernel Part
    virtual void user_kernel_config( Scalfmm_Kernel_Descriptor userKernel, void * userDatas){
        FAssertLF(1,"No user kernel defined, exiting ...\n");
    }

    virtual void scalfmm_execute_fmm(){
        FAssertLF(1,"No kernel set, exiting ...\n");
    }

};


extern "C" scalfmm_init(int TreeHeight,double BoxWidth,double* BoxCenter, scalfmm_kernel_type KernelType){
    ScalFmmCoreHandle * handle = new ScalFmmCoreHandle();

    switch(KernelType){
    case 0:
        //handle->engine = new FUserKernelEngine(...);
        break;
    case 1:
        //handle->engine = new FChebEngine(...);
        break;
    case 2:
        //handle->engine = new FLagrangeEngine(...);
        break;
    default:
        cout<< "Kernel type unsupported" << std::endl;
        exit(0);
        break;
    }
}

extern "C" void scalfmm_tree_insert_particles(scalfmm_handle Handle, int NbPositions, double * arrayX, double * arrayY, double * arrayZ){
    ((ScalFmmCoreHandle *) Handle)->engine->tree_insert_particles(NbPositions, arrayX, arrayY, arrayZ);
}

extern "C" void scalfmm_tree_insert_particles_xyz(scalfmm_handle Handle, int NbPositions, double * XYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->tree_insert_particles_xyz(NbPositions, XYZ);
}

extern "C" void scalfmm_set_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_physical_values(nbPhysicalValues, physicalValues);
}

extern "C" void scalfmm_get_physical_values(scalfmm_handle Handle, int nbPhysicalValues, double * physicalValues){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_physical_values(nbPhysicalValues, physicalValues);
}

extern "C" void scalfmm_set_physical_values_npart(scalfmm_handle Handle, int nbPhysicalValues,
                                                  int* idxOfParticles, double * physicalValues){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_physical_values_npart(nbPhysicalValues,
                                                                       idxOfParticles, physicalValues);
}
extern "C" void scalfmm_get_physical_values_npart(scalfmm_handle Handle, int nbPhysicalValues,
                                                  int* idxOfParticles, double * physicalValues){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_physical_values_npart(nbPhysicalValues,
                                                                       idxOfParticles, physicalValues);
}

//To get the result
extern "C" void scalfmm_get_forces(scalfmm_handle Handle, int nbParts, double * forcesToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_forces(nbParts, forcesToFill);
}

extern "C" void scalfmm_get_forces_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * forcesToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_forces_npart(nbParts, idxOfParticles, forcesToFill);
}
extern "C" void scalfmm_get_forces_xyz(scalfmm_handle Handle, int nbParts, double * fX, double* fY, double* fZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_forces_xyz(nbParts,fX, fY, fZ) ;
}

extern "C" void scalfmm_get_forces_xyz_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_forces_xyz_npart(nbParts, idxOfParticles, fX, fY, fZ) ;
}

//To set initial condition
extern "C" void scalfmm_set_forces(scalfmm_handle Handle, int nbParts, double * forcesToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_forces(nbParts, forcesToFill) ;
}

extern "C" void scalfmm_set_forces_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * forcesToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_forces_npart(nbParts, idxOfParticles, forcesToFill) ;
}

extern "C" void scalfmm_set_forces_xyz(scalfmm_handle Handle, int nbParts, double * fX, double* fY, double* fZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_forces_xyz(nbParts, fX, fY, fZ) ;
}

extern "C" void scalfmm_set_forces_xyz_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * fX, double* fY, double* fZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_forces_xyz_npart(nbParts, idxOfParticles, fX, fY, fZ) ;
}

//To deal with potential
extern "C" void scalfmm_get_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_potentials(nbParts, potentialsToFill) ;
}

extern "C" void scalfmm_set_potentials(scalfmm_handle Handle, int nbParts, double * potentialsToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_potentials(nbParts, potentialsToFill) ;
}

extern "C" void scalfmm_get_potentials_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * potentialsToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_potentials_npart(nbParts, idxOfParticles, potentialsToFill) ;
}

extern "C" void scalfmm_set_potentials_npart(scalfmm_handle Handle, int nbParts, int* idxOfParticles, double * potentialsToFill){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_potentials_npart(nbParts, idxOfParticles, potentialsToFill) ;
}


//To deal with positions
extern "C" void scalfmm_add_to_positions(scalfmm_handle Handle, int NbPositions, double * updatedXYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->add_to_positions(NbPositions, updatedXYZ);
}

extern "C" void scalfmm_add_to_positions_xyz(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z){
    ((ScalFmmCoreHandle * ) Handle)->engine->add_to_positions_xyz(NbPositions, X, Y , Z);
}

extern "C" void scalfmm_add_to_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->add_to_positions_npart(NbPositions, idxOfParticles, updatedXYZ);
}

extern "C" void scalfmm_add_to_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles,
                                                   double * X, double * Y , double * Z){
    ((ScalFmmCoreHandle * ) Handle)->engine->add_to_positions_xyz_npart(NbPositions, idxOfParticles,
                                                                        X, Y , Z);
}

extern "C" void scalfmm_set_positions(scalfmm_handle Handle, int NbPositions, double * updatedXYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_positions(NbPositions, updatedXYZ);
}

extern "C" void scalfmm_set_positions_xyz(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_positions_xyz(NbPositions, X, Y , Z);
}

extern "C" void scalfmm_set_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_positions_npart(NbPositions, idxOfParticles, updatedXYZ);
}

extern "C" void scalfmm_set_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles,
                                                double * X, double * Y , double * Z){
    ((ScalFmmCoreHandle * ) Handle)->engine->set_positions_xyz_npart(NbPositions, idxOfParticles,
                                                                     X, Y , Z);
}

extern "C" void scalfmm_get_positions(scalfmm_handle Handle, int NbPositions, double * updatedXYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_positions(NbPositions, updatedXYZ);
}

extern "C" void scalfmm_get_positions_xyz(scalfmm_handle Handle, int NbPositions, double * X, double * Y , double * Z){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_positions_xyz(NbPositions, X, Y , Z);
}

extern "C" void scalfmm_get_positions_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles, double * updatedXYZ){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_positions_npart(NbPositions, idxOfParticles, updatedXYZ);
}

extern "C" void scalfmm_get_positions_xyz_npart(scalfmm_handle Handle, int NbPositions, int* idxOfParticles,
                                                double * X, double * Y , double * Z){
    ((ScalFmmCoreHandle * ) Handle)->engine->get_positions_xyz_npart(NbPositions, idxOfParticles,
                                                                     X, Y , Z);
}
extern "C" void scalfmm_execute_fmm(scalfmm_handle Handle){
    ((ScalFmmCoreHandle * ) Handle)->engine->scalfmm_execute_fmm();
}
