#ifndef CALL_HPP
#define CALL_HPP


/** It should be compiled with C export */
extern "C" {
#include "CScalfmmApi.h"
}

#include "FInterEngine.hpp"
#include "FUserKernelEngine.hpp"

extern "C" scalfmm_handle scalfmm_init(int TreeHeight,double BoxWidth,double* BoxCenter, scalfmm_kernel_type KernelType){
    ScalFmmCoreHandle * handle = new ScalFmmCoreHandle();

    switch(KernelType){
    case 0:
        //handle->engine = new FUserKernelEngine(...);
        std::cout<< "User Kernel type unsupported yet" << std::endl;
        break;

    case 1:
        //TODO typedefs
        typedef FP2PParticleContainerIndexed<>                                 ContainerClass;
        typedef FChebCell<7>                                                         ChebCell;

        typedef FInterpMatrixKernelR                                        MatrixKernelClass;
        typedef FChebSymKernel<ChebCell,ContainerClass,MatrixKernelClass,7>        ChebKernel;

        handle->engine = new FInterEngine<ChebCell,ChebKernel>(TreeHeight,BoxWidth,BoxCenter, KernelType);
        break;
    case 2:
        //TODO typedefs
        typedef FP2PParticleContainerIndexed<>                                 ContainerClass;
        typedef FUnifCell<7>                                                         UnifCell;

        typedef FInterpMatrixKernelR                                        MatrixKernelClass;
        typedef FUnifKernel<UnifCell,ContainerClass,MatrixKernelClass,7>           UnifKernel;

        handle->engine = new FInterEngine<UnifCell,UnifKernel>(TreeHeight,BoxWidth,BoxCenter,KernelType);
        break;

    default:
        std::cout<< "Kernel type unsupported" << std::endl;
        exit(0);
        break;
    }
    return handle;
}


#endif
