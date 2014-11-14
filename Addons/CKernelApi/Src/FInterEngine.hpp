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
 * @file This file contains a class that inherits from FScalFMMEngine,
 * and will implement the API functions for Interpolations kernels.
 */

#ifndef FINTERENGINE_HPP
#define FINTERENGINE_HPP

#include "FScalFMMEngine.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif


/**
 * @class FInterEngine implements API for Interpolations kernels, its
 * templates can be ChebCell/ChebKernel or UnifCell/UnifKernel
 */
template<class InterCell,class InterKernel,
         class ContainerClass = FP2PParticleContainerIndexed<>,
         class LeafClass = FSimpleLeaf<FP2PParticleContainerIndexed<>>,
         class MatrixKernelClass = FInterpMatrixKernelR>
class FInterEngine : public FScalFMMEngine{
private:
    //Pointer to the kernel to be executed
    InterKernel * kernel;

    //Link to the tree
        FOctree<InterCell,ContainerClass,LeafClass> * octree;
public:
    /**
     * @brief Constructor : its only build the tree
     * @param TreeHeight Height of the tree
     * @param BoxWidth box Width
     * @param BoxCenter double[3] coordinate of the center of the
     * simulation box
     */
    FInterEngine(int TreeHeight, double BoxWidth , double * BoxCenter) : kernel(nullptr), octree(nullptr){
        octree = new FOctree<InterCell,ContainerClass,LeafClass>(TreeHeight,FMath::Min(3,TreeHeight-1),BoxWidth,FPoint(BoxCenter));
    }

    //TODO free kernel too
    ~FInterEngine(){
        free(octree);
    }
};


#endif
