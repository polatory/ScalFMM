// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================

// ==== CMAKE =====
// @FUSE_BLAS
// @SCALFMM_PRIVATE
// ================

#include "FUKernelTester.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Adaptive/FAdaptiveCell.hpp"
#include "Adaptive/FAdaptiveKernelWrapper.hpp"
#include "Adaptive/FAbstractAdaptiveKernel.hpp"
#include "Adaptive/FAdaptChebSymKernel.hpp"

/*
  In this test we compare the Chebyschev fmm results and the direct results.
 */

/** the test class
 *
 */
class TestChebyshevDirect : public FUKernelTester<TestChebyshevDirect> {

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////


    /** TestChebKernel */
    void TestChebKernel(){
        typedef double FReal;
        const unsigned int ORDER = 6;
        typedef FChebCell<FReal,ORDER>                                        CellClass;
        typedef FP2PParticleContainerIndexed<FReal>            ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass>    LeafClass;
        typedef FInterpMatrixKernelR<FReal>                               MatrixKernelClass;
        typedef FAdaptiveChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FAdaptiveCell< CellClass, ContainerClass >                                        CellWrapperClass;
        typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >   KernelWrapperClass;
        typedef FOctree< FReal, CellWrapperClass, ContainerClass , LeafClass >                  OctreeClass;
        // FFmmAlgorithmTask FFmmAlgorithmThread
        typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelWrapperClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(
                          [&](int NbLevels, FReal boxWidth, FPoint<FReal> centerOfBox, const MatrixKernelClass *const MatrixKernel){
             return std::unique_ptr<KernelWrapperClass>(new KernelWrapperClass(NbLevels, boxWidth, centerOfBox, MatrixKernel, 4, 4));
        });
    }



    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestChebyshevDirect::TestChebKernel,"Test Chebyshev Kernel with one big SVD");
    }
};


// You must do this
TestClass(TestChebyshevDirect)





