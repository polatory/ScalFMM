// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
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

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include "FUKernelTester.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../Src/Adaptive/FAdaptiveCell.hpp"
#include "../Src/Adaptive/FAdaptiveKernelWrapper.hpp"
#include "../Src/Adaptive/FAbstractAdaptiveKernel.hpp"
#include "../Src/Adaptive/FAdaptChebSymKernel.hpp"

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
        const unsigned int ORDER = 6;
        typedef FChebCell<ORDER>                                        CellClass;
        typedef FP2PParticleContainerIndexed<>            ContainerClass;
        typedef FSimpleLeaf<ContainerClass>    LeafClass;
        typedef FInterpMatrixKernelR                               MatrixKernelClass;
        typedef FAdaptiveChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FAdaptiveCell< CellClass, ContainerClass >                                        CellWrapperClass;
        typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >   KernelWrapperClass;
        typedef FOctree< CellWrapperClass, ContainerClass , LeafClass >                  OctreeClass;
        // FFmmAlgorithmTask FFmmAlgorithmThread
        typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

        // run test
        RunTest<CellClass,ContainerClass,KernelWrapperClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(
                          [&](int NbLevels, FReal boxWidth, FPoint centerOfBox, const MatrixKernelClass *const MatrixKernel){
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





