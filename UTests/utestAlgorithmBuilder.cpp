// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#include "FUTester.hpp"

#include "Core/FAlgorithmBuilder.hpp"

#include "Containers/FOctree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"
#include "Components/FBasicKernels.hpp"


/** This class test the core algorithm builder */
#ifndef SCALFMM_USE_MPI
class TestBuilder : public FUTester<TestBuilder> {
#else
class TestBuilder : public FUTesterMpi<TestBuilder> {
#endif
    void Test(){
        typedef double FReal;
        typedef FTestCell                   CellClass;
        typedef FTestParticleContainer<FReal>      ContainerClass;

        typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

        const int height = 5;
        const FReal dim = 1.0;
        const FPoint<FReal> center(0.0,0.0,0.0);

        OctreeClass tree(height, 2, dim, center);
        KernelClass kernel;

        {
            FAlgorithmBuilder<FReal>::SimulationProperties properties = FAlgorithmBuilder<FReal>::BuildKernelSimulationProperties(height,center,dim,false);
            uassert(properties.centerOfBox.getX() == center.getX() && properties.centerOfBox.getY() == center.getY() &&
                    properties.centerOfBox.getZ() == center.getZ() );
            uassert(properties.dimOfBox == dim);
            uassert(properties.height == height);

            FAbstractAlgorithm*const algo = FAlgorithmBuilder<FReal>::BuildAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>(&tree, &kernel,
#ifndef SCALFMM_USE_MPI
                                                                                                                                                     0,
#else
                                                                                                                                                     MPI_COMM_WORLD,
#endif
                                                                                                                                                     false);
#ifndef SCALFMM_USE_MPI
            uassert(dynamic_cast<FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr ||
                    dynamic_cast<FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr);
#else
            uassert(dynamic_cast<FFmmAlgorithmThreadProc<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr);
#endif
            delete algo;
        }
        {
            FAlgorithmBuilder<FReal>::SimulationProperties properties = FAlgorithmBuilder<FReal>::BuildKernelSimulationProperties(height,center,dim,true);
            uassert(properties.dimOfBox != dim);
            uassert(properties.height      != height);

            FAbstractAlgorithm*const algo = FAlgorithmBuilder<FReal>::BuildAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>(&tree, &kernel,
 #ifndef SCALFMM_USE_MPI
                                                                                                                                                      0,
 #else
                                                                                                                                                      MPI_COMM_WORLD,
 #endif
                                                                                                                                                      true);
#ifndef SCALFMM_USE_MPI
            uassert(dynamic_cast<FFmmAlgorithmPeriodic<FReal,OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr );
#else
            uassert(dynamic_cast<FFmmAlgorithmThreadProcPeriodic<FReal, OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr);
#endif

            delete algo;
        }

    }

    // set test
    void SetTests(){
        AddTest(&TestBuilder::Test,"Test Algo Creation");
    }

#ifdef SCALFMM_USE_MPI
public:
    TestBuilder(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }
#endif
};

// You must do this
#ifndef SCALFMM_USE_MPI
TestClass(TestBuilder)
#else
TestClassMpi(TestBuilder)
#endif


