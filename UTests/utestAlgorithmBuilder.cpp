// ===================================================================================
// Copyright ScalFmm 2011 INRIA
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
#include "FUTester.hpp"

#include "../Src/Core/FAlgorithmBuilder.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Components/FBasicKernels.hpp"


/** This class test the core algorithm builder */
class TestBuilder : public FUTester<TestBuilder> {

    void Test(){
        typedef FTestCell                   CellClass;
        typedef FTestParticleContainer      ContainerClass;

        typedef FSimpleLeaf< ContainerClass >                     LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

        const int height = 5;
        const FReal dim = 1.0;
        const FPoint center(0.0,0.0,0.0);

        OctreeClass tree(height, 2, dim, center);
        KernelClass kernel;

        {
            FAlgorithmBuilder::SimulationProperties properties = FAlgorithmBuilder::BuildKernelSimulationProperties(height,center,dim,false);
            uassert(properties.centerOfBox.getX() == center.getX() && properties.centerOfBox.getY() == center.getY() &&
                    properties.centerOfBox.getZ() == center.getZ() );
            uassert(properties.dimOfBox == dim);
            uassert(properties.height == height);

            FAbstractAlgorithm*const algo = FAlgorithmBuilder::BuildAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>(&tree, &kernel, 0, false);
#ifndef ScalFMM_USE_MPI
            uassert(dynamic_cast<FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr ||
                    dynamic_cast<FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr);
#else
            uassert(dynamic_cast<FFmmAlgorithmThreadProc<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr);
#endif
            delete algo;
        }
        {
            FAlgorithmBuilder::SimulationProperties properties = FAlgorithmBuilder::BuildKernelSimulationProperties(height,center,dim,true);
            uassert(properties.dimOfBox != dim);
            uassert(properties.height != height);

            FAbstractAlgorithm*const algo = FAlgorithmBuilder::BuildAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>(&tree, &kernel, 0, true);
#ifndef ScalFMM_USE_MPI
            uassert(dynamic_cast<FFmmAlgorithmPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr );
#else
            uassert(dynamic_cast<FFmmAlgorithmThreadProcPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>*>(algo) != nullptr);
#endif

            delete algo;
        }

    }

    // set test
    void SetTests(){
        AddTest(&TestBuilder::Test,"Test Algo Creation");
    }
};

// You must do this
TestClass(TestBuilder)



