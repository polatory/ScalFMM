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

// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE
// @FUSE_BLAS


#include <string>

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FRandomLoader.hpp"
#include "Containers/FOctree.hpp"

// Cell
//#include "Components/FBasicCell.hpp"
//#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Components/FTestCell.hpp"

// Particle Container
//#include "Components/FBasicParticleContainer.hpp"
#include "Components/FTestParticleContainer.hpp"

// Leaf
#include "Components/FSimpleLeaf.hpp"

// Kernel
//#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Components/FTestKernels.hpp"
#include "BalanceTree/FChebSymCostKernel.hpp"
#include "BalanceTree/FDumbCostKernel.hpp"

// Algorithm
#include "Core/FFmmAlgorithmThread.hpp"
#include "BalanceTree/FFmmAlgorithmThreadBalanced.hpp"

// Other
#include "BalanceTree/FCostZones.hpp"
#include "testFmmAlgorithmBalancedArgs.hpp"
#include "testFmmAlgorithmBalancedUtils.hpp" // last include, to shorten main file

#define ORDER 7

using FReal = double;

using CellClass         = FCostCell<FTestCell>;
using ContainerClass    = FTestParticleContainer<FReal>;
using LeafClass         = FSimpleLeaf< FReal, ContainerClass >;
using OctreeClass       = FOctree< FReal, CellClass, ContainerClass, LeafClass >;

using MatrixKernelClass = FInterpMatrixKernelR<FReal>;
// using CostKernelClass   = FChebSymCostKernel<FReal, CellClass, ContainerClass,
//                                              MatrixKernelClass, ORDER,
//                                              OctreeClass>;
using CostKernelClass   = FDumbCostKernel<FReal, CellClass, ContainerClass,
                                             OctreeClass>;

using KernelClass       = FTestKernels< CellClass, ContainerClass>;


const FReal epsilon = 1e-4;


int main(int argc, char** argv)
{
    // Handle arguments
    loadFMAAndRunFMMArgs args(argc, argv);
    omp_set_num_threads(args.zoneCount());

    /* Creating tree and insterting particles *********************************/
    FFmaGenericLoader<FReal> loader(args.inFileName().c_str());
    //FRandomLoader loader(20, 1, FPoint(0.5,0.5,0.5), 1);
    OctreeClass tree(args.treeHeight(),
                     args.subTreeHeight(),
                     loader.getBoxWidth(),
                     loader.getCenterOfBox());

    loadTree(tree, loader);
    /**************************************************************************/


    /* Compute the cost of each tree cell *************************************/
    CostKernelClass costKernel(&tree, epsilon);
    FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, CostKernelClass, LeafClass > costAlgo(&tree, &costKernel);

    costAlgo.execute();

    if (args.verboseLevel() > 1) {
        //costKernel.printResults(std::cout);
    }
    /**************************************************************************/

    std::cerr << "Running the costzones algorithm" << std::endl;
    /* Run the costzone algorithm *********************************************/
    FCostZones<OctreeClass, CellClass> costzones(&tree, args.zoneCount());
    costzones.run();

    /**************************************************************************/
    std::cerr << "Done" << std::endl;


    /* Run the balanced algorithm *********************************************/

    std::cout << "Running kernel" << std::endl;
    KernelClass computeKernel;
    FFmmAlgorithmThreadBalanced<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >
            fmmAlgo(&tree, &computeKernel, costzones.getZoneBounds(), costzones.getLeafZoneBounds());
    //FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >
    //     fmmAlgo(&tree, &computeKernel);
    
    fmmAlgo.execute();
    /**************************************************************************/


    /* Check the results ******************************************************/
    ValidateFMMAlgo<OctreeClass, CellClass, ContainerClass, LeafClass>(&tree);


    return EXIT_SUCCESS;
}


