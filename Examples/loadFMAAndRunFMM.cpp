// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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

#include <fstream>
#include <memory>
#include <string>


//#include "Utils/FMath.hpp"
//#include "Utils/FParameters.hpp"
//#include "Utils/FParameterNames.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Core/FFmmAlgorithm.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Containers/FOctree.hpp"
#include "Components/FBasicCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Components/FBasicParticleContainer.hpp"

#include "FFmmAlgorithmThreadBalanced.hpp"
#include "loadFMAAndRunFMMArgs.hpp"
#include "loadFMAAndRunFMMUtils.hpp"
#include "FChebBalanceSymKernel.hpp"
#include "CostZones.hpp"

using CellClass         = FCostCell;
using ContainerClass    = FBasicParticleContainer<0>;
using LeafClass         = FSimpleLeaf< ContainerClass >;
using OctreeClass       = FOctree< CellClass, ContainerClass, LeafClass >;
using MatrixKernelClass = FInterpMatrixKernelR;
using BalanceKernelClass= FChebBalanceSymKernel<CellClass, ContainerClass,
                                                MatrixKernelClass, 5,
                                                OctreeClass>;
// using KernelClass       = FChebSymKernel<CellClass, ContainerClass,
//                                          MatrixKernelClass, 5,
//                                          OctreeClass>;

template < template <typename...> class T>
using AlgoClass = T <OctreeClass, CellClass, ContainerClass, BalanceKernelClass, LeafClass >;


const FReal epsilon = 1e-4;


int main(int argc, char** argv)
{
    loadFMAAndRunFMMArgs args(argc, argv);

    FFmaGenericLoader loader(args.inFileName().c_str());
    OctreeClass tree(args.treeHeight(),
                     args.subTreeHeight(),
                     loader.getBoxWidth(),
                     loader.getCenterOfBox());

    loadTree(tree, loader);
   
    BalanceKernelClass kernel(&tree, epsilon);
    AlgoClass<FFmmAlgorithm> costAlgo(&tree, &kernel);

    costAlgo.execute();

    if (args.verboseLevel() > 1) {
        kernel.printResults(std::cout);
    }

    CostZones<OctreeClass, CellClass> costzones(&tree, args.zoneCount());
    costzones.run();

    writeZones(args, costzones);
    
    AlgoClass<FFmmAlgorithmThreadBalanced> fmmAlgo(&tree, &kernel, costzones.getZoneBounds());
    
    

    return EXIT_SUCCESS;
}


