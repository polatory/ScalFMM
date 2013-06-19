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

#include <string>

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FDebug.hpp"
#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"


int main(int /*argc*/,char** /*argv[]*/){
const int P = 10;
const int order = 1;
FPoint centerBox = FPoint(0,0,0);

typedef FTaylorCell<P,order> CellClass;
typedef FP2PParticleContainer ContainerClass;
typedef FTaylorKernel<CellClass,ContainerClass,P,order> KernelClass;
typedef FSimpleLeaf< ContainerClass > LeafClass;
typedef FOctree< CellClass, ContainerClass , LeafClass > OctreeClass;

KernelClass kernel(9,1.0,centerBox);

return 0;
}
