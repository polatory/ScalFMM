// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
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

//
#include <iostream>
#include <cstdio>


#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "Adaptative/FAdaptiveCell.hpp"
#include "Adaptative/FAdaptiveKernelWrapper.hpp"
#include "Adaptative/FAbstractAdaptiveKernel.hpp"

template< class CellClass, class ContainerClass>
class FAdaptiveTestKernel : public FTestKernels<CellClass, ContainerClass>, public FAbstractAdaptiveKernel<CellClass, ContainerClass> {
public:
    using FTestKernels<CellClass, ContainerClass>::P2M;
    using FTestKernels<CellClass, ContainerClass>::M2M;
    using FTestKernels<CellClass, ContainerClass>::M2L;
    using FTestKernels<CellClass, ContainerClass>::L2L;
    using FTestKernels<CellClass, ContainerClass>::L2P;
    using FTestKernels<CellClass, ContainerClass>::P2P;

    void P2M(CellClass* const pole, const int /*cellLevel*/, const ContainerClass* const particles) override {
        pole->setDataUp(pole->getDataUp() + particles->getNbParticles());
    }

    void M2M(CellClass* const pole, const int /*poleLevel*/, const CellClass* const subCell, const int /*subCellLevel*/) override {
        pole->setDataUp(pole->getDataUp() + subCell->getDataUp());
    }

    void P2L(CellClass* const local, const int /*localLevel*/, const ContainerClass* const particles) override {
        local->setDataDown(local->getDataDown() + particles->getNbParticles());
    }

    void M2L(CellClass* const local, const int /*localLevel*/, const CellClass* const aNeighbor, const int /*neighborLevel*/) override {
        local->setDataDown(local->getDataDown() + aNeighbor->getDataUp());
    }

    void M2P(const CellClass* const pole, const int /*poleLevel*/, ContainerClass* const particles) override {
        long long int*const particlesAttributes = particles->getDataDown();
        for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += pole->getDataUp();
        }
    }

    void L2L(const CellClass* const local, const int /*localLevel*/, CellClass* const subCell, const int /*subCellLevel*/) override {
        subCell->setDataDown(local->getDataDown() + subCell->getDataDown());
    }

    void L2P(const CellClass* const local, const int /*cellLevel*/, ContainerClass* const particles)  override {
        long long int*const particlesAttributes = particles->getDataDown();
        for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += local->getDataDown();
        }
    }

    void P2P(ContainerClass* target, const ContainerClass* sources)  override {
        long long int*const particlesAttributes = target->getDataDown();
        for(int idxPart = 0 ; idxPart < target->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += sources->getNbParticles();
        }
    }

    bool preferP2M(const ContainerClass* const particles) override {
        return particles->getNbParticles() < 10;
    }
    bool preferP2M(const int /*atLevel*/, const ContainerClass*const particles[], const int nbContainers) override {
        int counterParticles = 0;
        for(int idxContainer = 0 ; idxContainer < nbContainers ; ++idxContainer){
            counterParticles += particles[idxContainer]->getNbParticles();
        }
        return counterParticles < 10;
    }
};


/** This program show an example of use of the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer      ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FAdaptiveTestKernel< CellClass, ContainerClass >         KernelClass;
    typedef FAdaptiveCell< CellClass, ContainerClass >         CellWrapperClass;
    typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >         KernelWrapperClass;
    typedef FOctree< CellWrapperClass, ContainerClass , LeafClass >  OctreeClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const int NbPart       = FParameters::getValue(argc,argv,"-nb", 2000000);
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoader loader(NbPart, 1, FPoint(0.5,0.5,0.5), 1);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    {
        FPoint particlePosition;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition);
            tree.insert(particlePosition);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelWrapperClass kernels;            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    tree.forEachCellLeaf([&](CellWrapperClass*, LeafClass* leaf){
        long long int*const particlesAttributes = leaf->getTargets()->getDataDown();
        for(int idxPart = 0 ; idxPart < leaf->getTargets()->getNbParticles() ; ++idxPart){
            if(particlesAttributes[idxPart] != (NbPart-1)){
                printf("Incorrect %lld instead of %d\n", particlesAttributes[idxPart], (NbPart-1));
            }
        }
    });

    return 0;
}




