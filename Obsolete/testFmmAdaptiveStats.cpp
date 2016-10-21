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

// Keep in private GIT
// @SCALFMM_PRIVATE
#include <iostream>
#include <cstdio>


#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"
#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Adaptive/FAdaptiveCell.hpp"
#include "../../Src/Adaptive/FAdaptiveKernelWrapper.hpp"
#include "../../Src/Adaptive/FAbstractAdaptiveKernel.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

template< class CellClass, class ContainerClass>
class FAdaptiveStatsKernel : public FAbstractKernels<CellClass, ContainerClass>, public FAbstractAdaptiveKernel<CellClass, ContainerClass> {
    const int nbLevels;
    const int particlesThresh;

    int nbFmmOperations;
    int nbAdaptiveFmmOperations;
    int nbFmmOperationsVerbose[6];
    int nbAdaptiveFmmOperationsVerbose[8];
    int (*nbAdaptiveFmmOperationsVerbosePerLevel)[7];
    int (*nbAdaptiveFmmOperationsVerboseJump)[4];

public:
    FAdaptiveStatsKernel(const int inNbLevels, const int inParticlesThresh = 10)
        : nbLevels(inNbLevels), particlesThresh(inParticlesThresh) {

        nbFmmOperations = 0;
        nbAdaptiveFmmOperations = 0;

        for(int idx = 0 ; idx < 6 ; ++idx){
            nbFmmOperationsVerbose[idx] = 0;
        }
        for(int idx = 0 ; idx < 8 ; ++idx){
            nbAdaptiveFmmOperationsVerbose[idx] = 0;
        }

        nbAdaptiveFmmOperationsVerbosePerLevel = new int[nbLevels][7];
        memset(nbAdaptiveFmmOperationsVerbosePerLevel, 0, sizeof(int) * nbLevels * 7);

        nbAdaptiveFmmOperationsVerboseJump = new int[nbLevels][4];
        memset(nbAdaptiveFmmOperationsVerboseJump, 0, sizeof(int) * nbLevels * 4);
    }

    ~FAdaptiveStatsKernel(){
        delete[] nbAdaptiveFmmOperationsVerbosePerLevel;
        delete[] nbAdaptiveFmmOperationsVerboseJump;
    }

    void printStats() const {
        std::cout << "Nb Fmm Operations " << nbFmmOperations << "\n";
        std::cout << "Nb Adaptive Fmm Operations " << nbAdaptiveFmmOperations << "\n";

        std::cout << "Details Fmm :\n";
        std::cout << "\t P2M " << nbFmmOperationsVerbose[0] << "\n";
        std::cout << "\t M2M " << nbFmmOperationsVerbose[1] << "\n";
        std::cout << "\t M2L " << nbFmmOperationsVerbose[2] << "\n";
        std::cout << "\t L2L " << nbFmmOperationsVerbose[3] << "\n";
        std::cout << "\t L2P " << nbFmmOperationsVerbose[4] << "\n";
        std::cout << "\t P2P " << nbFmmOperationsVerbose[5] << "\n";

        std::cout << "Details Adaptive Fmm :\n";
        std::cout << "\t P2M " << nbAdaptiveFmmOperationsVerbose[0] << "\n";
        std::cout << "\t M2M " << nbAdaptiveFmmOperationsVerbose[1] << "\n";
        std::cout << "\t P2L " << nbAdaptiveFmmOperationsVerbose[2] << "\n";
        std::cout << "\t M2L " << nbAdaptiveFmmOperationsVerbose[3] << "\n";
        std::cout << "\t M2P " << nbAdaptiveFmmOperationsVerbose[4] << "\n";
        std::cout << "\t L2L " << nbAdaptiveFmmOperationsVerbose[5] << "\n";
        std::cout << "\t L2P " << nbAdaptiveFmmOperationsVerbose[6] << "\n";
        std::cout << "\t P2P " << nbAdaptiveFmmOperationsVerbose[7] << "\n";

        std::cout << "Details Adaptive Fmm Per Level:\n";
        for(int idxLevel = 1 ; idxLevel < nbLevels ; ++idxLevel){
            std::cout << "\t Level " << idxLevel << "\n";
            std::cout << "\t\t P2M " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][0] << "\n";
            std::cout << "\t\t M2M " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][1] << "\n";
            std::cout << "\t\t P2L " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][2] << "\n";
            std::cout << "\t\t M2L " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][3] << "\n";
            std::cout << "\t\t M2P " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][4] << "\n";
            std::cout << "\t\t L2L " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][5] << "\n";
            std::cout << "\t\t L2P " << nbAdaptiveFmmOperationsVerbosePerLevel[idxLevel][6] << "\n";
        }

        std::cout << "Details Adaptive Fmm Level difference for M2M:\n";
        for(int idxLevel = 1 ; idxLevel < nbLevels ; ++idxLevel){
            std::cout << "\t " << idxLevel << " \t" << nbAdaptiveFmmOperationsVerboseJump[idxLevel][0] << "\n";
        }
        std::cout << "Details Adaptive Fmm Level difference for M2L:\n";
        for(int idxLevel = 1 ; idxLevel < nbLevels ; ++idxLevel){
            std::cout << "\t " << idxLevel << " \t" <<  nbAdaptiveFmmOperationsVerboseJump[idxLevel][1] << "\t";
            std::cout << "\t " << -idxLevel << " \t" <<  nbAdaptiveFmmOperationsVerboseJump[idxLevel][2] << "\n";
        }
        std::cout << "Details Adaptive Fmm Level difference for L2L:\n";
        for(int idxLevel = 1 ; idxLevel < nbLevels ; ++idxLevel){
            std::cout << "\t " << idxLevel << " \t" <<  nbAdaptiveFmmOperationsVerboseJump[idxLevel][3] << "\n";
        }
    }

    void P2M(CellClass* const /*pole*/, const ContainerClass* const /*particles*/) override {
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[0] += 1;
    }

    void M2M(CellClass* const FRestrict /*pole*/, const CellClass *const FRestrict *const FRestrict /*child*/, const int /*level*/) override {
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[1] += 1;
    }

    void M2L(CellClass* const FRestrict /*pole*/, const CellClass* /*distantNeighbors*/[], const int /*neighborPositions*/[], const int /*size*/, const int /*level*/) override {
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[2] += 1;
    }

    void L2L(const CellClass*const FRestrict /*local*/, CellClass* FRestrict *const FRestrict /*child*/, const int /*level*/) override {
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[3] += 1;
    }

    void L2P(const CellClass* const  /*local*/, ContainerClass*const /*particles*/) override{
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[4] += 1;
    }

    void P2P(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict /*targets*/, const ContainerClass* const FRestrict /*sources*/,
                 ContainerClass* const /*directNeighborsParticles*/[], const int /*neighborPosition*/[], const int /*size*/) override{
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[5] += 1;
    }

    void P2POuter(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict /*targets*/,
                 ContainerClass* const /*directNeighborsParticles*/[], const int /*neighborPosition*/[], const int /*size*/) override{
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[5] += 1;
    }

    void P2PRemote(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict /*targets*/, const ContainerClass* const FRestrict /*sources*/,
                 const ContainerClass* const /*directNeighborsParticles*/[],  const int /*neighborPosition*/[], const int /*size*/) override{
        nbFmmOperations += 1;
        nbFmmOperationsVerbose[5] += 1;
    }

    void P2M(CellClass* const /*pole*/, const int cellLevel, const ContainerClass* const /*particles*/) override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[0] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[cellLevel][0] += 1;
    }

    void M2M(CellClass* const /*pole*/, const int poleLevel, const CellClass* const /*subCell*/, const int subCellLevel) override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[1] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[poleLevel][1] += 1;
        nbAdaptiveFmmOperationsVerboseJump[subCellLevel-poleLevel][0] += 1;
    }

    void P2L(CellClass* const /*local*/, const int localLevel, const ContainerClass* const /*particles*/) override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[2] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[localLevel][2] += 1;
    }

    void M2L(CellClass* const /*local*/, const int localLevel, const CellClass* const /*aNeighbor*/, const int neighborLevel) override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[3] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[localLevel][3] += 1;
        if(neighborLevel-localLevel >= 0 ){
            nbAdaptiveFmmOperationsVerboseJump[(neighborLevel-localLevel)][1] += 1;
        }
        else{
            nbAdaptiveFmmOperationsVerboseJump[(localLevel-neighborLevel)][2] += 1;
        }
    }

    void M2P(const CellClass* const /*pole*/, const int poleLevel, ContainerClass* const /*particles*/) override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[4] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[poleLevel][4] += 1;
    }

    void L2L(const CellClass* const /*local*/, const int localLevel, CellClass* const /*subCell*/, const int subCellLevel) override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[5] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[localLevel][5] += 1;
        nbAdaptiveFmmOperationsVerboseJump[subCellLevel-localLevel][3] += 1;
    }

    void L2P(const CellClass* const /*local*/, const int cellLevel, ContainerClass* const /*particles*/)  override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[6] += 1;
        nbAdaptiveFmmOperationsVerbosePerLevel[cellLevel][6] += 1;
    }

    void P2P(ContainerClass* /*target*/, const ContainerClass* /*sources*/)  override {
        nbAdaptiveFmmOperations += 1;
        nbAdaptiveFmmOperationsVerbose[7] += 1;
    }

    bool preferP2M(const ContainerClass* const particles) override {
        nbAdaptiveFmmOperations += 1;
        return particles->getNbParticles() < particlesThresh;
    }

    bool preferP2M(const int /*atLevel*/, const ContainerClass*const particles[], const int nbContainers) override {
        nbAdaptiveFmmOperations += 1;
        FSize counterParticles = 0;
        for(int idxContainer = 0 ; idxContainer < nbContainers ; ++idxContainer){
            counterParticles += particles[idxContainer]->getNbParticles();
        }
        return counterParticles < particlesThresh;
    }
};


/** This program show an example of use of the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test the adaptive FMM and print information about the real computation which is done.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;
    typedef FBasicCell                   CellClass;
    typedef FBasicParticleContainer<FReal, 0, FReal>   ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FAdaptiveStatsKernel< CellClass, ContainerClass >         KernelClass;
    typedef FAdaptiveCell< CellClass, ContainerClass >         CellWrapperClass;
    typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >         KernelWrapperClass;
    typedef FOctree< FReal, CellWrapperClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const int ParticlesThresh = FParameters::getValue(argc,argv,"-thresh", 10);
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoader<FReal> loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000000), 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    {
        FPoint<FReal> particlePosition;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
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

    KernelWrapperClass kernels(NbLevels, ParticlesThresh);
    FmmClass algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    kernels.getKernel().printStats();

    return 0;
}




