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

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/Adaptive/FAdaptiveCell.hpp"
#include "../../Src/Adaptive/FAdaptiveKernelWrapper.hpp"
#include "../../Src/Adaptive/FAbstractAdaptiveKernel.hpp"
#include "../../Src/Adaptive/FAdaptiveTestKernel.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


/** This program show an example of use of the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    const FParameterNames PrintTree =  { {"-print-trees"}, "Print the details of the trees."};
    FHelpDescribeAndExit(argc, argv,
                         "Test the adaptive FMM.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,PrintTree);

    typedef double FReal;
    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FAdaptiveTestKernel< CellClass, ContainerClass >         KernelClass;
    typedef FAdaptiveCell< CellClass, ContainerClass >         CellWrapperClass;
    typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >         KernelWrapperClass;
    typedef FOctree< FReal, CellWrapperClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

    typedef FTestCell                   CellClassTest;
    typedef FTestParticleContainer<FReal>      ContainerClassTest;
    typedef FSimpleLeaf<FReal, ContainerClassTest >                     LeafClassTest;
    typedef FOctree<FReal, CellClassTest, ContainerClassTest , LeafClassTest >  OctreeClassTest;
    typedef FTestKernels< CellClassTest, ContainerClassTest >         KernelClassTest;
    typedef FFmmAlgorithm<OctreeClassTest, CellClass, ContainerClassTest, KernelClassTest, LeafClassTest >     FmmClassTest;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    //FRandomLoader<FReal> loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000000), 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    FFmaGenericLoader<FReal> loader(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options,"../Data/noDist/prolate50-ref.fma"));
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    OctreeClassTest treeTest(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    {
        FPoint<FReal> particlePosition;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FReal pv;
            loader.fillParticle(&particlePosition, &pv);
            //loader.fillParticle(&particlePosition);
            tree.insert(particlePosition);
            treeTest.insert(particlePosition);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelWrapperClass kernels(4);            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    KernelClassTest kernelsTest;            // FTestKernels FBasicKernels
    FmmClassTest algoTest(&treeTest,&kernelsTest);  //FFmmAlgorithm FFmmAlgorithmThread
    algoTest.execute();

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    tree.forEachCellLeaf([&](CellWrapperClass*, LeafClass* leaf){
        long long int*const particlesAttributes = leaf->getTargets()->getDataDown();
        for(FSize idxPart = 0 ; idxPart < leaf->getTargets()->getNbParticles() ; ++idxPart){
            if(particlesAttributes[idxPart] != (loader.getNumberOfParticles()-1)){
                std::cout << "Incorrect " << particlesAttributes[idxPart] << " instead of " << (loader.getNumberOfParticles()-1) << "\n";
            }
        }
    });


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    tree.forEachCellWithLevel([&](CellWrapperClass* cell, int level){
        if(cell->hasDevelopment()){
            FAssertLF(cell->isAdaptive());
            CellClass* testcell = cell->getRealCell();
            FAssertLF(testcell);
            CellClassTest* testcelltest = treeTest.getCell(testcell->getMortonIndex(), level);
            FAssertLF(testcelltest);
            if(testcell->getDataUp() != testcelltest->getDataUp()){
                std::cout << "Error at index " << testcell->getMortonIndex() << " level " << level << " up is " << testcell->getDataUp()
                             << " should be " << testcelltest->getDataUp() << "\n";
            }
            if(testcell->getDataDown() != testcelltest->getDataDown()){
                std::cout << "Error at index " << testcell->getMortonIndex() << " level " << level << " up is " << testcell->getDataDown()
                             << " should be " << testcelltest->getDataDown() << "\n";
            }
        }
    });


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    if(FParameters::existParameter(argc, argv, PrintTree.options)){
        int previousLevel = 0;
        tree.forEachCellWithLevel([&](CellWrapperClass* cell, int level){
            if(cell->hasDevelopment()){
                if(previousLevel != level){
                    previousLevel = level;
                    std::cout << "\n" << level << "] ";
                }
                FAssertLF(cell->isAdaptive());
                CellClass* testcell = cell->getRealCell();
                std::cout << "Idx:" << testcell->getMortonIndex() << " Up " << testcell->getDataUp() << "\t";
            }
        });
        std::cout << "\n";

        previousLevel = 0;
        treeTest.forEachCellWithLevel([&](CellClassTest* cell, int level){
            if(previousLevel != level){
                previousLevel = level;
                std::cout << "\n" << level << "] ";
            }
            std::cout << "Idx:" << cell->getMortonIndex() << " Up " << cell->getDataUp() << "\t";
        });
        std::cout << "\n";
    }

    return 0;
}




