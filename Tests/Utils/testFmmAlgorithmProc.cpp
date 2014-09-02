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
// @FUSE_MPI
// ================

#include "../../Src/Utils/FMpi.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Components/FTestParticleContainer.hpp"

//#include "../../Src/Core/FFmmAlgorithmProcMpi.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"

#include <iostream>
#include <cstdio>
#include <cstdlib>


/** This program show an example of use of the fmm threaded + mpi algo
 * it also check that each particles is impacted each other particles
 */

/////////////////////////////////////////////////////////////////////////////
// Test function
/////////////////////////////////////////////////////////////////////////////

// Check if tree is built correctly
template<class OctreeClass>
void ValidateTree(OctreeClass& realTree,
                  OctreeClass& treeValide){

    typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
    octreeIteratorValide.gotoBottomLeft();

    typename OctreeClass::Iterator octreeIterator(&realTree);
    octreeIterator.gotoBottomLeft();

    while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
        if(octreeIteratorValide.moveRight() == false){
            std::cout << "Error the real tree smaller than the parallel one\n";
            std::cout << "Valide tree stop at " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";
            std::cout << "Other at " << octreeIterator.getCurrentGlobalIndex() << "\n";
            return;
        }
    }

    std::cout << "The tree starts at " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";

    while(true){
        if(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            std::cout << "Error the trees do not have the same idx.\n";
            std::cout << "Correct one is " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";
            std::cout << "Incorrect one is " << octreeIterator.getCurrentGlobalIndex() << "\n";
            return;
        }

        if(octreeIteratorValide.getCurrentListSrc()->getNbParticles() != octreeIterator.getCurrentListSrc()->getNbParticles()){
            std::cout << "Error the trees do not have the same number of particles at idx " << octreeIteratorValide.getCurrentGlobalIndex() << ".\n";
            std::cout << "Correct one is " << octreeIteratorValide.getCurrentListSrc()->getNbParticles() << "\n";
            std::cout << "Incorrect one is " << octreeIterator.getCurrentListSrc()->getNbParticles() << "\n";
            return;
        }

        if(octreeIterator.moveRight() == false){
            break;
        }

        if(octreeIteratorValide.moveRight() == false){
            std::cout << "Error the real tree smaller than the parallel one\n";
        }
    }

    std::cout << "The tree stops at " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";
}



/** This function tests the octree to be sure that the fmm algorithm
 * has worked completly.
 */
template<class OctreeClass, class ContainerClass, class FmmClassProc>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree,
                         FmmClassProc* const fmm){
    std::cout << "The working interval is from " << fmm->getWorkingInterval(badTree->getHeight()-1).leftIndex << "\n";
    std::cout << "The working interval is to " << fmm->getWorkingInterval(badTree->getHeight()-1).rightIndex << "\n";
    std::cout << "\tValidate L2L M2M...\n";
    const int OctreeHeight = badTree->getHeight();
    {
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 0 && fmm->hasWorkAtLevel(level) ; --level){

            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()) {
                octreeIteratorValide.moveRight();
            }

            while(octreeIteratorValide.getCurrentGlobalIndex() != fmm->getWorkingInterval(level).leftIndex){
                octreeIteratorValide.moveRight();
                octreeIterator.moveRight();
            }

            FSize countCheck = 0;
            do{
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Problem Error index are not equal! " << octreeIterator.getCurrentGlobalIndex() << " " << octreeIteratorValide.getCurrentGlobalIndex() << std::endl;
                }
                else{
                    if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                        std::cout << "Problem M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                  << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                    if(octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorValide.getCurrentCell()->getDataDown()){
                        std::cout << "Problem L2L error at level " << level << " down bad " << octreeIterator.getCurrentCell()->getDataDown()
                                  << " good " << octreeIteratorValide.getCurrentCell()->getDataDown() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                }
                ++countCheck;
            } while(octreeIteratorValide.moveRight() && octreeIterator.moveRight());

            // Check that each particle has been summed with all other

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    std::cout << "\tValidate L2P P2P...\n";
    {
        FSize NbPart = 0;
        FSize NbLeafs = 0;
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(valideTree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSrc()->getNbParticles();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
        }
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

                ContainerClass* container = (octreeIterator.getCurrentListTargets());
                const long long int*const dataDown = container->getDataDown();

                for(FSize idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( (!isUsingTsm && dataDown[idxPart] != NbPart - 1) ||
                            (isUsingTsm && dataDown[idxPart] != NbPart) ){
                        std::cout << "Problem L2P + P2P, value on particle is : " << dataDown[idxPart] <<
                                     " at pos " << idxPart << " index is " << octreeIterator.getCurrentGlobalIndex() << "\n";
                    }
                }
            } while( octreeIterator.moveRight());
        }
    }
    std::cout << "\tValidate P2M...\n";
    {
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                if(octreeIterator.getCurrentListSrc()->getNbParticles() != octreeIterator.getCurrentCell()->getDataUp()){
                    printf("P2M Problem nb part %d data up %lld \n",
                           octreeIterator.getCurrentListSrc()->getNbParticles(), octreeIterator.getCurrentCell()->getDataUp());
                }
            } while( octreeIterator.moveRight() );
        }
    }
    std::cout << "\tValidate Particles...\n";
    {
        // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator valideOctreeIterator(valideTree);
        valideOctreeIterator.gotoBottomLeft();
        while(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            valideOctreeIterator.moveRight();
        }

        do {
            if(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                printf("Problem Do not have the same index valide %lld invalide %lld \n",
                       valideOctreeIterator.getCurrentGlobalIndex(), octreeIterator.getCurrentGlobalIndex());
                break;
            }

            if(octreeIterator.getCurrentListTargets()->getNbParticles() != valideOctreeIterator.getCurrentListTargets()->getNbParticles()){
                printf("Problem Do not have the same number of particle at leaf id %lld, valide %d invalide %d \n",
                       octreeIterator.getCurrentGlobalIndex(), valideOctreeIterator.getCurrentListTargets()->getNbParticles(),
                       octreeIterator.getCurrentListTargets()->getNbParticles());
            }
            else {
                ContainerClass* container = (octreeIterator.getCurrentListTargets());
                const long long int*const dataDown = container->getDataDown();

                ContainerClass* containerValide = (valideOctreeIterator.getCurrentListTargets());
                const long long int*const dataDownValide = containerValide->getDataDown();

                for(int idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( dataDown[idxPart] != dataDownValide[idxPart]){
                        std::cout << "Problem on leaf " << octreeIterator.getCurrentGlobalIndex() <<
                                     " part " << idxPart << " valide data down " << dataDownValide[idxPart] <<
                                     " invalide " << dataDown[idxPart] << "\n";
                        std::cout << "Data down for leaf is: valide " << valideOctreeIterator.getCurrentCell()->getDataDown()
                                  << " invalide " << octreeIterator.getCurrentCell()->getDataDown()
                                  << " size is: valide " <<  valideOctreeIterator.getCurrentListTargets()->getNbParticles()
                                  << " invalide " << octreeIterator.getCurrentListTargets()->getNbParticles() << std::endl;
                    }
                }
            }

        }while( octreeIterator.moveRight() && valideOctreeIterator.moveRight());
    }
    std::cout << "\tDone!\n";
}


/** To print an octree
 * used to debug and understand how the values were passed
 */
template<class OctreeClass>
void print(OctreeClass* const valideTree){
    typename OctreeClass::Iterator octreeIterator(valideTree);
    for(int idxLevel = valideTree->getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
        do{
            std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
        } while(octreeIterator.moveRight());
        std::cout << "\n";
        octreeIterator.gotoLeft();
        octreeIterator.moveDown();
    }
}


/////////////////////////////////////////////////////////////////////
// Define the classes to use
/////////////////////////////////////////////////////////////////////

typedef FTestCell                  CellClass;
typedef FTestParticleContainer     ContainerClass;

typedef FSimpleLeaf< ContainerClass >                     LeafClass;
typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
typedef FFmmAlgorithmThreadProc<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassProc;

/////////////////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
                "../../Data/test20k.bin.fma.single":
                "../../Data/test20k.bin.fma.double";
    const char* const filename = FParameters::getStr(argc,argv,"-f", defaultFilename);
    std::cout << "Opening : " << filename << "\n";

    std::ifstream loader;
    loader.open(filename, std::ios::in | std::ios::binary);

    int nbParts;
    loader.read((char*)&nbParts, sizeof(nbParts));
    FReal boxWidthX, boxWidthY, boxWidthZ;
    loader.read((char*)&boxWidthX, sizeof(boxWidthX));
    loader.read((char*)&boxWidthY, sizeof(boxWidthY));
    loader.read((char*)&boxWidthZ, sizeof(boxWidthZ));
    FReal boxCenterX, boxCenterY, boxCenterZ;
    loader.read((char*)&boxCenterX, sizeof(boxCenterX));
    loader.read((char*)&boxCenterY, sizeof(boxCenterY));
    loader.read((char*)&boxCenterZ, sizeof(boxCenterZ));


    std::cout << "Simulation properties :\n";
    std::cout << "Nb Particles " << nbParts << "\n";
    std::cout << "Box Width : " << boxWidthX << "\n";
    std::cout << "Box Center : " << FPoint(boxCenterX, boxCenterY, boxCenterZ) << "\n";

    // The real tree to work on
    OctreeClass realTree(NbLevels, SizeSubLevels,boxWidthX,FPoint(boxCenterX, boxCenterY, boxCenterZ));

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;
        counter.tic();

        struct TestParticle{
            FPoint position;
            const FPoint& getPosition(){
                return position;
            }
        };

        TestParticle* particles = new TestParticle[nbParts];
        memset(particles, 0, sizeof(TestParticle) * nbParts);
        for(int idxPart = 0 ; idxPart < nbParts ; ++idxPart){
            FReal px, py, pz, pp;
            loader.read((char*)&px, sizeof(px));
            loader.read((char*)&py, sizeof(py));
            loader.read((char*)&pz, sizeof(pz));
            loader.read((char*)&pp, sizeof(pp));
            particles[idxPart].position = FPoint(px,py,pz);
        }

        const int leftPart = app.global().getLeft(nbParts);
        const int rightPart = app.global().getRight(nbParts);

        std::cout << "Go from " << leftPart << " to " << rightPart << "\n";

        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;
        FMpiTreeBuilder< TestParticle >::DistributeArrayToContainer(app.global(),&particles[leftPart],
                                                                    rightPart-leftPart,
                                                                    realTree.getBoxCenter(),
                                                                    realTree.getBoxWidth(),realTree.getHeight(),
                                                                    &finalParticles, &balancer);
        std::cout << "I have now " << finalParticles.getSize() << " particles\n";

        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            realTree.insert(finalParticles[idx].position);
        }

        delete[] particles;

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        const FSize nbParticles = nbParts;
        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            FReal px, py, pz, pp;
            loader.read((char*)&px, sizeof(px));
            loader.read((char*)&py, sizeof(py));
            loader.read((char*)&pz, sizeof(pz));
            loader.read((char*)&pp, sizeof(pp));
            realTree.insert(FPoint(px,py,pz));
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Create real tree
    //////////////////////////////////////////////////////////////////////////////////

    OctreeClass treeValide(NbLevels, SizeSubLevels,boxWidthX,FPoint(boxCenterX, boxCenterY, boxCenterZ));
    {
        std::ifstream myfile;
        myfile.open(filename, std::ios::in | std::ios::binary);

        int nbParts2;
        myfile.read((char*)&nbParts2, sizeof(nbParts2));
        FReal boxWidthX2, boxWidthY2, boxWidthZ2;
        myfile.read((char*)&boxWidthX2, sizeof(boxWidthX2));
        myfile.read((char*)&boxWidthY2, sizeof(boxWidthY2));
        myfile.read((char*)&boxWidthZ2, sizeof(boxWidthZ2));
        FReal boxCenterX2, boxCenterY2, boxCenterZ2;
        myfile.read((char*)&boxCenterX2, sizeof(boxCenterX2));
        myfile.read((char*)&boxCenterY2, sizeof(boxCenterY2));
        myfile.read((char*)&boxCenterZ2, sizeof(boxCenterZ2));

        std::cout << "Simulation properties :\n";
        std::cout << "Nb Particles " << nbParts2 << "\n";
        std::cout << "Box Width : " << boxWidthX2 << "\n";
        std::cout << "Box Center : " << FPoint(boxCenterX2, boxCenterY2, boxCenterZ2) << "\n";

        for(FSize idxPart = 0 ; idxPart < nbParts ; ++idxPart){
            FReal px, py, pz, pp;
            myfile.read((char*)&px, sizeof(px));
            myfile.read((char*)&py, sizeof(py));
            myfile.read((char*)&pz, sizeof(pz));
            myfile.read((char*)&pp, sizeof(pp));
            treeValide.insert(FPoint(px,py,pz));
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Validate tree ..." << std::endl;
    counter.tic();

    ValidateTree(realTree, treeValide);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working parallel particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;

    FmmClassProc algo(app.global(),&realTree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working sequential particles ..." << std::endl;
    counter.tic();

    FmmClass algoValide(&treeValide,&kernels);
    algoValide.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Checking data ..." << std::endl;
    counter.tic();

    ValidateFMMAlgoProc<OctreeClass,ContainerClass, FmmClassProc>(&realTree,&treeValide,&algo);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



