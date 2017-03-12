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

// ==== CMAKE =====
// @FUSE_MPI
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE


#include <iostream>

#include <cstdlib>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cstdio>

#include "ScalFmmConfig.h"
#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FMpi.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../Src/Utils/FLeafBalance.hpp"
#include "../Src/Containers/FTreeCoordinate.hpp"
#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"
#include "../../Src/Files/FTreeBuilder.hpp"

#include "Utils/FTic.hpp"
#include "Utils/FQuickSort.hpp"
#include "Utils/FParameters.hpp"
#include "../../Src/Utils/FParameterNames.hpp"

#include "Containers/FOctree.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif

#include "Utils/FTemplate.hpp"


/**
 * This program build a tree and insert the parts inside.
 * Time needed for the insert is outputed
 *
 */
int main(int argc, char** argv){

    typedef double FReal;

    struct TestParticle{
        MortonIndex index;
        FSize indexInFile;
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition()const{
            return position;
        }
        TestParticle& operator=(const TestParticle& other){
            index=other.index;
            indexInFile=other.indexInFile;
            position=other.position;
            physicalValue=other.physicalValue;
            return *this;
        }
        bool operator<=(const TestParticle& rhs)const{
            if(rhs.index < this->index){return false;}
            else{
                if(rhs.index > this->index){return true;}
                else{
                    if(rhs.indexInFile == this->indexInFile){
                        return true;
                    }
                    else {
                        return rhs.indexInFile> this->indexInFile ;
                    }
                }
            }
        }
    };

    FMpi app( argc, argv);

    static const int P = 9;

    typedef FRotationCell<FReal,P>               CellClass;
    typedef FP2PParticleContainer<FReal>          ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);

    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    const unsigned int NbThreads = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, omp_get_max_threads());

    omp_set_num_threads(NbThreads);
    std::cout << "Using " << omp_get_max_threads() <<" threads" << std::endl;
    std::cout << "Opening : " << filename << "\n";


    FMpiFmaGenericLoader<FReal> loaderRef(filename,app.global());
    if(!loaderRef.isOpen()){
        std::cout << "LoaderRef Error, " << filename << " is missing\n";
        return 1;
    }
    FTic regInsert;
    OctreeClass treeRef(NbLevels, SizeSubLevels, loaderRef.getBoxWidth(), loaderRef.getCenterOfBox());

    // -----------------------------------------------------
    {


        std::cout << "Creating & Inserting " << loaderRef.getNumberOfParticles() << " particles ..." << std::endl;
        std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;


        //Store the parts
        FmaRWParticle<FReal, 4,4>* particles = new FmaRWParticle<FReal,4,4>[loaderRef.getNumberOfParticles()];
        memset(particles, 0, sizeof(FmaRWParticle<FReal,4,4>) * loaderRef.getNumberOfParticles());

        for(FSize idxPart = 0 ; idxPart < loaderRef.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> pos;
            loaderRef.fillParticle(&pos,particles[idxPart].setPhysicalValue());
            particles[idxPart].setPosition(pos);
        }

        FVector<FmaRWParticle<FReal,4,4>> finalParticles;
        FLeafBalance balancer;
        FMpiTreeBuilder< FReal,FmaRWParticle<FReal,4,4> >::DistributeArrayToContainer(app.global(),particles,
                                                                                      loaderRef.getMyNumberOfParticles(),
                                                                                      treeRef.getBoxCenter(),
                                                                                      treeRef.getBoxWidth(),treeRef.getHeight(),
                                                                                      &finalParticles, &balancer);

        regInsert.tic();
        for(FSize idxPart = 0 ; idxPart < finalParticles.getSize() ; ++idxPart){
            treeRef.insert(finalParticles[idxPart].getPosition(),finalParticles[idxPart].getPhysicalValue() );
        }

        regInsert.tac();
        std::cout << "Time needed for regular insert : " << regInsert.elapsed() << " secondes" << std::endl;
    }
    //Second solution, parts must be sorted for that

    //Timers :
    FTic sortTimer;
    FTic insertTimer;
    FTic fillTimer;
    FTic enumTimer;
    FTic counter;
    FTic leavesOffset;
    FTic leavesPtr;

    FMpiFmaGenericLoader<FReal> loader(filename,app.global());
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    //Get the needed informations
    FSize nbOfParticles = loader.getNumberOfParticles();

    //Temporary TreeCoordinate
    FTreeCoordinate host;
    FmaRWParticle<FReal,4,4> * arrayOfParts = new FmaRWParticle<FReal,4,4>[nbOfParticles];
    memset(arrayOfParts,0,sizeof(FmaRWParticle<FReal,4,4>)*nbOfParticles);

    fillTimer.tic();

    for(FSize idxPart = 0 ; idxPart < nbOfParticles ; ++idxPart){
        FPoint<FReal> pos;
        loader.fillParticle(&pos,arrayOfParts[idxPart].setPhysicalValue());
        arrayOfParts[idxPart].setPosition(pos);
    }


    fillTimer.tac();
    std::cout << "Time needed for filling the array : "<< fillTimer.elapsed() << " secondes !" << std::endl;

    OctreeClass tree(NbLevels, SizeSubLevels, loaderRef.getBoxWidth(), loaderRef.getCenterOfBox());
    FVector< FmaRWParticle<FReal,4,4> > finalParticles;
    FLeafBalance balancer;
    FTic paraSort;
    paraSort.tic();
    FMpiTreeBuilder< FReal,FmaRWParticle<FReal,4,4> >::DistributeArrayToContainer(app.global(),arrayOfParts,
                                                                                  loader.getMyNumberOfParticles(),
                                                                                  tree.getBoxCenter(),
                                                                                  tree.getBoxWidth(),tree.getHeight(),
                                                                                  &finalParticles, &balancer);
    paraSort.tac();
    std::cout << "Time needed for FMpiTreeBuilder part : "<< paraSort.elapsed() << " secondes !" << std::endl;

    ContainerClass parts;
    parts.reserve(finalParticles.getSize());

    //Convert ouput of DistributeArrayToContainer to a ContainerClass
    for(FSize idxPart = 0; idxPart < finalParticles.getSize(); ++idxPart){

        parts.push(finalParticles[idxPart].getPosition(),finalParticles[idxPart].getPhysicalValue());
    }

    FTic treeBuilder;
    treeBuilder.tic();
    FTreeBuilder<FReal,OctreeClass,LeafClass>::BuildTreeFromArray(&tree,parts,true);
    treeBuilder.tac();
    std::cout << "Time needed for TreeBuilder : "<< treeBuilder.elapsed() << " secondes !" << std::endl;
#define CHECK_TREE
    //Check the datas
#ifdef CHECK_TREE

    FILE * fd = fopen("Tree","a+");
    FSize acc = 0;
    int accL = 0;
    tree.forEachLeaf([&](LeafClass* leaf){
        ContainerClass * partSrc = leaf->getSrc();
        FSize nbPa = partSrc->getNbParticles();
        acc += nbPa;
        accL++;
        const FReal * xPos = partSrc->getPositions()[0];
        const FReal * yPos = partSrc->getPositions()[1];
        const FReal * zPos = partSrc->getPositions()[2];
        const FReal * phyVal = partSrc->getAttribute(0);
        fprintf(fd,"Proc : %d, containing %lld \n",
                app.global().processId(),nbPa);
        for(int idP = 0 ; idP < nbPa ; ++idP){
            fprintf(fd,"[%d] : %e,%e,%e, phy : %e\n",
                    app.global().processId(), xPos[idP],yPos[idP],zPos[idP],phyVal[idP]);
        }
    });
    printf("Number of parts in MY Tree : %lld, Number of leaves : %d \n",acc,accL);
    fclose(fd);

    FSize accRef = 0;
    int accRefL = 0;
    FILE * fd2 = fopen("TreeRef","a+");
    treeRef.forEachLeaf([&](LeafClass* leaf){
        ContainerClass * partSrc = leaf->getSrc();
        FSize nbPa = partSrc->getNbParticles();
        accRef+=nbPa;
        accRefL++;
        const FReal * xPos = partSrc->getPositions()[0];
        const FReal * yPos = partSrc->getPositions()[1];
        const FReal * zPos = partSrc->getPositions()[2];
        const FReal * phyVal = partSrc->getAttribute(0);
        fprintf(fd2,"Proc : %d, containing %lld \n",
                app.global().processId(),nbPa);
        for(int idP = 0 ; idP < nbPa ; ++idP){
            fprintf(fd2,"[%d] : %e,%e,%e, phy : %e\n",app.global().processId(), xPos[idP],yPos[idP],zPos[idP],phyVal[idP]);
        }
    });
    fclose(fd2);
    printf("Number of parts in REG Tree : %lld, number of leaves %d\n",accRef,accRefL);
#endif //CHECK_TREE
    delete [] arrayOfParts;

    return 0;
}
