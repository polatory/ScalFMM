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

#include <cstdlib>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cstdio>

#include "ScalFmmConfig.h"

#include "Files/FFmaGenericLoader.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../Src/BalanceTree/FLeafBalance.hpp"
#include "../Src/Containers/FTreeCoordinate.hpp"

#include "Utils/FTic.hpp"
#include "Utils/FQuickSort.hpp"
#include "Utils/FParameters.hpp"
#include "../../Src/Utils/FParameterNames.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FCoordinateComputer.hpp"

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


    typedef double FReal;
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


    FFmaGenericLoader<FReal> loaderRef(filename);
    if(!loaderRef.isOpen()){
        std::cout << "LoaderRef Error, " << filename << " is missing\n";
        return 1;
    }
    FTic regInsert;
    // -----------------------------------------------------
    {
        OctreeClass treeRef(NbLevels, SizeSubLevels, loaderRef.getBoxWidth(), loaderRef.getCenterOfBox());


        std::cout << "Creating & Inserting " << loaderRef.getNumberOfParticles() << " particles ..." << std::endl;
        std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

        regInsert.tic();

        for(FSize idxPart = 0 ; idxPart < loaderRef.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> particlePosition;
            FReal physicalValue;
            loaderRef.fillParticle(&particlePosition,&physicalValue);
            treeRef.insert(particlePosition, physicalValue );
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

    FFmaGenericLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    //Get the needed informations
    FReal boxWidth = loader.getBoxWidth();
    FReal boxWidthAtLeafLevel = boxWidth/FReal(1 << (NbLevels - 1));
    FPoint<FReal> centerOfBox = loader.getCenterOfBox();
    FPoint<FReal> boxCorner   = centerOfBox - boxWidth/2;
    FSize nbOfParticles = loader.getNumberOfParticles();

    //Temporary TreeCoordinate
    FTreeCoordinate host;
    TestParticle * arrayOfParts = new TestParticle[nbOfParticles];
    memset(arrayOfParts,0,sizeof(TestParticle)*nbOfParticles);

    fillTimer.tic();

    for(FSize idxPart = 0 ; idxPart < nbOfParticles ; ++idxPart){
        loader.fillParticle(&arrayOfParts[idxPart].position,&arrayOfParts[idxPart].physicalValue);
        //Build temporary TreeCoordinate
        host.setX( FCoordinateComputer::GetTreeCoordinate<FReal>( arrayOfParts[idxPart].getPosition().getX() - boxCorner.getX(), boxWidth, boxWidthAtLeafLevel, NbLevels ));
        host.setY( FCoordinateComputer::GetTreeCoordinate<FReal>( arrayOfParts[idxPart].getPosition().getY() - boxCorner.getY(), boxWidth, boxWidthAtLeafLevel, NbLevels ));
        host.setZ( FCoordinateComputer::GetTreeCoordinate<FReal>( arrayOfParts[idxPart].getPosition().getZ() - boxCorner.getZ(), boxWidth, boxWidthAtLeafLevel, NbLevels ));

        //Set Morton index from Tree Coordinate
        arrayOfParts[idxPart].index = host.getMortonIndex();
        arrayOfParts[idxPart].indexInFile = idxPart;

    }


    fillTimer.tac();
    std::cout << "Time needed for filling the array : "<< fillTimer.elapsed() << " secondes !" << std::endl;

    sortTimer.tic();
    //std::sort(arrayOfParts,&arrayOfParts[nbOfParticles-1]);
    FQuickSort<TestParticle,MortonIndex>::QsOmp(arrayOfParts,nbOfParticles);
    sortTimer.tac();

    std::cout << "Time needed for sorting the array : "<< sortTimer.elapsed() << " secondes !" << std::endl;


    //Enumerate the different leaves
    unsigned int numberOfLeaves = 0;
    enumTimer.tic();

    for(FSize idxParts = 1 ; idxParts < nbOfParticles ; ++idxParts){
        if(arrayOfParts[idxParts].index != arrayOfParts[idxParts-1].index){
            numberOfLeaves++;
        }
    }
    enumTimer.tac();
    std::cout << "Time needed for enumerate the leaves : "<< enumTimer.elapsed() << " secondes !" << std::endl;
    std::cout << "Found " << numberOfLeaves << " leaves differents." << std::endl;

    //Store the size of each subOctree
    FSize * arrayOfSizeNbLeaves = new FSize[numberOfLeaves];
    memset(arrayOfSizeNbLeaves,0,sizeof(FSize)*(numberOfLeaves));
    //Init
    int indexInLeafArray = 0;
    arrayOfSizeNbLeaves[0] = 1;

    leavesOffset.tic();
    for(FSize idxParts = 1 ; idxParts < nbOfParticles ; ++idxParts){
        if(arrayOfParts[idxParts].index == arrayOfParts[idxParts-1].index){
            arrayOfSizeNbLeaves[indexInLeafArray]++;
        }
        else{
            indexInLeafArray++;
        }
    }
    leavesOffset.tac();

    std::cout << "Time needed for setting the offset of each leaves : "<< leavesOffset.elapsed() << " secondes !" << std::endl;

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    //Then, we create the leaves inside the tree

    //Idx of the first part in this leaf in the array of part
    FSize idxOfFirstPartInLeaf =	0;

    struct LeafToFill{
        LeafClass * leaf;
        FSize idxOfLeafInPartArray;
    };

    leavesPtr.tic();

    LeafToFill * leavesToFill = new LeafToFill[numberOfLeaves];
    memset(leavesToFill,0,sizeof(struct LeafToFill)*numberOfLeaves);
    insertTimer.tic();


    for(FSize idxLeaf = 0; idxLeaf < numberOfLeaves ; ++idxLeaf){
        leavesToFill[idxLeaf].leaf = tree.createLeaf(arrayOfParts[idxOfFirstPartInLeaf].index);
        idxOfFirstPartInLeaf += arrayOfSizeNbLeaves[idxLeaf];
        leavesToFill[idxLeaf].idxOfLeafInPartArray += idxOfFirstPartInLeaf;
    }

    leavesPtr.tac();
    std::cout << "Time needed for creating empty leaves : "<< leavesPtr.elapsed() << " secondes !" << std::endl;


    #pragma omp parallel for schedule(auto)
    for(FSize idxLeaf=0 ; idxLeaf<numberOfLeaves ; ++idxLeaf ){
        //Task consists in copy the parts inside the leaf
        leavesToFill[idxLeaf].leaf->pushArray(&arrayOfParts[leavesToFill[idxLeaf].idxOfLeafInPartArray].position,
                arrayOfSizeNbLeaves[idxLeaf],
                &arrayOfParts[leavesToFill[idxLeaf].idxOfLeafInPartArray].physicalValue);

    }


    insertTimer.tac();

    std::cout << "Time needed for inserting the parts " << insertTimer.elapsed() << " secondes !" << std::endl;
    double totalTimeForThatStrat = insertTimer.elapsed() + sortTimer.elapsed()
            + fillTimer.elapsed() + enumTimer.elapsed()
            + leavesPtr.elapsed() + leavesOffset.elapsed();
    std::cout << "Total time for that strat : " << totalTimeForThatStrat << " secondes !" << std::endl;
    std::cout << "---------Total------------> Difference : " << regInsert.elapsed() - totalTimeForThatStrat << " secondes !" << std::endl;
    std::cout << "-------On Insert----------> Difference : " << regInsert.elapsed() - insertTimer.elapsed() << " secondes !" << std::endl;


    delete [] arrayOfSizeNbLeaves;
    delete [] arrayOfParts;

    return 0;
}
