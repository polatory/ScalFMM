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
#ifndef FTREEBUILDER_H
#define FTREEBUILDER_H
#include <omp.h>

#include "../Utils/FLog.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FAssert.hpp"
#include "../Containers/FOctree.hpp"
#include "ScalFmmConfig.h"

/**
* @author Cyrille Piacibello
* @class FTreeBuilder
* @brief
* Please read the license
*
* This class provides a way to insert efficiently large amount of
* particles inside a tree.
*
* This is a static class. It's useless to instance it.  This class use
* the Threaded QuickSort or the output of FMpiTreeBuilder in order to
* sort the parts and insert them.
*
*/

template<class ParticleClass, class OctreeClass, class LeafClass>
class FTreeBuilder{

private:
    /**
     * This method has been taken from the octree class,
     * it computes a tree coordinate (x or y or z) from real cartesian position
     */
    static int GetTreeCoordinate(const FReal inRelativePosition, const FReal boxWidthAtLeafLevel,
				 const FReal boxWidth, const int height) {
	FAssertLF( (inRelativePosition >= 0 && inRelativePosition <= boxWidth), "inRelativePosition : ",inRelativePosition );
	if(inRelativePosition == boxWidth){
	    return FMath::pow2(height-1)-1;
	}
	const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
	return static_cast<int>(indexFReal);
    }


    /**
     * A particle may not have a MortonIndex Method (set/get morton
     * index) But in this algorithm they are sorted based on their
     * morton indexes.  So an IndexedParticle is storing a real
     * particle + its index.
     *
     */
    struct IndexedParticle{
    public:
	MortonIndex index;
	ParticleClass particle;

	operator MortonIndex() const {
	    return this->index;
	}
	bool operator<=(const IndexedParticle& rhs){
	    return this->index <= rhs.index;
	}
    };

public:
    /**
     *
     * Fill the tree with the parts from the array.
     * The Particles must provide a getPosition and a getPhysicalValue(), method.
     */
    template<class ContainerClass>
    static void BuildTreeFromArray(ContainerClass* arrayToBeInserted, const FSize numberOfParticle, OctreeClass * tree,
				   bool isAlreadySorted=false){
	//If the parts are already sorted, no need to sort again
	FLOG(FTic enumTimer);
	FLOG(FTic leavesPtr);
	FLOG(FTic leavesOffset);
	FLOG(FTic insertTimer);
	FLOG(FTic copyTimer);
	FLOG(FTic sortTimer);

	//General values needed
	int NbLevels = tree->getHeight();
	FPoint centerOfBox = tree->getBoxCenter();
	FReal boxWidth = tree->getBoxWidth();
	FReal boxWidthAtLeafLevel = boxWidth/FReal(1 << (NbLevels - 1));
	FPoint boxCorner   = centerOfBox - boxWidth/2;

	//We need to sort it in order to insert efficiently

	//First, copy datas into an array that will be sorted and
	//set Morton index for each particle

	//Temporary FTreeCoordinate
	FTreeCoordinate host;
	IndexedParticle * toBeSorted = new IndexedParticle[numberOfParticle];
	//Things to keep
	FPoint* posToBeIn = new FPoint[numberOfParticle];
	FReal * phyToBeIn = new  FReal[numberOfParticle];

	FLOG(copyTimer.tic());
	for(int idxParts=0; idxParts<numberOfParticle ; ++idxParts ){
	    toBeSorted[idxParts].particle = arrayToBeInserted->data()[idxParts];

	    host.setX( GetTreeCoordinate(arrayToBeInserted->data()[idxParts].getPosition().getX() - boxCorner.getX(),
					 boxWidthAtLeafLevel, boxWidth, NbLevels));
	    host.setY( GetTreeCoordinate(arrayToBeInserted->data()[idxParts].getPosition().getY() - boxCorner.getY(),
					 boxWidthAtLeafLevel, boxWidth, NbLevels ));
	    host.setZ( GetTreeCoordinate(arrayToBeInserted->data()[idxParts].getPosition().getZ() - boxCorner.getZ(),
					 boxWidthAtLeafLevel, boxWidth, NbLevels ));


	    toBeSorted[idxParts].index = host.getMortonIndex(NbLevels-1);

	}

	FLOG(copyTimer.tac());
	FLOG(FLog::Controller<<"Time needed for copying "<< numberOfParticle<<" particles : "<<copyTimer.elapsed() << " secondes !\n");

	if(!isAlreadySorted){
	    //Sort dat array
	    FLOG(sortTimer.tic());
	    FQuickSort<IndexedParticle,MortonIndex>::QsOmp(toBeSorted,numberOfParticle);
	    FLOG(sortTimer.tac());
	    FLOG(FLog::Controller << "Time needed for sorting the particles : "<< sortTimer.elapsed() << " secondes !\n");
	}



	//Enumerate the different leaves AND copy the positions
	unsigned int numberOfLeaves = 1;
	FLOG(enumTimer.tic());

	//First Values are copied
	posToBeIn[0] = toBeSorted[0].particle.getPosition();
	phyToBeIn[0] = toBeSorted[0].particle.getPhysicalValue();


	for(int idxParts = 1 ; idxParts < numberOfParticle ; ++idxParts){
	    //Copy
	    posToBeIn[idxParts] = toBeSorted[idxParts].particle.getPosition();
	    phyToBeIn[idxParts] = toBeSorted[idxParts].particle.getPhysicalValue();

	    if(toBeSorted[idxParts].index != toBeSorted[idxParts-1].index){
		numberOfLeaves++;
	    }
	}

	FLOG(enumTimer.tac());
	FLOG(FLog::Controller << "Time needed for enumerate the leaves : "<< enumTimer.elapsed() << " secondes !\n");
	FLOG(FLog::Controller << "Found " << numberOfLeaves << " leaves differents. \n");

	//Store the size of each leaves
	int * arrayOfSizeNbLeaves = new int[numberOfLeaves];
	memset(arrayOfSizeNbLeaves,0,sizeof(int)*(numberOfLeaves));
	//Init
	int indexInLeafArray = -1;
	arrayOfSizeNbLeaves[0] = 1;

	FLOG(leavesOffset.tic());
	MortonIndex currIndex = -1;

	for(int idxParts = 0 ; idxParts < numberOfParticle ; ++idxParts){
	    if(toBeSorted[idxParts].index == currIndex){
		arrayOfSizeNbLeaves[indexInLeafArray]++;
	    }
	    else{
		FAssertLF(FSize(indexInLeafArray)<numberOfLeaves,"Problem there : ",indexInLeafArray);
		indexInLeafArray++;
		//There is alway at least 1 part in each leaf
		arrayOfSizeNbLeaves[indexInLeafArray] = 1;
		currIndex = toBeSorted[idxParts].index;
	    }
	}
	//Debug
	int acc = 0;
	for(unsigned int i=0 ; i<numberOfLeaves ; ++i){
	    //printf("acc : %d arrayofsize[%d] = %d \n",acc,i,arrayOfSizeNbLeaves[i]);
	    acc += arrayOfSizeNbLeaves[i];
	}
	printf("Tot : %d/%lld\n",acc,numberOfParticle);

	FLOG(leavesOffset.tac());
	FLOG(FLog::Controller << "Time needed for setting the offset of each leaves : "<< leavesOffset.elapsed() << " secondes !\n");
	//Then, we create the leaves inside the tree

	//Idx of the first part in this leaf in the array of part
	int idxOfFirstPartInLeaf = 0;

	//struct to store leaves and idx
	struct LeafToFill{
	    LeafClass * leaf;
	    FSize idxOfLeafInPartArray;
	};

	FLOG(leavesPtr.tic());

	LeafToFill * leavesToFill = new LeafToFill[numberOfLeaves];
	memset(leavesToFill,0,sizeof(struct LeafToFill)*numberOfLeaves);


	for(FSize idxLeaf = 0; idxLeaf < numberOfLeaves ; ++idxLeaf){
	    leavesToFill[idxLeaf].leaf = tree->createLeaf(toBeSorted[idxOfFirstPartInLeaf].index);
	    leavesToFill[idxLeaf].idxOfLeafInPartArray = idxOfFirstPartInLeaf;
	    idxOfFirstPartInLeaf += arrayOfSizeNbLeaves[idxLeaf];
	}

	FLOG(leavesPtr.tac());
	FLOG(FLog::Controller << "Time needed for creating empty leaves : "<< leavesPtr.elapsed() << " secondes !\n");


	FLOG(insertTimer.tic());

	//Copy each parts into corresponding Leaf
#pragma omp parallel for schedule(auto)
	for(FSize idxLeaf=0 ; idxLeaf<numberOfLeaves ; ++idxLeaf ){
	    //Task consists in copy the parts inside the leaf
	    leavesToFill[idxLeaf].leaf->pushArray(&posToBeIn[leavesToFill[idxLeaf].idxOfLeafInPartArray],
						  arrayOfSizeNbLeaves[idxLeaf],
						  &phyToBeIn[leavesToFill[idxLeaf].idxOfLeafInPartArray]);
	}

	FLOG(insertTimer.tac());
	FLOG(FLog::Controller << "Time needed for inserting the parts into the leaves : "<< insertTimer.elapsed() << " secondes !\n");


	//Clean the mess
	delete [] leavesToFill;
	delete [] arrayOfSizeNbLeaves;
	delete [] posToBeIn;
	delete [] phyToBeIn;
	delete [] toBeSorted;
    }


};




#endif //FTREEBUILDER_H
