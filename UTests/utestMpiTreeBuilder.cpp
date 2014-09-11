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

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include "ScalFmmConfig.h"
#include <cstdlib>
#include <string.h>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include "FUTester.hpp"

#include "../Src/Utils/FMpi.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Files/FFmaGenericLoader.hpp"
#include "../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../Src/BalanceTree/FLeafBalance.hpp"
#include "../Src/Containers/FTreeCoordinate.hpp"


#include "../Src/Utils/FQuickSortMpi.hpp"
#include "../Src/Utils/FBitonicSort.hpp"
#include "../Src/Files/FMpiTreeBuilder.hpp"
#include "../Src/Core/FCoreCommon.hpp"

#include "../Src/Utils/FPoint.hpp"
#include "../Src/Utils/FMath.hpp"


class TestMpiTreeBuilder :  public FUTesterMpi< class TestMpiTreeBuilder> {
    /** This method has been tacken from the octree
     * it computes a tree coordinate (x or y or z) from real position
     */
    static int getTreeCoordinate(const FReal inRelativePosition, const FReal boxWidthAtLeafLevel) {
	const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
	const int index = int(FMath::dfloor(indexFReal));
	if( index && FMath::LookEqual(inRelativePosition, boxWidthAtLeafLevel * FReal(index) ) ){
	    return index - 1;
	}
	return index;
    }

    struct TestParticle{
	MortonIndex index;
	FSize indexInFile;
	FPoint position;
	FReal physicalValue;
	const FPoint& getPosition()const{
	    return position;
	}
	TestParticle& operator=(const TestParticle& other){
	    index=other.index;
	    indexInFile=other.indexInFile;
	    position=other.position;
	    physicalValue=other.physicalValue;
	    return *this;
	}
	bool operator<(const TestParticle& rhs)const{
	    if(rhs.index < this->index){return false;}
	    else{
		if(rhs.index > this->index){return true;}
		else{
		    return this->indexInFile < rhs.indexInFile;
		}
	    }
	}

    };

    void RunTest(){
	//
	// Load particles
	//
	if(sizeof(FReal) == sizeof(float) ) {
	    std::cerr << "No input data available for Float "<< std::endl;
	    exit(EXIT_FAILURE);
	}
	const std::string parFile( (sizeof(FReal) == sizeof(float))?
				   "Test/DirectFloatbfma":
				   "UTest/DirectDouble.bfma");
	//Let the choice there to test
	std::string filename(SCALFMMDataPath+parFile);
	//std::string filename("../Data/unitCubeXYZQ100.bfma");

	int TreeHeight =3;


	//First part is Sequential
	//Since there is no MPI  Loader with ascii datas, the file need to be in binary format.
	FFmaGenericLoader loaderSeq(filename,true);
	if(!loaderSeq.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;

	//Get the needed informations
	FReal boxWidth = loaderSeq.getBoxWidth();
	FReal boxWidthAtLeafLevel = boxWidth/FReal(1 << (TreeHeight - 1));

	FPoint centerOfBox = loaderSeq.getCenterOfBox();
	FPoint boxCorner   = centerOfBox - boxWidth/2;

	FTreeCoordinate host;

	//Copy from the file to the array that will be sorted
	FSize nbOfParticles = loaderSeq.getNumberOfParticles();
	printf("nbOfParticles : %lld \n",nbOfParticles);
	struct TestParticle* arrayOfParticles = new TestParticle[nbOfParticles];

	memset(arrayOfParticles,0,sizeof(struct TestParticle)*nbOfParticles);

	for(FSize idxParts=0 ; idxParts<nbOfParticles ; ++idxParts){

	    //Fill automatically position AND physicalValue attributes
	    loaderSeq.fillParticle(&(arrayOfParticles[idxParts].position),&(arrayOfParticles[idxParts].physicalValue));

	    //We store the index in the file
	    arrayOfParticles[idxParts].indexInFile = idxParts;

	    //Build temporary TreeCoordinate
	    host.setX( getTreeCoordinate( arrayOfParticles[idxParts].getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
	    host.setY( getTreeCoordinate( arrayOfParticles[idxParts].getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
	    host.setZ( getTreeCoordinate( arrayOfParticles[idxParts].getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));

	    //Set Morton index from Tree Coordinate
	    arrayOfParticles[idxParts].index = host.getMortonIndex(TreeHeight - 1);
	}
	//Save the original array
	struct TestParticle * originalArray   =  new TestParticle[nbOfParticles];
	memcpy(originalArray,arrayOfParticles,sizeof(struct TestParticle)*nbOfParticles);
	//Sort the array

	std::sort(arrayOfParticles,arrayOfParticles+nbOfParticles);
	// for(int k=0 ; k< nbOfParticles ; ++k){
	//   printf("arrayOfParticles[].index %lld \n",arrayOfParticles[k].index);
	// }
	//Start of the parallel part :
	MortonIndex ref = -1;
	int numMort = 0;
	if(app.global().processId()==0){
	    for(int i=0 ; i<loaderSeq.getNumberOfParticles() ; ++i){
		if (arrayOfParticles[i].index !=ref){
		    numMort++;
		    ref = arrayOfParticles[i].index;
		}
	    }
	    printf("Total leaf : %d Last : %lld \n",numMort,ref);
	}

	FSize outputSize;

	//Refer to ChebyshevInterpolationAlgorithmProc to know how to FMpiFmaLoader + index
	FMpiFmaGenericLoader loader(filename,app.global());
	if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;
	//Now, we sort again the particles with MPI QuickSort
	int idxStart = loader.getStart();

	FMpiTreeBuilder<TestParticle>::IndexedParticle * arrayToBeSorted = new FMpiTreeBuilder<TestParticle>::IndexedParticle[loader.getMyNumberOfParticles()];
	//Copy the TestParticles into an array of indexedParticle
	for(int i=0 ; i<loader.getMyNumberOfParticles() ; ++i){
	    arrayToBeSorted[i].particle = originalArray[i+idxStart];
	    arrayToBeSorted[i].index = arrayToBeSorted[i].particle.index;
	}
	FMpiTreeBuilder<TestParticle>::IndexedParticle* outputArray = nullptr;
	FQuickSortMpi<FMpiTreeBuilder<TestParticle>::IndexedParticle,MortonIndex,FSize>::QsMpi(arrayToBeSorted,loader.getMyNumberOfParticles(),&outputArray,&outputSize,app.global());

	//FBitonicSort<FMpiTreeBuilder<TestParticle>::IndexedParticle,MortonIndex, FSize>::Sort(arrayToBeSorted,loader.getMyNumberOfParticles(),app.global());
	//Sum the outputSize of every body for knowing where to start inside the sorted array
	FSize starter = 0;
	//We use a prefix sum


	MPI_Exscan(&outputSize,&starter,1,MPI_LONG_LONG_INT,MPI_SUM,app.global().getComm());

	//We sort the output array relatvely to line number in origin file
	FSize inc = 0;
	FMpiTreeBuilder<TestParticle>::IndexedParticle * saveForSort = outputArray;
	int nbOfPartsInLeaf = 0;
	while(inc < outputSize){
	    while(outputArray[inc].index == saveForSort->index && inc < outputSize){
		inc++;
		nbOfPartsInLeaf++;
	    }
	    std::sort(saveForSort,saveForSort+nbOfPartsInLeaf,
		      [&](FMpiTreeBuilder<TestParticle>::IndexedParticle a,FMpiTreeBuilder<TestParticle>::IndexedParticle b) -> bool {
			  return (a.particle.indexInFile)<(b.particle.indexInFile);
		      });
	    nbOfPartsInLeaf = 0;
	    saveForSort = &outputArray[inc];
	}

	bool resultQsMpi = true; //passed..
	//Test

	//Gather size of output
	int * nbPartPerProcess = new int[app.global().processCount()];
	nbPartPerProcess[app.global().processId()] = outputSize;
	MPI_Gather(&nbPartPerProcess[app.global().processId()],1,MPI_INT,nbPartPerProcess,1,MPI_INT,0,app.global().getComm());

	int * toSend = new int[outputSize];
	int * displ = nullptr;
	int * recvParts = nullptr;

	//Prepare the indexInFile to send
	for(int idPart=0 ; idPart<outputSize ; ++idPart){
	    toSend[idPart] = int(outputArray[idPart].particle.indexInFile);
	}


	if(app.global().processId() == 0){
	    //There, we build the array of displacement
	    displ = new int[app.global().processCount()];
	    displ[0] = 0;
	    for(int idProc = 1 ; idProc < app.global().processCount() ; ++idProc){
		displ[idProc] = nbPartPerProcess[idProc-1] + displ[idProc-1];
	    }
	    //Buffer to recv into
	    recvParts = new int[loader.getNumberOfParticles()];
	    MPI_Gatherv(toSend,outputSize,MPI_INT,recvParts,nbPartPerProcess,displ,MPI_INT,0,app.global().getComm());

	    //Buffer to put result into
	    int * myPart = new int[loader.getNumberOfParticles()];
	    memset(myPart,0,sizeof(int)*loader.getNumberOfParticles());
	    for(int idP = 0 ; idP < loader.getNumberOfParticles() ; ++idP){
		myPart[recvParts[idP]] += 1;
	    }
	    //Check if everything is set to 1
	    for(int idP = 0 ; idP < loader.getNumberOfParticles() ; ++idP){
		if(myPart[idP] != 1){
		    std::cout << "Part number "<< idP << " in file is lost or duplicated : "<< myPart[idP]<< std::endl;
		    resultQsMpi = false;
		}
	    }
	}
	else{
	    MPI_Gatherv(toSend,outputSize,MPI_INT,recvParts,nbPartPerProcess,displ,MPI_INT,0,app.global().getComm());
	}



	//Test
	// printf("Fin tri : %d hold %d my first [%lld,%lld] last [%lld,%lld]\n",app.global().processId(),outputSize,
	//        outputArray[0].index,outputArray[0].particle.indexInFile,
	//        outputArray[outputSize-1].index,outputArray[outputSize-1].particle.indexInFile);
	// int tot = 0;
	// MPI_Reduce((int*)&outputSize,&tot,1,MPI_INT,MPI_SUM,0,app.global().getComm());
	// if(app.global().processId() == 0){
	//     printf("Total Particule après tri = %d \n",tot);
	// }
	// for(int i=0 ; i<outputSize ; ++i){
	//     if(outputArray[i].particle.indexInFile != arrayOfParticles[i+starter].indexInFile){
	//	resultQsMpi = false;
	//	printf("QsMpi :: %d finalParticle  [%lld,%lld] : sortedArray [%lld,%lld] \n",app.global().processId(),
	//	       outputArray[i].particle.indexInFile, outputArray[i].index,
	//	       arrayOfParticles[i+starter].indexInFile,arrayOfParticles[i+starter].index);
	//     }
	// }


	Print("Test 1 : is QsMpi really sorting the array");
	uassert(resultQsMpi);


	//Test MergeLeaves
	bool resultMergeLeaves= true;

	//inputs needed
	TestParticle * leavesArray = nullptr;
	FSize * leavesIndices = nullptr;
	FSize leaveSize = 0;

	FMpiTreeBuilder<TestParticle>::MergeSplitedLeaves(app.global(),outputArray,&outputSize,&leavesIndices,&leavesArray,&leaveSize);

	//Compare again the results with the output of std::qsort

	//we need to know how many parts still remains
	FSize CounterStart = 0;
	//We use a prefix sum
	MPI_Exscan(&outputSize,&CounterStart,1,MPI_LONG_LONG_INT,MPI_SUM,app.global().getComm());

	//Test if no problems
	for(FSize k=0 ; k<outputSize ; ++k){
	    if(leavesArray[k].indexInFile != arrayOfParticles[k+CounterStart].indexInFile){
		printf("MergeLeaves :: Proc %d, finalParticles : %lld,%lld, sortedArray %lld,%lld \n",
		       app.global().processId(),
		       leavesArray[k].index,leavesArray[k].indexInFile,
		       arrayOfParticles[k+CounterStart].index,arrayOfParticles[k+CounterStart].indexInFile);
		resultMergeLeaves = false;
	    }
	}
	Print("Test 2 : Output of merging leaves is tested");
	uassert(resultMergeLeaves);

	//Test the Equalize and Fill tree
	FLeafBalance balancer;
	FVector<TestParticle> finalParticles;

	bool resultEqualize = true;

	FMpiTreeBuilder<TestParticle>::EqualizeAndFillContainer(app.global(),&finalParticles,leavesIndices,leavesArray,leaveSize,outputSize,&balancer);
	//Ok now count the Particles at the end of the Equalize
	int finalNbPart = finalParticles.getSize();
	int finalStart = 0;

	MPI_Exscan(&finalNbPart,&finalStart,1,MPI_INT,MPI_SUM,app.global().getComm());
	for (int k=0; k<finalNbPart ; k++){
	    if(finalParticles[k].indexInFile != arrayOfParticles[k+finalStart].indexInFile){
		printf("Equalize :: Proc %d, k=[%d+%d] finalParticles : %lld,%lld, sortedArray %lld,%lld \n",
		       app.global().processId(),k,finalStart,
		       finalParticles[k].index,finalParticles[k].indexInFile,
		       arrayOfParticles[k+finalStart].index,arrayOfParticles[k+finalStart].indexInFile);
		resultEqualize = false;
	    }
	}

	Print("Test 3 : Output of Equalize Tree is tested");
	MPI_Barrier(MPI_COMM_WORLD);
	uassert(resultEqualize);

	delete [] originalArray;
	delete [] arrayOfParticles;
	delete [] arrayToBeSorted;
	delete [] outputArray;

    }
    /** If memstas is running print the memory used */
    void PostTest() {
	if( FMemStats::controler.isUsed() ){
	    std::cout << app.global().processId() << "-> Memory used at the end " << FMemStats::controler.getCurrentAllocated()
		      << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
	    std::cout << app.global().processId() << "-> Max memory used " << FMemStats::controler.getMaxAllocated()
		      << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
	    std::cout << app.global().processId() << "-> Total memory used " << FMemStats::controler.getTotalAllocated()
		      << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
	}
    }

    void SetTests(){
	AddTest(&TestMpiTreeBuilder::RunTest,"Load a File, sort it, merge it, and Equalize it (4 steps)");
    }

public:
    TestMpiTreeBuilder(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }


};

TestClassMpi(TestMpiTreeBuilder);
