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
      if(rhs.index < this->index){
	return false;
      }
      else{
	if(rhs.index > this->index){
	  return true;
	}
	else{
	  if(rhs.getPosition().getX() < this->getPosition().getX()){
	    return false;
	  }
	  else{
	    if(rhs.getPosition().getX() > this->getPosition().getX()){
	      return true;
	    }
	    else{
	      if(rhs.getPosition().getY() < this->getPosition().getY()){
		return false;
	      }
	      else{
		if(rhs.getPosition().getY() > this->getPosition().getY()){
		  return true;
		  }
		else{
		  if(rhs.getPosition().getZ() < this->getPosition().getZ()){
		    return false;
		  }
		  else{
		    if(rhs.getPosition().getZ() > this->getPosition().getZ()){
		      return true;
		    }
		    else{
		      return false;
		    }
		  }
		}
	      }
	    }
	  }
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
    // const std::string parFile( (sizeof(FReal) == sizeof(float))?
    // 			       "Test/DirectFloatbfma":
    // 			       "UTest/DirectDouble.bfma");
    // //
    // std::string filename(SCALFMMDataPath+parFile);
    std::string filename("../Data/unitCubeXYZQ100.bfma");
    int TreeHeight =4;
    
    
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
    FQuickSortMpi<FMpiTreeBuilder<TestParticle>::IndexedParticle,MortonIndex,FSize>::QsMpi(arrayToBeSorted,loader.getMyNumberOfParticles(),outputArray,outputSize,app.global());
    //Sum the outputSize of every body for knowing where to start inside the sorted array
    FSize starter = 0; 
    //We use a prefix sum
    MPI_Exscan(&outputSize,&starter,1,MPI_LONG_LONG_INT,MPI_SUM,app.global().getComm());
    
    //We sort the output array relatvely to position inside leafs
    FSize inc = 0;
    FMpiTreeBuilder<TestParticle>::IndexedParticle * saveForSort = outputArray;
    int nbOfPartsInLeaf = 0;
    while(inc < outputSize){
      while(outputArray[inc].index == saveForSort->index){
	inc++;
	nbOfPartsInLeaf++;
      }
      std::sort(saveForSort,saveForSort+nbOfPartsInLeaf,
		[&](FMpiTreeBuilder<TestParticle>::IndexedParticle a,FMpiTreeBuilder<TestParticle>::IndexedParticle b) -> bool {
		  return (a.particle)<(b.particle);
		});
      nbOfPartsInLeaf = 0;
      saveForSort = &outputArray[inc];
    }
  
    bool resultQsMpi = true; //passed..
    //Test
    for(int i=0 ; i<outputSize ; ++i){
      if(outputArray[i].particle.indexInFile != arrayOfParticles[i+starter].indexInFile){
	printf("QsMPI :: Proc %d i=%d Particles file : [%lld,%lld] with Mpi Morton : %lld != %lld\n",app.global().processId(),i,
	       arrayOfParticles[i+idxStart].indexInFile,
	       outputArray[i].particle.indexInFile,
	       outputArray[i].index,
	       arrayOfParticles[i+idxStart].index);
	resultQsMpi = false;
      }
    }
    Print("Test 1 : is QsMpi really sorting the array");
    uassert(resultQsMpi);
    

    //Test MergeLeaves
    bool resultMergeLeaves= true;
    
    //inputs needed
    TestParticle * leavesArray = nullptr;
    FSize * leavesIndices = nullptr;
    FSize leaveSize = 0;
    
    FMpiTreeBuilder<TestParticle>::testMergeLeaves(app.global(),outputArray,&outputSize,&leavesIndices,&leavesArray,&leaveSize);
    
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

    FMpiTreeBuilder<TestParticle>::testEqualizeAndFillTree(app.global(),&finalParticles,leavesIndices,leavesArray,leaveSize,outputSize,&balancer);
    //Ok now count the Particles at the end of the Equalize
    int finalNbPart = finalParticles.getSize();
    int finalStart = 0;
    MPI_Exscan(&finalNbPart,&finalStart,1,MPI_INT,MPI_SUM,app.global().getComm());
    
    for (int k=0; k<finalNbPart ; k++){
      if(finalParticles[k].indexInFile != arrayOfParticles[k+finalStart].indexInFile){
    	printf("Equalize :: Proc %d, finalParticles : %lld,%lld, sortedArray %lld,%lld \n",
    	       app.global().processId(),
    	       finalParticles[k].index,finalParticles[k].indexInFile,
    	       arrayOfParticles[k+finalStart].index,arrayOfParticles[k+finalStart].indexInFile);
    	resultEqualize = false;
      }
    }
    
    Print("Test 3 : Output of Equalize Tree is tested");
    uassert(resultEqualize);
  
    delete [] originalArray;
    delete [] arrayOfParticles;
    delete [] arrayToBeSorted;
    delete [] outputArray;

  }
  
  void SetTests(){
    AddTest(&TestMpiTreeBuilder::RunTest,"Load a File, sort it, merge it, and Equalize it (4 steps)");
  }
  
public:
  TestMpiTreeBuilder(int argc,char ** argv) : FUTesterMpi(argc,argv){
  }


};

TestClassMpi(TestMpiTreeBuilder);
