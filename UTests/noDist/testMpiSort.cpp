#include <cstdlib>
#include <string.h>
#include <stdexcept>

#include "../../Src/Utils/FMpi.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiFmaLoader.hpp"
#include "../../Src/BalanceTree/FLeafBalance.hpp"
#include "../../Src/Containers/FTreeCoordinate.hpp"


#include "../../Src/Utils/FQuickSortMpi.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"
#include "../../Src/Core/FCoreCommon.hpp"

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Utils/FMath.hpp"


#include "../../Src/Utils/FParameters.hpp"

// ==== CMAKE =====
// @FUSE_MPI
// ================


/* Struct Test particle : classic particle but with index In file,
 * which is the number of line where the part was found. It is a true
 * constant.
 */
struct TestParticle{
  MortonIndex index;
  FSize indexInFile;
  FPoint position;
  FReal physicalValue;
  FPoint& getPosition(){
    return position;
  }
};



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

/**
 * This function compare IndexedParticles that are in the same leaf by
 * sorting them along X,Y and then Z axis. The aim is to get the same
 * sorting output for each leaf.
 */
int compareLeafI(const void *a,const void *b){
  FMpiTreeBuilder<TestParticle>::IndexedParticle partA = *(reinterpret_cast<const FMpiTreeBuilder<TestParticle>::IndexedParticle*>(a));
  FMpiTreeBuilder<TestParticle>::IndexedParticle partB = *(reinterpret_cast<const FMpiTreeBuilder<TestParticle>::IndexedParticle*>(b));
  if(partA.particle.getPosition().getX() < partB.particle.getPosition().getX()){
    return -1;
  }
  else{
    if(partA.particle.getPosition().getX() > partB.particle.getPosition().getX()){
      return 1;
    }
    else{
      if(partA.particle.getPosition().getY() < partB.particle.getPosition().getY()){
	return -1;
      }
      else{
	if(partA.particle.getPosition().getY() > partB.particle.getPosition().getY()){
	  return 1;
	}
	else{
	  if(partA.particle.getPosition().getZ() < partB.particle.getPosition().getZ()){
	    return -1;
	  }
	  else{
	    if(partA.particle.getPosition().getZ() > partB.particle.getPosition().getZ()){
	      return 1;
	    }
	    else{
	      return 0;
	    }
	  }
	}
      }
    }
  }
}

/**
 * Same function as above but with TestParticle instead of
 * IndexedParticles (pretty much same thing).
 */
int compareLeaf(const void *a,const void *b){
  TestParticle partA = *(reinterpret_cast<const TestParticle*>(a));
  TestParticle partB = *(reinterpret_cast<const TestParticle*>(b));
  if(partA.getPosition().getX() < partB.getPosition().getX()){
    return -1;
  }
  else{
    if(partA.getPosition().getX() > partB.getPosition().getX()){
      return 1;
    }
    else{
      if(partA.getPosition().getY() < partB.getPosition().getY()){
	return -1;
      }
      else{
	if(partA.getPosition().getY() > partB.getPosition().getY()){
	  return 1;
	}
	else{
	  if(partA.getPosition().getZ() < partB.getPosition().getZ()){
	    return -1;
	  }
	  else{
	    if(partA.getPosition().getZ() > partB.getPosition().getZ()){
	      return 1;
	    }
	    else{
	      return 0;
	      printf("Should not happen\n");
	    }
	  }
	}
      }
    }
  }
}

/**
 * Compare two articles by MortonIndex, if same morton index, call to
 * compareLeaf
 */
int cmpTestParticles(const void* a,const void* b){

  const TestParticle partA = *(reinterpret_cast<const TestParticle*>(a));
  const TestParticle partB = *(reinterpret_cast<const TestParticle*>(b));

  if(partA.index < partB.index){
    return -1;
  }
  else{
    if(partA.index > partB.index){
      return 1;
    }
    else{
      return compareLeaf(a,b);
    }
  }
}


/** This test aimed to validate the Tree Builder algorithms 
 *
 *  Each process runs a sort (from std::qsort) and then a parallel
 *  sort and compare its particles with the global sort.
 *
 * Then, each process runs a MergeLeaves, and compare again the output
 * with the one from sstd::qsort.
 *
 * Finally, each process runs EqualizeAndFillTree, and compare again
 * the output with the one from std::sort.
 *
 */
int main(int argc,char* argv[]){
  //File to sort
  const char* const filename          = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma.double");
  const unsigned int TreeHeight       = FParameters::getValue(argc, argv, "-depth", 5);
  

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
  struct TestParticle* arrayOfParticles = new struct TestParticle[nbOfParticles];
  
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
  qsort((void*) arrayOfParticles,nbOfParticles,sizeof(TestParticle),cmpTestParticles);
  

  //Start of the parallel part : 
  FMpi app(argc,argv);
  
  FSize outputSize;

  //Refer to ChebyshevInterpolationAlgorithmProc to know how to FMpiFmaLoader + index
  FMpiFmaLoader loader(filename,app.global());
  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;
  

  //Now, we sort again the particles with MPI QuickSort 
  int idxStart = loader.getStart();

  FMpiTreeBuilder<TestParticle>::IndexedParticle * arrayToBeSorted = new FMpiTreeBuilder<TestParticle>::IndexedParticle[loader.getNumberOfParticles()];
  //Copy the TestParticles into an array of indexedParticle
  for(int i=0 ; i<loader.getNumberOfParticles() ; ++i){
    arrayToBeSorted[i].particle = originalArray[i+idxStart];
    arrayToBeSorted[i].index = arrayToBeSorted[i].particle.index;
  }
  
  FMpiTreeBuilder<TestParticle>::IndexedParticle* outputArray = 0;
  FQuickSortMpi<FMpiTreeBuilder<TestParticle>::IndexedParticle,MortonIndex,FSize>::QsMpi(arrayToBeSorted,loader.getNumberOfParticles(),outputArray,outputSize,app.global());
  
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
    qsort((void*)saveForSort,nbOfPartsInLeaf,sizeof(FMpiTreeBuilder<TestParticle>::IndexedParticle),compareLeafI);
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
  if(!resultQsMpi){
    printf("Proc :: %d ====> QsMpi didn't make it !!\n",app.global().processId());
  }
  else{
    printf("Proc :: %d ====> QsMpi successfull !!\n",app.global().processId());
  }
  
  //Test MergeLeaves
  bool resultMergeLeaves= true;
  
  //inputs needed
  TestParticle * leavesArray = 0;
  FSize * leavesIndices = 0;
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
  
  if(!resultMergeLeaves){
    printf("Proc :: %d ====> MergeLeaves didn't make it !!\n",app.global().processId());
  }
  else{
    printf("Proc :: %d ====> MergeLeaves successfull !!\n",app.global().processId());
  }
    

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

  if(!resultEqualize){
    printf("Proc :: %d ====> Equalize didn't make it !!\n",app.global().processId());
  }
  else{
    printf("Proc :: %d ====> Equalize successfull !!\n",app.global().processId());
  }
    
 
  delete [] originalArray;
  delete [] arrayOfParticles;
  delete [] arrayToBeSorted;
  delete [] outputArray;

  return EXIT_SUCCESS;
}
