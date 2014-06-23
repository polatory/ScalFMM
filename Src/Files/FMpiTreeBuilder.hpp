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
#ifndef FMPITREEBUILDER_H
#define FMPITREEBUILDER_H

#include "../Utils/FMpi.hpp"
#include "../Utils/FQuickSortMpi.hpp"
#include "../Utils/FBitonicSort.hpp"

#include "../Utils/FMemUtils.hpp"
#include "../Utils/FTrace.hpp"

#include "../BalanceTree/FLeafBalance.hpp"
#include "../BalanceTree/FEqualize.hpp"

/** This class manage the loading of particles for the mpi version.
 * It use a BinLoader and then sort the data with a parallel quick sort
 * or bitonic sort.
 * After it needs to merge the leaves to finally equalize the data.
 */
template<class ParticleClass>
class FMpiTreeBuilder{
public:
  /** What sorting algorithm to use */
  enum SortingType{
    QuickSort,
    BitonicSort,
  };


  /** A particle may not have a MortonIndex Method
   * but they are sorted on their morton index
   * so this struct store a particle + its index
   */
  struct IndexedParticle{
  public:
    MortonIndex index;
    ParticleClass particle;

    operator MortonIndex(){
      return this->index;
    }
  };

  
private:
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


  //////////////////////////////////////////////////////////////////////////
  // To sort tha particles we hold
  //////////////////////////////////////////////////////////////////////////


  template <class LoaderClass>
  static void SortParticles( const FMpi::FComm& communicator, LoaderClass& loader, const SortingType type,
			     const int TreeHeight, IndexedParticle**const outputArray, FSize* const outputSize){
    FTRACE( FTrace::FFunction functionTrace(__FUNCTION__ , "Loader to Tree" , __FILE__ , __LINE__) );

    // create particles
    IndexedParticle*const realParticlesIndexed = new IndexedParticle[loader.getNumberOfParticles()];
    FMemUtils::memset(realParticlesIndexed, 0, sizeof(IndexedParticle) * loader.getNumberOfParticles());
    FPoint boxCorner(loader.getCenterOfBox() - (loader.getBoxWidth()/2));
    FTreeCoordinate host;

    const FReal boxWidthAtLeafLevel = loader.getBoxWidth() / FReal(1 << (TreeHeight - 1) );
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
      loader.fillParticle(realParticlesIndexed[idxPart].particle);
      host.setX( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
      host.setY( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
      host.setZ( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));

      realParticlesIndexed[idxPart].index = host.getMortonIndex(TreeHeight - 1);
    }

    // sort particles
    if(type == QuickSort){
      FQuickSortMpi<IndexedParticle,MortonIndex, FSize>::QsMpi(realParticlesIndexed, loader.getNumberOfParticles(), *outputArray, *outputSize,communicator);
      delete [] (realParticlesIndexed);
    }
    else {
      FBitonicSort<IndexedParticle,MortonIndex, FSize>::Sort( realParticlesIndexed, loader.getNumberOfParticles(), communicator );
      *outputArray = realParticlesIndexed;
      *outputSize = loader.getNumberOfParticles();
    }
  }

  static void SortParticlesFromArray( const FMpi::FComm& communicator, const ParticleClass array[], const FSize size, const SortingType type,
				      const FPoint& centerOfBox, const FReal boxWidth,
				      const int TreeHeight, IndexedParticle**const outputArray, FSize* const outputSize){
    FTRACE( FTrace::FFunction functionTrace(__FUNCTION__ , "Loader to Tree" , __FILE__ , __LINE__) );

    // create particles
    IndexedParticle*const realParticlesIndexed = new IndexedParticle[size];
    FMemUtils::memset(realParticlesIndexed, 0, sizeof(IndexedParticle) * size);

    FPoint boxCorner(centerOfBox - (boxWidth/2));
    FTreeCoordinate host;

    const FReal boxWidthAtLeafLevel = boxWidth / FReal(1 << (TreeHeight - 1) );

    for(int idxPart = 0 ; idxPart < size ; ++idxPart){
      realParticlesIndexed[idxPart].particle = array[idxPart];
      host.setX( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidthAtLeafLevel ));
      host.setY( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidthAtLeafLevel ));
      host.setZ( getTreeCoordinate( realParticlesIndexed[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidthAtLeafLevel ));

      realParticlesIndexed[idxPart].index = host.getMortonIndex(TreeHeight - 1);
    }

    // sort particles
    if(type == QuickSort){
      FQuickSortMpi<IndexedParticle,MortonIndex, FSize>::QsMpi(realParticlesIndexed, size, *outputArray, *outputSize,communicator);
      delete [] (realParticlesIndexed);
    }
    else {
      FBitonicSort<IndexedParticle,MortonIndex, FSize>::Sort( realParticlesIndexed, size, communicator );
      *outputArray = realParticlesIndexed;
      *outputSize = size;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // To merge the leaves
  //////////////////////////////////////////////////////////////////////////

  static void MergeLeaves(const FMpi::FComm& communicator, IndexedParticle*& workingArray, FSize* workingSize,
			  FSize ** leavesIndices, ParticleClass** leavesArray, FSize* const leavesSize){
    FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader to Tree" , __FILE__ , __LINE__) );
    const int rank = communicator.processId();
    const int nbProcs = communicator.processCount();
    if(nbProcs == 1){
      //Nothing to do there : there is no need to verify if leaves are split, if there is one process...
    }
    else{
      // be sure there is no splited leaves
      // to do that we exchange the first index with the left proc
      {
	FTRACE( FTrace::FRegion regionTrace("Remove Splited leaves", __FUNCTION__ , __FILE__ , __LINE__) );

	MortonIndex otherFirstIndex = -1;
	if((*workingSize) != 0 && rank != 0 && rank != nbProcs - 1){
	  MPI_Sendrecv(&workingArray[0].index, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs,
		       &otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs,
		       communicator.getComm(), MPI_STATUS_IGNORE);
	}
	else if( rank == 0){
	  MPI_Recv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs, communicator.getComm(), MPI_STATUS_IGNORE);
	}
	else if( rank == nbProcs - 1){
	  MPI_Send( &workingArray[0].index, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs, communicator.getComm());
	}
	else {
	  MPI_Recv(&otherFirstIndex, 1, MPI_LONG_LONG, rank + 1, FMpi::TagExchangeIndexs, communicator.getComm(), MPI_STATUS_IGNORE);
	  MPI_Send(&otherFirstIndex, 1, MPI_LONG_LONG, rank - 1, FMpi::TagExchangeIndexs, communicator.getComm());
	}

	// at this point every one know the first index of his right neighbors
	const bool needToRecvBeforeSend = (rank != 0 && (((*workingSize) && otherFirstIndex == workingArray[0].index ) || !(*workingSize)));
	MPI_Request requestSendLeaf;

	IndexedParticle* sendBuffer = nullptr;
	if(rank != nbProcs - 1 && needToRecvBeforeSend == false){
	  FSize idxPart = (*workingSize) - 1 ;
	  while(idxPart >= 0 && workingArray[idxPart].index == otherFirstIndex){
	    --idxPart;
	  }
	  const int particlesToSend = int((*workingSize) - 1 - idxPart);
	  if(particlesToSend){
	    (*workingSize) -= particlesToSend;
	    sendBuffer = new IndexedParticle[particlesToSend];
	    memcpy(sendBuffer, &workingArray[idxPart + 1], particlesToSend * sizeof(IndexedParticle));

	    MPI_Isend( sendBuffer, particlesToSend * int(sizeof(IndexedParticle)), MPI_BYTE,
		       rank + 1, FMpi::TagSplittedLeaf, communicator.getComm(), &requestSendLeaf);
	  }
	  else{
	    MPI_Isend( nullptr, 0, MPI_BYTE, rank + 1, FMpi::TagSplittedLeaf, communicator.getComm(), &requestSendLeaf);
	  }
	}

	if( rank != 0 ){
	  int sendByOther = 0;

	  MPI_Status probStatus;
	  MPI_Probe(rank - 1, FMpi::TagSplittedLeaf, communicator.getComm(), &probStatus);
	  MPI_Get_count( &probStatus,  MPI_BYTE, &sendByOther);

	  if(sendByOther){
	    sendByOther /= int(sizeof(IndexedParticle));

	    const IndexedParticle* const reallocOutputArray = workingArray;
	    const FSize reallocOutputSize = (*workingSize);

	    (*workingSize) += sendByOther;
	    workingArray = new IndexedParticle[(*workingSize)];
	    FMemUtils::memcpy(&workingArray[sendByOther], reallocOutputArray, reallocOutputSize * sizeof(IndexedParticle));
	    delete[] reallocOutputArray;

	    MPI_Recv(workingArray, int(sizeof(IndexedParticle)) * sendByOther, MPI_BYTE,
		     rank - 1, FMpi::TagSplittedLeaf, communicator.getComm(), MPI_STATUS_IGNORE);
	  }
	  else{
	    MPI_Recv( nullptr, 0, MPI_BYTE, rank - 1, FMpi::TagSplittedLeaf, communicator.getComm(), MPI_STATUS_IGNORE);
	  }
	}

	if(rank != nbProcs - 1 && needToRecvBeforeSend == true){
	  MPI_Send( workingArray, int((*workingSize) * sizeof(IndexedParticle)), MPI_BYTE,
		    rank + 1, FMpi::TagSplittedLeaf, communicator.getComm());
	  delete[] workingArray;
	  workingArray = nullptr;
	  (*workingSize)  = 0;
	}
	else if(rank != nbProcs - 1){
	  MPI_Wait( &requestSendLeaf, MPI_STATUS_IGNORE);
	  delete[] sendBuffer;
	  sendBuffer = nullptr;
	}
      }
    }
    {//Filling the Array with leaves and parts //// COULD BE MOVED IN AN OTHER FUCTION
	  
      (*leavesSize)    = 0; //init ptr
      (*leavesArray)   = nullptr; //init ptr
      (*leavesIndices) = nullptr; //init ptr
	  
      if((*workingSize)){
	    
	//Array of particles
	*leavesArray = new ParticleClass[(*workingSize)];
	    
	//Temporary array, we don't know yet how many leaves there will be 
	FSize * tempIndicesArray = new FSize[(*workingSize)];
	memset(tempIndicesArray,0,sizeof(FSize)*(*workingSize));
	    
	FSize idxInIndices = 0;
	MortonIndex previousIndex = -1;
	for(FSize idxPart = 0 ; idxPart < (*workingSize) ; ++idxPart){
	  if(workingArray[idxPart].index != previousIndex){
	    previousIndex = workingArray[idxPart];
	    tempIndicesArray[idxInIndices] = idxPart;
	    idxInIndices++;
	  }
	  memcpy(&(*leavesArray)[idxPart],&workingArray[idxPart].particle,sizeof(ParticleClass));
	}
	*leavesSize = idxInIndices;
	    
	*leavesIndices = new FSize[idxInIndices];
	memcpy(*leavesIndices,tempIndicesArray,(*leavesSize)*sizeof(FSize));

      }
	  
      delete []  workingArray;
	  
      workingArray = nullptr;
	
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // To equalize (same number of leaves among the procs)
  //////////////////////////////////////////////////////////////////////////

  /** Put the interval into a tree */
  template <class ContainerClass>
  static void EqualizeAndFillTree(const FMpi::FComm& communicator,  ContainerClass* particlesSaver,
				  FSize * leavesIndices, ParticleClass* leavesArray, FSize& nbLeavesInIntervals, FSize nbParts, 
				  FAbstractBalanceAlgorithm * balancer){


    FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Loader to Tree" , __FILE__ , __LINE__) );
    const int myRank = communicator.processId();
    const int nbProcs = communicator.processCount();
    const FSize myCurrentsParts = nbParts;
	
    if(nbProcs == 1){ //I'm the only one !
      //Just copy each part into the Particle Saver
      for(FSize idxPart =0 ; idxPart < myCurrentsParts ; ++idxPart){
	particlesSaver->push(leavesArray[idxPart]);
      }
    }
    else{
            // We have to know the number of leaves each procs holds
            FSize*const numberOfLeavesPerProc = new FSize[nbProcs];
            memset(numberOfLeavesPerProc, 0, sizeof(FSize) * nbProcs);
            FSize intNbLeavesInIntervals = nbLeavesInIntervals;
            MPI_Allgather(&intNbLeavesInIntervals, 1, MPI_LONG_LONG_INT, numberOfLeavesPerProc, 1, MPI_LONG_LONG_INT, communicator.getComm());

            //We need the max number of leafs over th procs
            FSize*const diffNumberOfLeavesPerProc = new FSize[nbProcs+1];
            diffNumberOfLeavesPerProc[0] = 0;
            FSize totalNumberOfLeaves  = 0;
            for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc ){
                totalNumberOfLeaves += numberOfLeavesPerProc[idxProc];
                diffNumberOfLeavesPerProc[idxProc+1] = diffNumberOfLeavesPerProc[idxProc] + numberOfLeavesPerProc[idxProc];
            }

            //Building of counter send buffer
            FSize * toSend = new FSize[nbProcs];
            memset(toSend,0,sizeof(FSize)*nbProcs);

            //Building of array of indices to send
            FSize * idxToSend = new FSize[nbProcs];
            memset(idxToSend,0,sizeof(FSize)*nbProcs);

            std::vector< std::pair<size_t,size_t> > allObjectives;
            allObjectives.resize(nbProcs);
            for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
                allObjectives[idxProc].first  = balancer->getLeft(totalNumberOfLeaves,NULL,0,0,nbProcs,idxProc);
                allObjectives[idxProc].second = balancer->getRight(totalNumberOfLeaves,NULL,0,0,nbProcs,idxProc);
            }

            std::pair<size_t, size_t> myCurrentInter = {diffNumberOfLeavesPerProc[myRank], diffNumberOfLeavesPerProc[myRank+1]};
            const std::vector<FEqualize::Package> packsToSend = FEqualize::GetPackToSend(myCurrentInter, allObjectives);

            for(const FEqualize::Package& pack : packsToSend){
                printf("%d] to %d from %llu to %llu\n", myRank, pack.idProc, pack.elementFrom, pack.elementTo);
                idxToSend[pack.idProc] = pack.elementFrom;
                toSend[pack.idProc]    = (pack.elementTo ? leavesIndices[pack.elementTo-1] : 0)
                                        - (pack.elementFrom ? leavesIndices[pack.elementFrom-1] : 0);
            }

            //Then, we exchange the datas to send
            FSize * globalSendRecvMap = new FSize[nbProcs*nbProcs];
            memset(globalSendRecvMap,0,sizeof(FSize)*nbProcs*nbProcs);
            //This could be replace by an array toRecv buildt in the same way as toSend
            MPI_Allgather(toSend,nbProcs,MPI_LONG_LONG,globalSendRecvMap,nbProcs,MPI_LONG_LONG,communicator.getComm());

            //Then, we have our global recv map.
            //We just need to send and recv for real.

            //Finally, store the remaining parts, recv the parts, send my parts
            ParticleClass * finalPartBuffer;
            FSize finalTotParts;
            {
                finalTotParts = myCurrentsParts;     //We need to know how many particles we will have
                FSize finalCurrentParts = myCurrentsParts; //Parts that I had and that belongs to me

                for(int idxProc=0 ; idxProc<nbProcs ; ++idxProc){
                    finalCurrentParts -= toSend[idxProc]; //substract the parts sent
                    finalTotParts -= toSend[idxProc];     //substract the parts sent
                    finalTotParts += globalSendRecvMap[idxProc*nbProcs+myRank]; //add the parts received
                }
                finalPartBuffer = new ParticleClass[finalTotParts];
                memset(finalPartBuffer,0,sizeof(ParticleClass)*finalTotParts);

                //Copy of the parts we already hold
                FSize finalIdxToStart = 0; //idx of the start of my parts inside leavesArray
                FSize idxToWrite = 0;      //idx to write my parts
                {//we go from idxProc==0 to idxProc==myRank to increment the self starter
                    for(int idxProc=0 ; idxProc<myRank ; ++idxProc){
                        idxToWrite += globalSendRecvMap[idxProc*nbProcs+myRank];
                        finalIdxToStart += toSend[idxProc];
                    }
                    memcpy(&finalPartBuffer[idxToWrite],&leavesArray[finalIdxToStart],sizeof(ParticleClass)*finalCurrentParts);
                }


                //Second, receive in place:

                MPI_Request* requests = new MPI_Request[nbProcs * 2];
                int counterRequest = 0;
                int tag = 99;
                FSize idxToWriteRecvedDatas = 0;
                //While I received from left, i write the datas at the start of the buffer
                for(int idxProc=0 ; idxProc<nbProcs ; ++idxProc){
                    if(idxProc == myRank){ //When I myRank==idxProc, I increment the idxToWrite of What I kept to avoid erasing my parts with received parts
                        idxToWriteRecvedDatas += finalCurrentParts;
                    }
                    else{ //I received and inc the write index of what I get
                        if(globalSendRecvMap[idxProc*nbProcs+myRank]){//If i expect something from idxProc
                            MPI_Irecv(&finalPartBuffer[idxToWriteRecvedDatas],int(sizeof(ParticleClass))*int(globalSendRecvMap[idxProc*nbProcs+myRank]),MPI_BYTE,
                                    idxProc,tag,communicator.getComm(),&requests[counterRequest++]);
                            idxToWriteRecvedDatas += globalSendRecvMap[idxProc*nbProcs+myRank];
                        }
                    }
                }

                //Third, send
                for(int idxProc=0 ; idxProc<nbProcs ; ++idxProc){
                    if(toSend[idxProc]){ //If i have something for idxProc
                        MPI_Isend(&leavesArray[leavesIndices[idxToSend[idxProc]]],int(sizeof(ParticleClass))*int(globalSendRecvMap[myRank*nbProcs+idxProc]),MPI_BYTE,
                                idxProc,tag,communicator.getComm(),&requests[counterRequest++]);
                    }
                }
                //Wait for the comm :
                MPI_Waitall(counterRequest,requests,MPI_STATUSES_IGNORE);
                delete[] requests;
            }

            for(FSize idPartsToStore=0; idPartsToStore<finalTotParts ; ++idPartsToStore){
                particlesSaver->push(finalPartBuffer[idPartsToStore]);
            }

            delete[] finalPartBuffer;
            delete[] globalSendRecvMap;
            delete[] idxToSend;
            delete[] toSend;
            delete[] numberOfLeavesPerProc;
        }
    delete[] leavesArray;
    delete[] leavesIndices;
  }

 
public:
  
  /**
   * Those three function get through pubilc/private member issues for testing purpose
   */
  static void testSortParticlesFromArray( const FMpi::FComm& communicator, const ParticleClass array[], const FSize size, const SortingType type,
					  const FPoint& centerOfBox, const FReal boxWidth,
					  const int TreeHeight, IndexedParticle**const outputArray, FSize* const outputSize){
    SortParticlesFromArray(communicator, array, size, type, centerOfBox, boxWidth, TreeHeight,&outputArray, &outputSize);
  }

  
  static void testMergeLeaves(const FMpi::FComm& communicator, IndexedParticle*& workingArray, FSize* workingSize,
			      FSize ** leavesIndices, ParticleClass** leavesArray, FSize* const leavesSize){
    MergeLeaves(communicator,workingArray,workingSize,leavesIndices,leavesArray,leavesSize);
  }

  template<class ContainerClass>
  static void testEqualizeAndFillTree(const FMpi::FComm& communicator,  ContainerClass* particlesSaver,
				      FSize * leavesIndices, ParticleClass* leavesArray, FSize& nbLeavesInIntervals, FSize nbParts, 
				      FAbstractBalanceAlgorithm * balancer){
    EqualizeAndFillTree(communicator,particlesSaver,leavesIndices, leavesArray, nbLeavesInIntervals, nbParts, balancer);
  }

  //////////////////////////////////////////////////////////////////////////
  // The builder function
  //////////////////////////////////////////////////////////////////////////
  
  template <class ContainerClass>
  static void ArrayToTree(const FMpi::FComm& communicator, const ParticleClass array[], const FSize size,
			  const FPoint& boxCenter, const FReal boxWidth, const int treeHeight,
			  ContainerClass* particleSaver, FAbstractBalanceAlgorithm* balancer,const SortingType type = QuickSort){
      
    IndexedParticle* particlesArray = nullptr;
    FSize particlesSize = 0;
    SortParticlesFromArray(communicator, array, size, type, boxCenter, boxWidth, treeHeight,
			   &particlesArray, &particlesSize);

    ParticleClass* leavesArray = nullptr;
    FSize leavesSize = 0;
    FSize * leavesIndices = nullptr;
    
    MergeLeaves(communicator, particlesArray, &particlesSize, &leavesIndices, &leavesArray, &leavesSize);

    EqualizeAndFillTree(communicator, particleSaver, leavesIndices, leavesArray, leavesSize, particlesSize, balancer);

    /** To produce stats after the Equalize phase
     */
    int* nbPartsPerProc;
    int myParts = particleSaver->getSize();
    int nbProc = communicator.processCount();
	
    if(communicator.processId()==0){
      nbPartsPerProc = new int[nbProc];
    }
	
    MPI_Gather(&myParts,1,MPI_INT,nbPartsPerProc,1,MPI_INT,0,communicator.getComm());

    if(communicator.processId()==0){
      float av=0;
      int min=myParts, max=myParts;
      for(int id=0; id<nbProc ; id++){
	if(nbPartsPerProc[id] > max){
	  max=nbPartsPerProc[id];
	}
	if(nbPartsPerProc[id] < min){
	  min=nbPartsPerProc[id];
	}
	av += float(nbPartsPerProc[id]);
      }
      av /= float(nbProc);
      printf("End of Equalize Phase : \n \t Min number of parts : %d \n \t Max number of parts : %d \n \t Average number of parts : %e \n",
	     min,max,av);
      delete[] nbPartsPerProc;
    }
  }


};

#endif // FMPITREEBUILDER_H
