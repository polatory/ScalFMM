
#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Containers/FTreeCoordinate.hpp"

#include <list>
#include <functional>


//This is temporary. The idea is to avoid the situation where two
//cells are far away in the tree, but the constructor put them inside
//the same block, generating a huge header.
//Need to be greater than nbCellPerBlocks, obviously.
#define SIZEMAXFORHEAD 500

template <class CellClass>
class FBlockedTree {
  //< This value is for not used cells
  static const int CellIsEmptyFlag = -1;

protected:
  /**
   * @brief The FBlockOfCells class manages the cells in block allocation.
   */
  class FBlockOfCells {
    /** One header is allocated at the beginning of each block */
    struct BlockHeader{
      MortonIndex startingIndex;
      MortonIndex endingIndex;
      int numberOfCellsInBlock;
      int blockIndexesTableSize;
    };
        
  protected:
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    int*            blockIndexesTable;
    //< Pointer to the cells inside the block memory
    CellClass*      blockCells;
    //<Max Number of Cell per Block defined in the tree
    int nbCellsPerBlocks;
    //< This value is for not used cells
    static const MortonIndex CellIsEmptyFlag = -1;
    
  public:
    /**
     * @brief FBlockOfCells
     * @param inStartingIndex first cell morton index
     * @param inEndingIndex last cell morton index + 1
     * @param inNumberOfCells total number of cells in the interval (should be <= inEndingIndex-inEndingIndex)
     */
    FBlockOfCells(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfCells, const FBlockedTree *const tree )
      : memoryBuffer(0), blockHeader(0), blockIndexesTable(0), blockCells(0), nbCellsPerBlocks(tree->nbCellsPerBlocks){
      // Find the number of cell to allocate in the blocks
      const int blockIndexesTableSize = int(inEndingIndex-inStartingIndex);
      FAssertLF(inNumberOfCells <= blockIndexesTableSize);
      // Total number of bytes in the block
      const size_t memoryToAlloc = sizeof(BlockHeader) + (blockIndexesTableSize*sizeof(int)) + (inNumberOfCells*sizeof(CellClass));

      // Allocate
      memoryBuffer = new unsigned char[memoryToAlloc];
      FAssertLF(memoryBuffer);
      memset(memoryBuffer, 0, memoryToAlloc);

      // Move the pointers to the correct position
      blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
      blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
      blockCells          = reinterpret_cast<CellClass*>(memoryBuffer+sizeof(BlockHeader)+(blockIndexesTableSize*sizeof(int)));

      // Init header
      blockHeader->startingIndex = inStartingIndex;
      blockHeader->endingIndex   = inEndingIndex;
      blockHeader->numberOfCellsInBlock  = inNumberOfCells;
      blockHeader->blockIndexesTableSize = blockIndexesTableSize;

      // Set all index to not used
      for(int idxCellPtr = 0 ; idxCellPtr < blockIndexesTableSize ; ++idxCellPtr){
	blockIndexesTable[idxCellPtr] = CellIsEmptyFlag;
      }
    }
      
    /**@ init function allocate memory for a new Block
     */
    void init(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfCells){
      // Find the number of cell to allocate in the blocks
      const int blockIndexesTableSize = int(inEndingIndex-inStartingIndex);
      FAssertLF(inNumberOfCells <= blockIndexesTableSize);
      // Total number of bytes in the block
      const size_t memoryToAlloc = sizeof(BlockHeader) + (blockIndexesTableSize*sizeof(int)) + (inNumberOfCells*sizeof(CellClass));
	
      // Allocate
      memoryBuffer = new unsigned char[memoryToAlloc];
      FAssertLF(memoryBuffer);
      memset(memoryBuffer, 0, memoryToAlloc);
	
      // Move the pointers to the correct position
      blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
      blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
      blockCells          = reinterpret_cast<CellClass*>(memoryBuffer+sizeof(BlockHeader)+(blockIndexesTableSize*sizeof(int)));
	
      // Init header
      blockHeader->startingIndex = inStartingIndex;
      blockHeader->endingIndex   = inEndingIndex;
      blockHeader->numberOfCellsInBlock  = inNumberOfCells;
      blockHeader->blockIndexesTableSize = blockIndexesTableSize;
	
      // Set all index to not used
      for(int idxCellPtr = 0 ; idxCellPtr < blockIndexesTableSize ; ++idxCellPtr){
	blockIndexesTable[idxCellPtr] = CellIsEmptyFlag;
      }
    }
      
    /** @brief This constructor take an array of Morton index to
     * create the block and the cells !!
     */
    FBlockOfCells(MortonIndex head[],const FBlockedTree *const tree)
      : memoryBuffer(0), blockHeader(0), blockIndexesTable(0), blockCells(0),nbCellsPerBlocks(tree->nbCellsPerBlocks){
      //Store the start and end
      MortonIndex start = head[0];
      MortonIndex end = start;
      int count = 0;
      // Find the number of cell to allocate in the blocks
      for(int idxHead = 0 ; idxHead<nbCellsPerBlocks ; idxHead++){
	if (head[idxHead] != CellIsEmptyFlag){
	  count++;
	  end = head[idxHead];
	}
	else{
	  break;
	}
      }
      //allocation of memory
      init(start,end+1,count);
      //allocation of cells
      for(int idx=0 ; idx<count ; idx++){
	this->newCell(head[idx], idx);
	(this->getCell(head[idx]))->setMortonIndex(head[idx]);
	//(this->getCell(head[idx]))->setCoordinate();
	//printf("Done");
      }
    }
      
      
    /** Call the destructor of cells and dealloc block memory */
    ~FBlockOfCells(){
      for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
	if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
	  (&blockCells[blockIndexesTable[idxCellPtr]])->~CellClass();
	}
      }
      delete[] memoryBuffer;
    }

    /** The index of the fist cell (set from the constructor) */
    MortonIndex getStartingIndex() const {
      return blockHeader->startingIndex;
    }

    /** The index of the last cell + 1 (set from the constructor) */
    MortonIndex getEndingIndex() const {
      return blockHeader->endingIndex;
    }

    /** The number of cell (set from the constructor) */
    int getNumberOfCellsInBlock() const {
      return blockHeader->numberOfCellsInBlock;
    }

    /** The size of the interval endingIndex-startingIndex (set from the constructor) */
    int getSizeOfInterval() const {
      return blockHeader->blockIndexesTableSize;
    }

    /** Return true if inIndex should be located in the current block */
    bool isInside(const MortonIndex inIndex) const{
      return blockHeader->startingIndex <= inIndex && inIndex < blockHeader->endingIndex;
    }

    /** Return true if inIndex is located in the current block and is not empty */
    bool exists(const MortonIndex inIndex) const {
      return isInside(inIndex) && (blockIndexesTable[inIndex-blockHeader->startingIndex] != CellIsEmptyFlag);
    }

    /** Return the address of the cell if it exists (or NULL) */
    CellClass* getCell(const MortonIndex inIndex){
      if( exists(inIndex) ) return &blockCells[blockIndexesTable[inIndex-blockHeader->startingIndex]];
      else return 0;
    }

    /** Return the address of the cell if it exists (or NULL) */
    const CellClass* getCell(const MortonIndex inIndex) const {
      if( exists(inIndex) ) return &blockCells[blockIndexesTable[inIndex-blockHeader->startingIndex]];
      else return 0;
    }

    /** Allocate a new cell by calling its constructor */
    template<typename... CellConstructorParams>
    void newCell(const MortonIndex inIndex, const int id, CellConstructorParams... args){
      FAssertLF(isInside(inIndex));
      FAssertLF(!exists(inIndex));
      FAssertLF(id < blockHeader->blockIndexesTableSize);
      new((void*)&blockCells[id]) CellClass(args...);
      blockIndexesTable[inIndex-blockHeader->startingIndex] = id;
    }
      
    /** Iterate on each allocated cells */
    template<typename... FunctionParams>
    void forEachCell(std::function<void(CellClass*)> function, FunctionParams... args){
      for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
	if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
	  function(&blockCells[blockIndexesTable[idxCellPtr]], args...);
	}
      }
    }
  };

  //< height of the tree (1 => only the root)
  const int treeHeight;
  //< max number of cells in a block
  const int nbCellsPerBlocks;
  //< all the blocks of the tree
  std::list<FBlockOfCells*>* cellBlocksPerLevel;

public:
  /** This constructor create a blocked octree from a usual octree
   * The cell are allocated as in the usual octree (no copy constructor are called!)
   * Once allocated each cell receive its morton index and tree coordinate.
   * No blocks are allocated at level 0.
   */
  template<class OctreeClass>
  FBlockedTree(const int inTreeHeight, const int inNbCellsPerBlock, OctreeClass*const inOctreeSrc)
    : treeHeight(inTreeHeight), nbCellsPerBlocks(inNbCellsPerBlock), cellBlocksPerLevel(0){
    cellBlocksPerLevel = new std::list<FBlockOfCells*>[treeHeight];

    // Iterate on the tree and build
    typename OctreeClass::Iterator octreeIterator(inOctreeSrc);
    octreeIterator.gotoBottomLeft();

    // For each level from heigth - 1 to 1
    for(int idxLevel = treeHeight-1; idxLevel > 0 ; --idxLevel){
      typename OctreeClass::Iterator avoidGotoLeft = octreeIterator;
      // For each cell at this level
      do {
	typename OctreeClass::Iterator blockIteratorInOctree = octreeIterator;
	// Move the iterator per nbCellsPerBlocks (or until it cannot move right)
	int sizeOfBlock = 1;
	while(sizeOfBlock < nbCellsPerBlocks && octreeIterator.moveRight()){
	  sizeOfBlock += 1;
	}

	// Create a block with the apropriate parameters
	FBlockOfCells*const newBlock = new FBlockOfCells(blockIteratorInOctree.getCurrentGlobalIndex(),
							 octreeIterator.getCurrentGlobalIndex()+1,
							 sizeOfBlock,this);
	// Initialize each cell of the block
	int cellIdInBlock = 0;
	while(cellIdInBlock != sizeOfBlock){
	  const CellClass*const oldNode = blockIteratorInOctree.getCurrentCell();
	  newBlock->newCell(oldNode->getMortonIndex(), cellIdInBlock);

	  CellClass* newNode = newBlock->getCell(oldNode->getMortonIndex());
	  newNode->setMortonIndex(oldNode->getMortonIndex());
	  newNode->setCoordinate(oldNode->getCoordinate());

	  cellIdInBlock += 1;
	  blockIteratorInOctree.moveRight();
	}

	// Keep the block
	cellBlocksPerLevel[idxLevel].push_back(newBlock);

	// If we can move right then add another block
      } while(octreeIterator.moveRight());

      avoidGotoLeft.moveUp();
      octreeIterator = avoidGotoLeft;
    }
  }
  
  /** This constructor build the BlockOctree from an Octree, but only
   * the cells at leaf level are read. The cells ares constructed with
   * several walk through the leafs.
   */
  template<class OctreeClass>
  FBlockedTree(const int inTreeHeight, const int inNbCellsPerBlock, OctreeClass*const inOctreeSrc, int FLAG):
    treeHeight(inTreeHeight),nbCellsPerBlocks(inNbCellsPerBlock),cellBlocksPerLevel(0)
  {
    cellBlocksPerLevel = new std::list<FBlockOfCells*>[treeHeight];
    int *nbCellPerLevel = new int[treeHeight];
    inOctreeSrc->getNbCellsPerLevel(nbCellPerLevel);
    int nbLeaf = nbCellPerLevel[treeHeight-1];
    delete[] nbCellPerLevel;
     
    //We get the Morton index of all the cells at leaf level and store
    //them in order to build our blocks
    MortonIndex *leafsIdx = new MortonIndex[nbLeaf];
    int idxLeafs = 0;
    //Iterator on source octree
    typename OctreeClass::Iterator octreeIterator(inOctreeSrc);
    octreeIterator.gotoBottomLeft();
    do{
      leafsIdx[idxLeafs] = octreeIterator.getCurrentGlobalIndex();
      idxLeafs++;
    }
    while(octreeIterator.moveRight());
    
    //Now the work start !!1!
    idxLeafs = 0;
    
    //Temporary Header, will be filled
    MortonIndex * head = new MortonIndex[inNbCellsPerBlock];
    memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*inNbCellsPerBlock);
    //Creation and addition of all the cells at leaf level
    while(idxLeafs < nbLeaf){
      
      //Initialisation of header
      memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*inNbCellsPerBlock);
      memcpy(head,&leafsIdx[idxLeafs],FMath::Min(sizeof(MortonIndex)*inNbCellsPerBlock,sizeof(MortonIndex)*(nbLeaf-idxLeafs)));
      idxLeafs += inNbCellsPerBlock;
    
      //creation of the block and addition to the list
      FBlockOfCells * tempBlock = new FBlockOfCells(head,this);
      cellBlocksPerLevel[treeHeight-1].push_back(tempBlock);
    }
    printf("Cells at Leaf level created. Nb of Blocks used : %lu\n",cellBlocksPerLevel[treeHeight-1].size());
    
    typename std::list<FBlockOfCells*>::const_iterator curBlockPrint;
    int i =0;
    for(curBlockPrint = cellBlocksPerLevel[treeHeight-1].cbegin() ; curBlockPrint != cellBlocksPerLevel[treeHeight-1].end() ; ++curBlockPrint){
      printf("Block : %d, from %lld to %lld, %d cells \n",i,(*curBlockPrint)->getStartingIndex(),
	     (*curBlockPrint)->getEndingIndex(),(*curBlockPrint)->getNumberOfCellsInBlock());
      i++;
    }
    
    //Creation of the cells at others level
    //Loop over levels
    for(int idxLevel = treeHeight-2 ; idxLevel > 0 ; --idxLevel){
      
      //Set the storing system (a simple array) to CellIsEmptyFlag
      memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*inNbCellsPerBlock);
      int idxHead = -1;
      //Iterator over the list at a deeper level (READ)
      typename std::list<FBlockOfCells*>::const_iterator curBlockRead;
      for(curBlockRead = cellBlocksPerLevel[idxLevel+1].begin() ; curBlockRead != cellBlocksPerLevel[idxLevel+1].end() ; ++curBlockRead){
	//Loop over cells in READ list
	for(MortonIndex idxCell = (*curBlockRead)->getStartingIndex() ; idxCell < (*curBlockRead)->getEndingIndex() ; ++idxCell){

	  MortonIndex newIdx(idxCell >> 3);
	  
	  //For Head Increment, we need to know when parent change
	  if((newIdx != ((idxCell-1) >> 3))){
	    //test if (Read cell exists) 
	    if((*curBlockRead)->exists(idxCell)){
	      idxHead++; //cell at deeper level exist, so we increment the count to store it
	      
	      //Test if a new Block is necessary
	      if(idxHead < inNbCellsPerBlock){
		//No need for a Block, just the cell is created
		head[idxHead] = newIdx;
	      }
	      else{
		//Creation of the block from head, then reset head, and
		//storage of new idx in new head
		FBlockOfCells * tempBlock = new FBlockOfCells(head,this);
		cellBlocksPerLevel[idxLevel].push_back(tempBlock);
	      
		//Need a new block
		memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*nbCellsPerBlocks);
		idxHead = 0;
		head[idxHead] = newIdx;
	      }
	    }
	  }
	}
      }
      //Before changing Level, need to close current Block
      FBlockOfCells * tempBlock = new FBlockOfCells(head,this);
      cellBlocksPerLevel[idxLevel].push_back(tempBlock);
      //Printing Information about the blocks
      printf("Cells at level %d created. Nb of Blocks used : %lu\n",idxLevel,cellBlocksPerLevel[idxLevel].size());
     
      for(curBlockPrint = cellBlocksPerLevel[idxLevel].cbegin() ; curBlockPrint != cellBlocksPerLevel[idxLevel].end() ; ++curBlockPrint){
	printf("Block : %d, from %lld to %lld, %d cells \n",i,(*curBlockPrint)->getStartingIndex(),
	       (*curBlockPrint)->getEndingIndex(),(*curBlockPrint)->getNumberOfCellsInBlock());
	i++;
      }
    }
  }
    

  /** This function dealloc the tree by deleting each block */
  ~FBlockedTree(){
    for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
      std::list<FBlockOfCells*>& levelBlocks = cellBlocksPerLevel[idxLevel];
      for (FBlockOfCells* &block: levelBlocks){
	delete block;
      }
    }
    delete[] cellBlocksPerLevel;
  }


  /////////////////////////////////////////////////////////
  // Lambda function to apply to all member
  /////////////////////////////////////////////////////////

  /**
   * @brief forEachLeaf iterate on the leaf and apply the function
   * @param function
   */
  //void forEachLeaf(std::function<void(LeafClass*)> function){
  // TODO
  //}

  /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
  void forEachCell(std::function<void(CellClass*)> function){
    for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
      std::list<FBlockOfCells*>& levelBlocks = cellBlocksPerLevel[idxLevel];
      for (FBlockOfCells* &block: levelBlocks){
	block->forEachCell(function);
      }
    }
  }

  /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
  void forEachCellWithLevel(std::function<void(CellClass*,const int)> function){
    for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
      std::list<FBlockOfCells*>& levelBlocks = cellBlocksPerLevel[idxLevel];
      for (FBlockOfCells* &block: levelBlocks){
	block->forEachCell(function, idxLevel);
      }
    }
  }
  
  /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
  //void forEachCellLeaf(std::function<void(CellClass*,LeafClass*)> function){
  // TODO
  //}
};




#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Files/FFmaBinLoader.hpp"



int main(int argc, char* argv[]){
  static const int P = 9;
  typedef FRotationCell<P>               CellClass;
  typedef FP2PParticleContainer<>          ContainerClass;

  typedef FSimpleLeaf< ContainerClass >                     LeafClass;
  typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
  //typedef FRotationKernel< CellClass, ContainerClass , P>   KernelClass;
  typedef FBlockedTree< CellClass >  BlockedOctreeClass;

  FTic counter;
  const int NbLevels      = FParameters::getValue(argc,argv,"-h", 5);
  const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
  const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma.double");

  FFmaBinLoader loader(filename);
  FAssertLF(loader.isOpen());

  OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

  for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
    FPoint particlePosition;
    FReal physicalValue;
    loader.fillParticle(&particlePosition,&physicalValue);
    tree.insert(particlePosition, physicalValue );
  }

  std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.tacAndElapsed() << "s)." << std::endl;

  const int blockSize      = FParameters::getValue(argc,argv,"-bs", 250);

  counter.tic();
  BlockedOctreeClass blockedTree(NbLevels, blockSize, &tree);
  std::cout << "Done  " << "(@Converting the tree = " << counter.tacAndElapsed() << "s). Block size is " << blockSize << "." << std::endl;

  blockedTree.forEachCell([&](CellClass * cell){
      std::cout << "Cell " << cell->getMortonIndex() << std::endl;
    });
  
  return 0;
}
