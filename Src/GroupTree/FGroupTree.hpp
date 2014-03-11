#ifndef FGROUPTREE_HPP
#define FGROUPTREE_HPP

#include "../Utils/FAssert.hpp"
#include "../Containers/FTreeCoordinate.hpp"

#include "FGroupOfCells.hpp"
#include "FGroupOfParticles.hpp"

#include <list>
#include <functional>

template <class CellClass, unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FGroupTree {
    //< This value is for not used cells
    static const int CellIsEmptyFlag = -1;

protected:
    //< height of the tree (1 => only the root)
    const int treeHeight;
    //< max number of cells in a block
    const int nbElementsPerBlock;
    //< all the blocks of the tree
    std::list<FGroupOfCells<CellClass>*>* cellBlocksPerLevel;
    //< all the blocks of leaves
    std::list<FGroupOfParticles<NbAttributesPerParticle,AttributeClass>*> particleBlocks;


public:
    /** This constructor create a blocked octree from a usual octree
   * The cell are allocated as in the usual octree (no copy constructor are called!)
   * Once allocated each cell receive its morton index and tree coordinate.
   * No blocks are allocated at level 0.
   */
    template<class OctreeClass>
    FGroupTree(const int inTreeHeight, const int inNbElementsPerBlock, OctreeClass*const inOctreeSrc)
        : treeHeight(inTreeHeight), nbElementsPerBlock(inNbElementsPerBlock), cellBlocksPerLevel(0){
        cellBlocksPerLevel = new std::list<FGroupOfCells<CellClass>*>[treeHeight];

        // Iterate on the tree and build
        typename OctreeClass::Iterator octreeIterator(inOctreeSrc);
        octreeIterator.gotoBottomLeft();

        { // First leaf level, we create leaves and cells groups
            const int idxLevel = treeHeight-1;
            typename OctreeClass::Iterator avoidGotoLeft = octreeIterator;
            // For each cell at this level
            do {
                typename OctreeClass::Iterator blockIteratorInOctree = octreeIterator;
                // Move the iterator per nbElementsPerBlock (or until it cannot move right)
                int sizeOfBlock = 1;
                int nbParticlesInGroup = octreeIterator.getCurrentLeaf()->getSrc()->getNbParticles();
                while(sizeOfBlock < nbElementsPerBlock && octreeIterator.moveRight()){
                    sizeOfBlock += 1;
                    nbParticlesInGroup += octreeIterator.getCurrentLeaf()->getSrc()->getNbParticles();
                }

                // Create a block with the apropriate parameters
                FGroupOfCells<CellClass>*const newBlock = new FGroupOfCells<CellClass>(blockIteratorInOctree.getCurrentGlobalIndex(),
                                                                 octreeIterator.getCurrentGlobalIndex()+1,
                                                                 sizeOfBlock);
                FGroupOfParticles<NbAttributesPerParticle, AttributeClass>*const newParticleBlock = new FGroupOfParticles<NbAttributesPerParticle, AttributeClass>(blockIteratorInOctree.getCurrentGlobalIndex(),
                                                                 octreeIterator.getCurrentGlobalIndex()+1,
                                                                 sizeOfBlock, nbParticlesInGroup);

                // Initialize each cell of the block
                int cellIdInBlock = 0;
                int nbParticlesBeforeLeaf = 0;
                while(cellIdInBlock != sizeOfBlock){
                    // Add cell
                    const CellClass*const oldNode = blockIteratorInOctree.getCurrentCell();
                    newBlock->newCell(oldNode->getMortonIndex(), cellIdInBlock);

                    CellClass* newNode = newBlock->getCell(oldNode->getMortonIndex());
                    newNode->setMortonIndex(oldNode->getMortonIndex());
                    newNode->setCoordinate(oldNode->getCoordinate());

                    // Add leaf
                    newParticleBlock->newLeaf(oldNode->getMortonIndex(), cellIdInBlock,
                                              blockIteratorInOctree.getCurrentLeaf()->getSrc()->getNbParticles(),
                                              nbParticlesBeforeLeaf, blockIteratorInOctree.getCurrentLeaf()->getSrc(), 0);

                    nbParticlesBeforeLeaf += blockIteratorInOctree.getCurrentLeaf()->getSrc()->getNbParticles();
                    cellIdInBlock += 1;
                    blockIteratorInOctree.moveRight();
                }

                // Keep the block
                cellBlocksPerLevel[idxLevel].push_back(newBlock);
                particleBlocks.push_back(newParticleBlock);

                // If we can move right then add another block
            } while(octreeIterator.moveRight());

            avoidGotoLeft.moveUp();
            octreeIterator = avoidGotoLeft;
        }

        // For each level from heigth - 2 to 1
        for(int idxLevel = treeHeight-2; idxLevel > 0 ; --idxLevel){
            typename OctreeClass::Iterator avoidGotoLeft = octreeIterator;
            // For each cell at this level
            do {
                typename OctreeClass::Iterator blockIteratorInOctree = octreeIterator;
                // Move the iterator per nbElementsPerBlock (or until it cannot move right)
                int sizeOfBlock = 1;
                while(sizeOfBlock < nbElementsPerBlock && octreeIterator.moveRight()){
                    sizeOfBlock += 1;
                }

                // Create a block with the apropriate parameters
                FGroupOfCells<CellClass>*const newBlock = new FGroupOfCells<CellClass>(blockIteratorInOctree.getCurrentGlobalIndex(),
                                                                 octreeIterator.getCurrentGlobalIndex()+1,
                                                                 sizeOfBlock);
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


    ////////////////////////////////////////////////////////////////////
    // Work in progress part, build the tree from an array of particles
    ////////////////////////////////////////////////////////////////////

    /** @brief This private method take an array of Morton index to
   * create the block and the cells, using the constructor of
   * FGroupOfCells !!
   */
    FGroupOfCells<CellClass> * createBlockFromArray(MortonIndex head[]){
        //Store the start and end
        MortonIndex start = head[0];
        MortonIndex end = start;
        int count = 0;
        // Find the number of cell to allocate in the blocks
        for(int idxHead = 0 ; idxHead<nbElementsPerBlock ; idxHead++){
            if (head[idxHead] != CellIsEmptyFlag){
                count++;
                end = head[idxHead];
            }
            // else{
            // 	break;
            // }
        }
        //allocation of memory
        FGroupOfCells<CellClass> * newBlock = new FGroupOfCells<CellClass>(start,end+1,count);
        //allocation of cells
        for(int idx=0 ; idx<count ; idx++){
            newBlock->newCell(head[idx], idx);
            (newBlock->getCell(head[idx]))->setMortonIndex(head[idx]);
            //(this->getCell(head[idx]))->setCoordinate();
        }
        return newBlock;
    }

    /** This constructor build the BlockOctree from an Octree, but only
   * the cells at leaf level are read. The cells ares constructed with
   * several walk through the leafs.
   */
    template<class OctreeClass>
    FGroupTree(const int inTreeHeight, const int inNbElementsPerBlock, OctreeClass*const inOctreeSrc, int FLAG):
        treeHeight(inTreeHeight),nbElementsPerBlock(inNbElementsPerBlock),cellBlocksPerLevel(0)
    {
        cellBlocksPerLevel = new std::list<FGroupOfCells<CellClass>*>[treeHeight];
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
        MortonIndex * head = new MortonIndex[inNbElementsPerBlock];
        memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*inNbElementsPerBlock);
        //Creation and addition of all the cells at leaf level
        while(idxLeafs < nbLeaf){

            //Initialisation of header
            memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*inNbElementsPerBlock);
            memcpy(head,&leafsIdx[idxLeafs],FMath::Min(sizeof(MortonIndex)*inNbElementsPerBlock,sizeof(MortonIndex)*(nbLeaf-idxLeafs)));
            idxLeafs += inNbElementsPerBlock;

            //creation of the block and addition to the list
            FGroupOfCells<CellClass> * tempBlock = createBlockFromArray(head);
            cellBlocksPerLevel[treeHeight-1].push_back(tempBlock);
        }
        delete[] leafsIdx;

        //Creation of the cells at others level
        //Loop over levels
        for(int idxLevel = treeHeight-2 ; idxLevel >= 0 ; --idxLevel){

            //Set the storing system (a simple array) to CellIsEmptyFlag
            memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*inNbElementsPerBlock);
            int idxHead = -1;
            MortonIndex previous = -1;

            //Iterator over the list at a deeper level (READ)
            typename std::list<FGroupOfCells<CellClass>*>::const_iterator curBlockRead;
            for(curBlockRead = cellBlocksPerLevel[idxLevel+1].begin() ; curBlockRead != cellBlocksPerLevel[idxLevel+1].end() ; ++curBlockRead){
                //Loop over cells in READ list
                for(MortonIndex idxCell = (*curBlockRead)->getStartingIndex() ; idxCell < (*curBlockRead)->getEndingIndex() ; ++idxCell){

                    MortonIndex newIdx(idxCell >> 3);
                    //printf("Lvl : %d Current Idx %lld, parent's one : %lld \n",idxLevel,idxCell,newIdx);
                    //For Head Increment, we need to know if this parent has
                    //already been created
                    if(newIdx != previous){
                        //test if (Read cell exists)
                        if((*curBlockRead)->exists(idxCell)){
                            idxHead++; //cell at deeper level exist, so we increment the count to store it

                            //Test if a new Block is necessary
                            if(idxHead < inNbElementsPerBlock){
                                //No need for a Block, just the cell is created
                                head[idxHead] = newIdx;
                                previous = newIdx;
                            }
                            else{
                                //Creation of the block from head, then reset head, and
                                //storage of new idx in new head
                                FGroupOfCells<CellClass> * tempBlock = createBlockFromArray(head);
                                cellBlocksPerLevel[idxLevel].push_back(tempBlock);

                                //Need a new block
                                memset(head,CellIsEmptyFlag,sizeof(MortonIndex)*nbElementsPerBlock);
                                head[0] = newIdx;
                                idxHead = 0;
                                previous = newIdx;
                            }
                        }
                    }
                }
            }
            //Before changing Level, need to close current Block
            FGroupOfCells<CellClass> * tempBlock = createBlockFromArray(head);
            cellBlocksPerLevel[idxLevel].push_back(tempBlock);
        }
        printf("toto \n");
        delete[] head;
    }

    ////////////////////////////////////////////////////////////////////
    // End work in progress part, build the tree from an array of particles
    ////////////////////////////////////////////////////////////////////


    /** This function dealloc the tree by deleting each block */
    ~FGroupTree(){
        for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
            std::list<FGroupOfCells<CellClass>*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            for (FGroupOfCells<CellClass>* block: levelBlocks){
                delete block;
            }
        }
        delete[] cellBlocksPerLevel;

        for (FGroupOfParticles<NbAttributesPerParticle,AttributeClass>* block: particleBlocks){
            delete block;
        }
    }


    /////////////////////////////////////////////////////////
    // Lambda function to apply to all member
    /////////////////////////////////////////////////////////

    /**
   * @brief forEachLeaf iterate on the leaf and apply the function
   * @param function
   */
    template<class ParticlesAttachedClass>
    void forEachLeaf(std::function<void(ParticlesAttachedClass*)> function){
        for (FGroupOfParticles<NbAttributesPerParticle,AttributeClass>* block: particleBlocks){
            block->forEachLeaf(function);
        }
    }

    /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
    void forEachCell(std::function<void(CellClass*)> function){
        for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
            std::list<FGroupOfCells<CellClass>*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            for (FGroupOfCells<CellClass>* block: levelBlocks){
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
            std::list<FGroupOfCells<CellClass>*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            for (FGroupOfCells<CellClass>* block: levelBlocks){
                block->forEachCell(function, idxLevel);
            }
        }
    }

    /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
    template<class ParticlesAttachedClass>
    void forEachCellLeaf(std::function<void(CellClass*,ParticlesAttachedClass*)> function){
        typename std::list<FGroupOfCells<CellClass>*>::iterator iterCells = cellBlocksPerLevel[treeHeight-1].begin();
        const typename std::list<FGroupOfCells<CellClass>*>::iterator iterEndCells = cellBlocksPerLevel[treeHeight-1].end();

        typename std::list<FGroupOfParticles<NbAttributesPerParticle,AttributeClass>*>::iterator iterLeaves = particleBlocks.begin();
        const typename std::list<FGroupOfParticles<NbAttributesPerParticle,AttributeClass>*>::iterator iterEndLeaves = particleBlocks.end();

        while(iterCells != iterEndCells && iterLeaves != iterEndLeaves){
            (*iterCells)->forEachCell([&](CellClass* aCell){
                ParticlesAttachedClass*const aLeaf = (*iterLeaves)->getLeaf(aCell->getMortonIndex());
                FAssertLF(aLeaf);
                function(aCell, aLeaf);
                delete aLeaf;
            });

            ++iterCells;
            ++iterLeaves;
        }

        FAssertLF(iterCells == iterEndCells && iterLeaves == iterEndLeaves);
    }



    /** @brief, for statistic purpose, display each block with number of
   * cell, size of header, starting index, and ending index
   */
    void printInfoBlocks(){
        std::cout << "Group Tree information:\n";
        std::cout << "\t Group Size = " << nbElementsPerBlock << "\n";
        std::cout << "\t Tree height = " << treeHeight << "\n";
        for(int idxLevel = 1 ; idxLevel < treeHeight ; ++idxLevel){
            std::list<FGroupOfCells<CellClass>*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            std::cout << "Level " << idxLevel << ", there are " << levelBlocks.size() << " groups.\n";
            int idxGroup = 0;
            for (const FGroupOfCells<CellClass>* block: levelBlocks){
                std::cout << "\t Group " << (idxGroup++);
                std::cout << "\t Size = " << block->getNumberOfCellsInBlock();
                std::cout << "\t Starting Index = " << block->getStartingIndex();
                std::cout << "\t Ending Index = " << block->getEndingIndex();
                std::cout << "\t Ratio of usage = " << float(block->getNumberOfCellsInBlock())/float(block->getEndingIndex()-block->getStartingIndex()) << "\n";
            }
        }

        std::cout << "There are " << particleBlocks.size() << " leaf-groups.\n";
        int idxGroup = 0;
        int totalNbParticles = 0;
        for (const FGroupOfParticles<NbAttributesPerParticle,AttributeClass>* block: particleBlocks){
            std::cout << "\t Group " << (idxGroup++);
            std::cout << "\t Size = " << block->getNumberOfLeavesInBlock();
            std::cout << "\t Starting Index = " << block->getStartingIndex();
            std::cout << "\t Ending Index = " << block->getEndingIndex();
            std::cout << "\t Nb Particles = " << block->getNbParticlesInGroup();
            std::cout << "\t Ratio of usage = " << float(block->getNumberOfLeavesInBlock())/float(block->getEndingIndex()-block->getStartingIndex()) << "\n";
            totalNbParticles += block->getNbParticlesInGroup();
        }
        std::cout << "There are " << totalNbParticles << " particles.\n";
    }
};

#endif // FGROUPTREE_HPP
