
// Keep in private GIT
#ifndef FGROUPTREE_HPP
#define FGROUPTREE_HPP
#include <vector>
#include <functional>

#include "../../Utils/FAssert.hpp"
#include "../../Utils/FPoint.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Containers/FCoordinateComputer.hpp"
#include "FGroupOfCells.hpp"
#include "FGroupOfParticles.hpp"
#include "FGroupAttachedLeaf.hpp"



template <class FReal, class CompositeCellClass, class SymboleCellClass, class PoleCellClass, class LocalCellClass,
          class GroupAttachedLeafClass, unsigned NbSymbAttributes, unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FGroupTreeBerenger {
public:
    typedef GroupAttachedLeafClass BasicAttachedClass;
    typedef FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle,AttributeClass> ParticleGroupClass;
    typedef FGroupOfCells<CompositeCellClass, SymboleCellClass, PoleCellClass, LocalCellClass> CellGroupClass;

protected:
    //< height of the tree (1 => only the root)
    const int treeHeight;
    //< max number of cells in a block
    const int nbElementsPerBlock;
    //< all the blocks of the tree
    std::vector<CellGroupClass*>* cellBlocksPerLevel;
    //< all the blocks of leaves
    std::vector<ParticleGroupClass*> particleBlocks;

    //< the space system center
    const FPoint<FReal> boxCenter;
    //< the space system corner (used to compute morton index)
    const FPoint<FReal> boxCorner;
    //< the space system width
    const FReal boxWidth;
    //< the width of a box at width level
    const FReal boxWidthAtLeafLevel;

public:
    typedef typename std::vector<CellGroupClass*>::iterator CellGroupIterator;
    typedef typename std::vector<CellGroupClass*>::const_iterator CellGroupConstIterator;
    typedef typename std::vector<ParticleGroupClass*>::iterator ParticleGroupIterator;
    typedef typename std::vector<ParticleGroupClass*>::const_iterator ParticleGroupConstIterator;

    /**
     * This constructor create a group tree from a particle container index.
     * The morton index are computed and the particles are sorted in a first stage.
     * Then the leaf level is done.
     * Finally the other leve are proceed one after the other.
     * It should be easy to make it parallel using for and tasks.
     * If no limite give inLeftLimite = -1
     */
    template<class ParticleContainer>
    FGroupTreeBerenger(const int inTreeHeight, const FReal inBoxWidth, const FPoint<FReal>& inBoxCenter,
               const int inNbElementsPerBlock, ParticleContainer* inParticlesContainer,
               const bool particlesAreSorted, std::vector<MortonIndex> const& distributedMortonIndex, MortonIndex inLeftLimite = -1):
            treeHeight(inTreeHeight),nbElementsPerBlock(inNbElementsPerBlock),cellBlocksPerLevel(nullptr),
            boxCenter(inBoxCenter), boxCorner(inBoxCenter,-(inBoxWidth/2)), boxWidth(inBoxWidth),
            boxWidthAtLeafLevel(inBoxWidth/FReal(1<<(inTreeHeight-1))){

        cellBlocksPerLevel = new std::vector<CellGroupClass*>[treeHeight];

        MortonIndex* currentBlockIndexes = new MortonIndex[nbElementsPerBlock];
        // First we work at leaf level
        {
            // Build morton index for particles
            struct ParticleSortingStruct{
                FSize originalIndex;
                MortonIndex mindex;
            };
            // Convert position to morton index
            const FSize nbParticles = inParticlesContainer->getNbParticles();
            ParticleSortingStruct* particlesToSort = new ParticleSortingStruct[nbParticles];
            {
                const FReal* xpos = inParticlesContainer->getPositions()[0];
                const FReal* ypos = inParticlesContainer->getPositions()[1];
                const FReal* zpos = inParticlesContainer->getPositions()[2];

                for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
                    const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPositionAndCorner<FReal>(this->boxCorner, this->boxWidth,
                                                                                                       treeHeight,
                                                                                                       FPoint<FReal>(xpos[idxPart], ypos[idxPart], zpos[idxPart]) );
                    const MortonIndex particleIndex = host.getMortonIndex(treeHeight-1);
                    particlesToSort[idxPart].mindex = particleIndex;
                    particlesToSort[idxPart].originalIndex = idxPart;
                }
            }

            // Sort if needed
            if(particlesAreSorted == false){
                FQuickSort<ParticleSortingStruct, FSize>::QsOmp(particlesToSort, nbParticles, [](const ParticleSortingStruct& v1, const ParticleSortingStruct& v2){
                    return v1.mindex <= v2.mindex;
                });
            }

            FAssertLF(nbParticles == 0 || inLeftLimite < particlesToSort[0].mindex);

            // Convert to block
            const int idxLevel = (treeHeight - 1);
            FSize* nbParticlesPerLeaf = new FSize[nbElementsPerBlock];
            FSize firstParticle = 0;
            // We need to proceed each group in sub level
            while(firstParticle != nbParticles){
                int sizeOfBlock = 0;
                FSize lastParticle = firstParticle;
                // Count until end of sub group is reached or we have enough cells (or until it reach the next mortonIndex boundary) TODO
                while(sizeOfBlock < nbElementsPerBlock && lastParticle < nbParticles){
                    if(sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != particlesToSort[lastParticle].mindex){
                        currentBlockIndexes[sizeOfBlock] = particlesToSort[lastParticle].mindex;
                        nbParticlesPerLeaf[sizeOfBlock]  = 1;
                        sizeOfBlock += 1;
                    }
                    else{
                        nbParticlesPerLeaf[sizeOfBlock-1] += 1;
                    }
                    lastParticle += 1;
                }
                while(lastParticle < nbParticles && currentBlockIndexes[sizeOfBlock-1] == particlesToSort[lastParticle].mindex){
                    nbParticlesPerLeaf[sizeOfBlock-1] += 1;
                    lastParticle += 1;
                }

                // Create a group
                CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                                                                 currentBlockIndexes[sizeOfBlock-1]+1,
                                                                 sizeOfBlock);
                ParticleGroupClass*const newParticleBlock = new ParticleGroupClass(currentBlockIndexes[0],
                        currentBlockIndexes[sizeOfBlock-1]+1,
                        sizeOfBlock, lastParticle-firstParticle);

                // Init cells
                size_t nbParticlesOffsetBeforeLeaf = 0;
                FSize offsetParticles = firstParticle;
                for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                    newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                    CompositeCellClass newNode = newBlock->getCompleteCell(cellIdInBlock);
                    newNode.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                    FTreeCoordinate coord;
                    coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock], idxLevel);
                    newNode.setCoordinate(coord);

                    // Add leaf
                    nbParticlesOffsetBeforeLeaf = newParticleBlock->newLeaf(currentBlockIndexes[cellIdInBlock], cellIdInBlock,
                                              nbParticlesPerLeaf[cellIdInBlock], nbParticlesOffsetBeforeLeaf);

                    BasicAttachedClass attachedLeaf = newParticleBlock->template getLeaf<BasicAttachedClass>(cellIdInBlock);
                    // Copy each particle from the original position
                    for(FSize idxPart = 0 ; idxPart < nbParticlesPerLeaf[cellIdInBlock] ; ++idxPart){
                        attachedLeaf.setParticle(idxPart, particlesToSort[idxPart + offsetParticles].originalIndex, inParticlesContainer);
                    }
                    offsetParticles += nbParticlesPerLeaf[cellIdInBlock];
                }

                // Keep the block
                cellBlocksPerLevel[idxLevel].push_back(newBlock);
                particleBlocks.push_back(newParticleBlock);

                sizeOfBlock = 0;
                firstParticle = lastParticle;
            }
            delete[] nbParticlesPerLeaf;
            delete[] particlesToSort;
        }


        // For each level from heigth - 2 to 1
        for(int idxLevel = treeHeight-2; idxLevel > 0 ; --idxLevel){
            inLeftLimite = (inLeftLimite == -1 ? inLeftLimite : (inLeftLimite>>3));

            CellGroupConstIterator iterChildCells = cellBlocksPerLevel[idxLevel+1].begin();
            const CellGroupConstIterator iterChildEndCells = cellBlocksPerLevel[idxLevel+1].end();

            // Skip blocks that do not respect limit
            while(iterChildCells != iterChildEndCells
                  && ((*iterChildCells)->getEndingIndex()>>3) <= inLeftLimite){
                ++iterChildCells;
            }
            // If lower level is empty or all blocks skiped stop here
            if(iterChildCells == iterChildEndCells){
                break;
            }

            MortonIndex currentCellIndex = (*iterChildCells)->getStartingIndex();
            if((currentCellIndex>>3) <= inLeftLimite) currentCellIndex = ((inLeftLimite+1)<<3);
            int sizeOfBlock = 0;

            // We need to proceed each group in sub level
            while(iterChildCells != iterChildEndCells){
                // Count until end of sub group is reached or we have enough cells
                while(sizeOfBlock < nbElementsPerBlock && iterChildCells != iterChildEndCells ){
                    if((sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != (currentCellIndex>>3))
                            && (*iterChildCells)->exists(currentCellIndex)){
                        currentBlockIndexes[sizeOfBlock] = (currentCellIndex>>3);
                        sizeOfBlock += 1;
                        currentCellIndex = (((currentCellIndex>>3)+1)<<3);
                    }
                    else{
                        currentCellIndex += 1;
                    }
                    // If we are at the end of the sub group, move to next
                    while(iterChildCells != iterChildEndCells && (*iterChildCells)->getEndingIndex() <= currentCellIndex){
                        ++iterChildCells;
                        // Update morton index
                        if(iterChildCells != iterChildEndCells && currentCellIndex < (*iterChildCells)->getStartingIndex()){
                            currentCellIndex = (*iterChildCells)->getStartingIndex();
                        }
                    }
                }

                // If group is full
                if(sizeOfBlock == nbElementsPerBlock || (sizeOfBlock && iterChildCells == iterChildEndCells)){
                    // Create a group
                    CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                                                                     currentBlockIndexes[sizeOfBlock-1]+1,
                                                                     sizeOfBlock);
                    // Init cells
                    for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                        newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                        CompositeCellClass newNode = newBlock->getCompleteCell(cellIdInBlock);
                        newNode.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                        FTreeCoordinate coord;
                        coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock], idxLevel);
                        newNode.setCoordinate(coord);
                    }

                    // Keep the block
                    cellBlocksPerLevel[idxLevel].push_back(newBlock);

                    sizeOfBlock = 0;
                }
            }
        }
        delete[] currentBlockIndexes;
    }

    /** This function dealloc the tree by deleting each block */
    ~FGroupTreeBerenger(){
        for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
            std::vector<CellGroupClass*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            for (CellGroupClass* block: levelBlocks){
                delete block;
            }
        }
        delete[] cellBlocksPerLevel;

        for (ParticleGroupClass* block: particleBlocks){
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
        for (ParticleGroupClass* block: particleBlocks){
            block->forEachLeaf(function);
        }
    }

    /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
    void forEachCell(std::function<void(CompositeCellClass)> function){
        for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
            std::vector<CellGroupClass*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            for (CellGroupClass* block: levelBlocks){
                block->forEachCell(function);
            }
        }
    }

    /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
    void forEachCellWithLevel(std::function<void(CompositeCellClass,const int)> function){
        for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
            std::vector<CellGroupClass*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            for (CellGroupClass* block: levelBlocks){
                block->forEachCell(function, idxLevel);
            }
        }
    }

    /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
    template<class ParticlesAttachedClass>
    void forEachCellLeaf(std::function<void(CompositeCellClass,ParticlesAttachedClass*)> function){
        CellGroupIterator iterCells = cellBlocksPerLevel[treeHeight-1].begin();
        const CellGroupIterator iterEndCells = cellBlocksPerLevel[treeHeight-1].end();

        ParticleGroupIterator iterLeaves = particleBlocks.begin();
        const ParticleGroupIterator iterEndLeaves = particleBlocks.end();

        while(iterCells != iterEndCells && iterLeaves != iterEndLeaves){
            (*iterCells)->forEachCell([&](CompositeCellClass aCell){
                const int leafIdx = (*iterLeaves)->getLeafIndex(aCell.getMortonIndex());
                FAssertLF(leafIdx != -1);
                ParticlesAttachedClass aLeaf = (*iterLeaves)->template getLeaf <ParticlesAttachedClass>(leafIdx);
                FAssertLF(aLeaf.isAttachedToSomething());
                function(aCell, &aLeaf);
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
            std::vector<CellGroupClass*>& levelBlocks = cellBlocksPerLevel[idxLevel];
            std::cout << "Level " << idxLevel << ", there are " << levelBlocks.size() << " groups.\n";
            int idxGroup = 0;
            for (const CellGroupClass* block: levelBlocks){
                std::cout << "\t Group " << (idxGroup++);
                std::cout << "\t Size = " << block->getNumberOfCellsInBlock();
                std::cout << "\t Starting Index = " << block->getStartingIndex();
                std::cout << "\t Ending Index = " << block->getEndingIndex();
                std::cout << "\t Ratio of usage = " << float(block->getNumberOfCellsInBlock())/float(block->getEndingIndex()-block->getStartingIndex()) << "\n";
            }
        }

        std::cout << "There are " << particleBlocks.size() << " leaf-groups.\n";
        int idxGroup = 0;
        FSize totalNbParticles = 0;
        for (const ParticleGroupClass* block: particleBlocks){
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

    /////////////////////////////////////////////////////////
    // Algorithm function
    /////////////////////////////////////////////////////////

    int getHeight() const {
        return treeHeight;
    }

    CellGroupIterator cellsBegin(const int inLevel){
        FAssertLF(inLevel < treeHeight);
        return cellBlocksPerLevel[inLevel].begin();
    }

    CellGroupConstIterator cellsBegin(const int inLevel) const {
        FAssertLF(inLevel < treeHeight);
        return cellBlocksPerLevel[inLevel].begin();
    }

    CellGroupIterator cellsEnd(const int inLevel){
        FAssertLF(inLevel < treeHeight);
        return cellBlocksPerLevel[inLevel].end();
    }

    CellGroupConstIterator cellsEnd(const int inLevel) const {
        FAssertLF(inLevel < treeHeight);
        return cellBlocksPerLevel[inLevel].end();
    }

    int getNbCellGroupAtLevel(const int inLevel) const {
        FAssertLF(inLevel < treeHeight);
        return int(cellBlocksPerLevel[inLevel].size());
    }

    CellGroupClass* getCellGroup(const int inLevel, const int inIdx){
        FAssertLF(inLevel < treeHeight);
        FAssertLF(inIdx < int(cellBlocksPerLevel[inLevel].size()));
        return cellBlocksPerLevel[inLevel][inIdx];
    }

    const CellGroupClass* getCellGroup(const int inLevel, const int inIdx) const {
        FAssertLF(inLevel < treeHeight);
        FAssertLF(inIdx < int(cellBlocksPerLevel[inLevel].size()));
        return cellBlocksPerLevel[inLevel][inIdx];
    }

    ParticleGroupIterator leavesBegin(){
        return particleBlocks.begin();
    }

    ParticleGroupConstIterator leavesBegin() const {
        return particleBlocks.begin();
    }

    ParticleGroupIterator leavesEnd(){
        return particleBlocks.end();
    }

    ParticleGroupConstIterator leavesEnd() const {
        return particleBlocks.end();
    }

    int getNbParticleGroup() const {
        return int(particleBlocks.size());
    }

    ParticleGroupClass* getParticleGroup(const int inIdx){
        FAssertLF(inIdx < int(particleBlocks.size()));
        return particleBlocks[inIdx];
    }

    const ParticleGroupClass* getParticleGroup(const int inIdx) const {
        FAssertLF(inIdx < int(particleBlocks.size()));
        return particleBlocks[inIdx];
    }
};

#endif // FGROUPTREE_HPP
