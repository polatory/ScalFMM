// @SCALFMM_PRIVATE
#ifndef FCUDAGROUPOFCELLS_HPP
#define FCUDAGROUPOFCELLS_HPP

#include "FCudaGlobal.hpp"
#include "FCudaCompositeCell.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

/**
* @brief The FCudaGroupOfCells class manages the cells in block allocation.
*/
template <class SymboleCellClass, class PoleCellClass, class LocalCellClass>
class FCudaGroupOfCells {
    /** One header is allocated at the beginning of each block */
    struct alignas(FStarPUDefaultAlign::StructAlign) BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int numberOfCellsInBlock;
        int blockIndexesTableSize;
    };

protected:
    //< The size of the memoryBuffer
    int allocatedMemoryInByte;
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    int*            blockIndexesTable;
    //< Pointer to the cells inside the block memory
    SymboleCellClass*      blockCells;
    //< This value is for not used cells
    static const MortonIndex CellIsEmptyFlag = -1;

    //< The multipole data
    PoleCellClass* cellMultipoles;
    //< The local data
    LocalCellClass* cellLocals;

public:
    typedef FCudaCompositeCell<SymboleCellClass, PoleCellClass, LocalCellClass> CompleteCellClass;

    __device__ FCudaGroupOfCells()
        : allocatedMemoryInByte(0), memoryBuffer(nullptr),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr){
    }

    __device__ void reset(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                          unsigned char* inCellMultipoles, unsigned char* inCellLocals){
        // Move the pointers to the correct position
        allocatedMemoryInByte = (inAllocatedMemoryInByte);
        memoryBuffer = (inBuffer);
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        blockIndexesTable   = reinterpret_cast<int*>(inBuffer);
        inBuffer += (blockHeader->blockIndexesTableSize*sizeof(int));
        blockCells          = reinterpret_cast<SymboleCellClass*>(inBuffer);
        inBuffer += (sizeof(SymboleCellClass)*blockHeader->numberOfCellsInBlock);
        //FAssertLF(size_t(inBuffer-memoryBuffer) == allocatedMemoryInByte);

        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
    }

    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    __device__ FCudaGroupOfCells(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                                 unsigned char* inCellMultipoles, unsigned char* inCellLocals)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        blockIndexesTable   = reinterpret_cast<int*>(inBuffer);
        inBuffer += (blockHeader->blockIndexesTableSize*sizeof(int));
        blockCells          = reinterpret_cast<SymboleCellClass*>(inBuffer);
        inBuffer += (sizeof(SymboleCellClass)*blockHeader->numberOfCellsInBlock);
        //FAssertLF(size_t(inBuffer-memoryBuffer) == allocatedMemoryInByte);

        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
    }

    /** The index of the fist cell (set from the constructor) */
    __device__ MortonIndex getStartingIndex() const {
        return blockHeader->startingIndex;
    }

    /** The index of the last cell + 1 (set from the constructor) */
    __device__ MortonIndex getEndingIndex() const {
        return blockHeader->endingIndex;
    }

    /** The number of cell (set from the constructor) */
    __device__ int getNumberOfCellsInBlock() const {
        return blockHeader->numberOfCellsInBlock;
    }

    /** The size of the interval endingIndex-startingIndex (set from the constructor) */
    __device__ int getSizeOfInterval() const {
        return blockHeader->blockIndexesTableSize;
    }

    /** Return true if inIndex should be located in the current block */
    __device__ bool isInside(const MortonIndex inIndex) const{
        return blockHeader->startingIndex <= inIndex && inIndex < blockHeader->endingIndex;
    }

    /** Return true if inIndex is located in the current block and is not empty */
    __device__ bool exists(const MortonIndex inIndex) const {
        return isInside(inIndex) && (blockIndexesTable[inIndex-blockHeader->startingIndex] != CellIsEmptyFlag);
    }

    /** Return the address of the cell if it exists (or NULL) */
    __device__ CompleteCellClass getCompleteCell(const MortonIndex inIndex){
        //FAssertLF(cellMultipoles && cellLocals);
        if( exists(inIndex) ){
            CompleteCellClass cell;
            const int cellPos = blockIndexesTable[inIndex-blockHeader->startingIndex];
            cell.symb = &blockCells[cellPos];
            cell.up   = &cellMultipoles[cellPos];
            cell.down = &cellLocals[cellPos];
            return cell;
        }
        else{
            CompleteCellClass cell;
            cell.symb = nullptr;
            cell.up   = nullptr;
            cell.down = nullptr;
            return cell;
        }
    }

    /** Return the address of the cell if it exists (or NULL) */
    __device__ CompleteCellClass getUpCell(const MortonIndex inIndex){
        //FAssertLF(cellMultipoles);
        if( exists(inIndex) ){
            CompleteCellClass cell;
            const int cellPos = blockIndexesTable[inIndex-blockHeader->startingIndex];
            cell.symb = &blockCells[cellPos];
            cell.up   = &cellMultipoles[cellPos];
            cell.down = nullptr;
            return cell;
        }
        else{
            CompleteCellClass cell;
            cell.symb = nullptr;
            cell.up   = nullptr;
            cell.down = nullptr;
            return cell;
        }
    }

    /** Return the address of the cell if it exists (or NULL) */
    __device__ CompleteCellClass getDownCell(const MortonIndex inIndex){
        //FAssertLF(cellLocals);
        if( exists(inIndex) ){
            CompleteCellClass cell;
            const int cellPos = blockIndexesTable[inIndex-blockHeader->startingIndex];
            cell.symb = &blockCells[cellPos];
            cell.up   = nullptr;
            cell.down = &cellLocals[cellPos];
            return cell;
        }
        else{
            CompleteCellClass cell;
            cell.symb = nullptr;
            cell.up   = nullptr;
            cell.down = nullptr;
            return cell;
        }
    }

};

#endif // FCUDAGROUPOFCELLS_HPP

