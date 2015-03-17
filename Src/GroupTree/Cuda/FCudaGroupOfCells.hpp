// @SCALFMM_PRIVATE
#ifndef FCUDAGROUPOFCELLS_HPP
#define FCUDAGROUPOFCELLS_HPP

#include "FCudaGlobal.hpp"

#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

/**
* @brief The FCudaGroupOfCells class manages the cells in block allocation.
*/
template <class CellClass>
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
    unsigned char*      blockCells;
    //< This value is for not used cells
    static const MortonIndex CellIsEmptyFlag = -1;

public:
    __device__ FCudaGroupOfCells()
        : allocatedMemoryInByte(0), memoryBuffer(nullptr),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr){
    }

    __device__ void reset(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte){
        // Move the pointers to the correct position
        allocatedMemoryInByte = (inAllocatedMemoryInByte);
        memoryBuffer = (inBuffer);
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        blockCells          = reinterpret_cast<unsigned char*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));
    }

    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    __device__ FCudaGroupOfCells(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        blockCells          = reinterpret_cast<unsigned char*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));
    }

    /** Call the destructor of cells and dealloc block memory */
    __device__ ~FCudaGroupOfCells(){
    }

    /** Give access to the buffer to send the data */
    __device__ const unsigned char* getRawBuffer() const{
        return memoryBuffer;
    }

    /** The the size of the allocated buffer */
    __device__ int getBufferSizeInByte() const {
        return allocatedMemoryInByte;
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
    __device__ CellClass* getCell(const MortonIndex inIndex){
        if( exists(inIndex) ) return (CellClass*)(&blockCells[sizeof(CellClass)*blockIndexesTable[inIndex-blockHeader->startingIndex]]);
        else return nullptr;
    }

    /** Return the address of the cell if it exists (or NULL) */
    __device__ const CellClass* getCell(const MortonIndex inIndex) const {
        if( exists(inIndex) ) return (CellClass*)(&blockCells[sizeof(CellClass)*blockIndexesTable[inIndex-blockHeader->startingIndex]]);
        else return nullptr;
    }
};

#endif // FCUDAGROUPOFCELLS_HPP

