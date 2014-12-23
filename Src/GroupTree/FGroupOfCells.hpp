
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPOFCELLS_HPP
#define FGROUPOFCELLS_HPP

#include "../Utils/FAssert.hpp"
#include "../Utils/FAlignedMemory.hpp"
#include "../Containers/FTreeCoordinate.hpp"

#include <list>
#include <functional>


/**
* @brief The FGroupOfCells class manages the cells in block allocation.
*/
template <class CellClass>
class FGroupOfCells {
    /** One header is allocated at the beginning of each block */
    struct BlockHeader{
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
    CellClass*      blockCells;
    //< This value is for not used cells
    static const MortonIndex CellIsEmptyFlag = -1;

    //< To kown if the object has to delete the memory
    bool deleteBuffer;

public:
    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    FGroupOfCells(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          deleteBuffer(false){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        blockCells          = reinterpret_cast<CellClass*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));
    }

    /**
 * @brief FGroupOfCells
 * @param inStartingIndex first cell morton index
 * @param inEndingIndex last cell morton index + 1
 * @param inNumberOfCells total number of cells in the interval (should be <= inEndingIndex-inEndingIndex)
 */
    FGroupOfCells(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfCells)
        : allocatedMemoryInByte(0), memoryBuffer(nullptr), blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          deleteBuffer(true){
        // Find the number of cell to allocate in the blocks
        const int blockIndexesTableSize = int(inEndingIndex-inStartingIndex);
        FAssertLF(inNumberOfCells <= blockIndexesTableSize);
        // Total number of bytes in the block
        const size_t memoryToAlloc = sizeof(BlockHeader) + (blockIndexesTableSize*sizeof(int)) + (inNumberOfCells*sizeof(CellClass));

        // Allocate
        FAssertLF(0 <= int(memoryToAlloc) && int(memoryToAlloc) < std::numeric_limits<int>::max());
        allocatedMemoryInByte = memoryToAlloc;
        memoryBuffer = (unsigned char*)FAlignedMemory::Allocate32BAligned(memoryToAlloc);
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

    /** Call the destructor of cells and dealloc block memory */
    ~FGroupOfCells(){
        if(deleteBuffer){
            for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
                if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
                    (&blockCells[blockIndexesTable[idxCellPtr]])->~CellClass();
                }
            }
            FAlignedMemory::Dealloc32BAligned(memoryBuffer);
        }
    }

    /** Give access to the buffer to send the data */
    const unsigned char* getRawBuffer() const{
        return memoryBuffer;
    }

    /** The the size of the allocated buffer */
    int getBufferSizeInByte() const {
        return allocatedMemoryInByte;
    }

    /** To know if the object will delete the memory block */
    bool getDeleteMemory() const{
        return deleteBuffer;
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
        else return nullptr;
    }

    /** Return the address of the cell if it exists (or NULL) */
    const CellClass* getCell(const MortonIndex inIndex) const {
        if( exists(inIndex) ) return &blockCells[blockIndexesTable[inIndex-blockHeader->startingIndex]];
        else return nullptr;
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
    void forEachCell(std::function<void(CellClass*, FunctionParams...)> function, FunctionParams... args){
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
            if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
                function(&blockCells[blockIndexesTable[idxCellPtr]], args...);
            }
        }
    }

    void forEachCell(std::function<void(CellClass*)> function){
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
            if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
                function(&blockCells[blockIndexesTable[idxCellPtr]]);
            }
        }
    }
};


#endif // FGROUPOFCELLS_HPP
