
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPOFCELLS_HPP
#define FGROUPOFCELLS_HPP

#include "../../Utils/FAssert.hpp"
#include "../../Utils/FAlignedMemory.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

#include <list>
#include <functional>

/**
* @brief The FGroupOfCells class manages the cells in block allocation.
*/
template <class CompositeCellClass, class SymboleCellClass, class PoleCellClass, class LocalCellClass>
class FGroupOfCells {
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

    //< To kown if the object has to delete the memory
    bool deleteBuffer;

public:
    typedef CompositeCellClass CompleteCellClass;

    FGroupOfCells()
        : allocatedMemoryInByte(0), memoryBuffer(nullptr),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr), deleteBuffer(false){
    }

    void reset(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
               unsigned char* inCellMultipoles, unsigned char* inCellLocals){
        if(deleteBuffer){
            for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
                if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
                    (&blockCells[blockIndexesTable[idxCellPtr]])->~SymboleCellClass();
                    (&cellMultipoles[blockIndexesTable[idxCellPtr]])->~PoleCellClass();
                    (&cellLocals[blockIndexesTable[idxCellPtr]])->~LocalCellClass();
                }
            }
            FAlignedMemory::DeallocBytes(memoryBuffer);
            FAlignedMemory::DeallocBytes(cellMultipoles);
            FAlignedMemory::DeallocBytes(cellLocals);
        }
        // Move the pointers to the correct position
        allocatedMemoryInByte = (inAllocatedMemoryInByte);
        memoryBuffer = (inBuffer);
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        blockCells          = reinterpret_cast<SymboleCellClass*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));
        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
        deleteBuffer = (false);
    }

    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    FGroupOfCells(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                  unsigned char* inCellMultipoles, unsigned char* inCellLocals)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr), deleteBuffer(false){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        blockCells          = reinterpret_cast<SymboleCellClass*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));
        cellMultipoles = (PoleCellClass*)inCellMultipoles;
        cellLocals     = (LocalCellClass*)inCellLocals;
    }

    /**
 * @brief FGroupOfCells
 * @param inStartingIndex first cell morton index
 * @param inEndingIndex last cell morton index + 1
 * @param inNumberOfCells total number of cells in the interval (should be <= inEndingIndex-inEndingIndex)
 */
    FGroupOfCells(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfCells)
        : allocatedMemoryInByte(0), memoryBuffer(nullptr), blockHeader(nullptr), blockIndexesTable(nullptr), blockCells(nullptr),
          cellMultipoles(nullptr), cellLocals(nullptr), deleteBuffer(true){
        // Find the number of cell to allocate in the blocks
        const int blockIndexesTableSize = int(inEndingIndex-inStartingIndex);
        FAssertLF(inNumberOfCells <= blockIndexesTableSize);
        // Total number of bytes in the block
        const size_t memoryToAlloc = sizeof(BlockHeader) + (blockIndexesTableSize*sizeof(int)) + (inNumberOfCells*sizeof(SymboleCellClass));

        // Allocate
        FAssertLF(0 <= int(memoryToAlloc) && int(memoryToAlloc) < std::numeric_limits<int>::max());
        allocatedMemoryInByte = memoryToAlloc;
        memoryBuffer = (unsigned char*)FAlignedMemory::AllocateBytes<32>(memoryToAlloc);
        FAssertLF(memoryBuffer);
        memset(memoryBuffer, 0, memoryToAlloc);

        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        blockCells          = reinterpret_cast<SymboleCellClass*>(memoryBuffer+sizeof(BlockHeader)+(blockIndexesTableSize*sizeof(int)));

        // Init header
        blockHeader->startingIndex = inStartingIndex;
        blockHeader->endingIndex   = inEndingIndex;
        blockHeader->numberOfCellsInBlock  = inNumberOfCells;
        blockHeader->blockIndexesTableSize = blockIndexesTableSize;

        cellMultipoles = (PoleCellClass*)FAlignedMemory::AllocateBytes<32>(inNumberOfCells*sizeof(PoleCellClass));
        cellLocals     = (LocalCellClass*)FAlignedMemory::AllocateBytes<32>(inNumberOfCells*sizeof(LocalCellClass));
        for(int idxCell = 0 ; idxCell < inNumberOfCells ; ++idxCell){
            new (&cellMultipoles[idxCell]) PoleCellClass();
            new (&cellLocals[idxCell]) LocalCellClass();
        }

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
                    (&blockCells[blockIndexesTable[idxCellPtr]])->~SymboleCellClass();
                    (&cellMultipoles[blockIndexesTable[idxCellPtr]])->~PoleCellClass();
                    (&cellLocals[blockIndexesTable[idxCellPtr]])->~LocalCellClass();
                }
            }
            FAlignedMemory::DeallocBytes(memoryBuffer);
            FAlignedMemory::DeallocBytes(cellMultipoles);
            FAlignedMemory::DeallocBytes(cellLocals);
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

    /** Give access to the buffer to send the data */
    const PoleCellClass* getRawMultipoleBuffer() const{
        return cellMultipoles;
    }

    /** The the size of the allocated buffer */
    int getMultipoleBufferSizeInByte() const {
        return sizeof(PoleCellClass)*blockHeader->numberOfCellsInBlock;
    }

    /** Give access to the buffer to send the data */
    const LocalCellClass* getRawLocalBuffer() const{
        return cellLocals;
    }

    /** The the size of the allocated buffer */
    int getLocalBufferSizeInByte() const {
        return sizeof(LocalCellClass)*blockHeader->numberOfCellsInBlock;
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
    CompositeCellClass getCompleteCell(const MortonIndex inIndex){
        FAssertLF(cellMultipoles && cellLocals);
        if( exists(inIndex) ){
            const int cellPos = blockIndexesTable[inIndex-blockHeader->startingIndex];
            return CompositeCellClass(&blockCells[cellPos], &cellMultipoles[cellPos], &cellLocals[cellPos]);
        }
        else return CompositeCellClass();
    }

    /** Return the address of the cell if it exists (or NULL) */
    CompositeCellClass getUpCell(const MortonIndex inIndex){
        FAssertLF(cellMultipoles);
        if( exists(inIndex) ){
            const int cellPos = blockIndexesTable[inIndex-blockHeader->startingIndex];
            return CompositeCellClass(&blockCells[cellPos], &cellMultipoles[cellPos], nullptr);
        }
        else return CompositeCellClass();
    }

    /** Return the address of the cell if it exists (or NULL) */
    CompositeCellClass getDownCell(const MortonIndex inIndex){
        FAssertLF(cellLocals);
        if( exists(inIndex) ){
            const int cellPos = blockIndexesTable[inIndex-blockHeader->startingIndex];
            return CompositeCellClass(&blockCells[cellPos], nullptr, &cellLocals[cellPos]);
        }
        else return CompositeCellClass();
    }

    /** Allocate a new cell by calling its constructor */
    template<typename... CellConstructorParams>
    void newCell(const MortonIndex inIndex, const int id, CellConstructorParams... args){
        FAssertLF(isInside(inIndex));
        FAssertLF(!exists(inIndex));
        FAssertLF(id < blockHeader->blockIndexesTableSize);
        new((void*)&blockCells[id]) SymboleCellClass(args...);
        blockIndexesTable[inIndex-blockHeader->startingIndex] = id;
    }

    /** Iterate on each allocated cells */
    template<typename... FunctionParams>
    void forEachCell(std::function<void(CompositeCellClass, FunctionParams...)> function, FunctionParams... args){
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
            if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
                const int cellPos = blockIndexesTable[idxCellPtr];
                function(CompositeCellClass(&blockCells[cellPos], &cellMultipoles[cellPos], &cellLocals[cellPos]), args...);
            }
        }
    }

    void forEachCell(std::function<void(CompositeCellClass)> function){
        for(int idxCellPtr = 0 ; idxCellPtr < blockHeader->blockIndexesTableSize ; ++idxCellPtr){
            if(blockIndexesTable[idxCellPtr] != CellIsEmptyFlag){
                const int cellPos = blockIndexesTable[idxCellPtr];
                function(CompositeCellClass(&blockCells[cellPos], &cellMultipoles[cellPos], &cellLocals[cellPos]));
            }
        }
    }
};


#endif // FGROUPOFCELLS_HPP
