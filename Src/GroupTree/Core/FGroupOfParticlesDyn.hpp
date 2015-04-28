
// @SCALFMM_PRIVATE
#ifndef FGROUPOFPARTICLESDYN_HPP
#define FGROUPOFPARTICLESDYN_HPP


#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FAssert.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FAlignedMemory.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

#include <list>
#include <functional>


/**
* @brief The FGroupOfParticlesDyn class manages the leaves in block allocation.
*/
class FGroupOfParticlesDyn {
    /** One header is allocated at the beginning of each block */
    struct alignas(FStarPUDefaultAlign::StructAlign) BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int numberOfLeavesInBlock;
        int blockIndexesTableSize;
        size_t offsetSymbPart;
    };

    /** Information about a leaf */
    struct alignas(FStarPUDefaultAlign::StructAlign) LeafHeader {
        size_t offSetSymb;
        size_t offSetDown;
        size_t sizeSymb;
        size_t sizeDown;
    };


protected:
    static const int MemoryAlignementBytes     = FP2PDefaultAlignement;

    //< This value is for not used leaves
    static const int LeafIsEmptyFlag = -1;

    //< The size of memoryBuffer in byte
    size_t allocatedMemoryInByteSymb;
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    int*            blockIndexesTable;
    //< Pointer to leaves information
    LeafHeader*     leafHeader;

    //< Pointers to particle position x, y, z
    unsigned char* symbPart;

    //< Pointers to the particles data inside the block memory
    unsigned char* downPart;
    size_t allocatedMemoryInByteDown;

    /** To know if we have to delete the buffer */
    bool deleteBuffer;

public:
    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    FGroupOfParticlesDyn(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                      unsigned char* inAttributes)
        : allocatedMemoryInByteSymb(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), leafHeader(nullptr), symbPart(nullptr),
          downPart(nullptr), allocatedMemoryInByteDown(0), deleteBuffer(false){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(inBuffer);
        inBuffer += sizeof(BlockHeader);
        blockIndexesTable   = reinterpret_cast<int*>(inBuffer);
        inBuffer += (blockHeader->blockIndexesTableSize*sizeof(int));
        leafHeader          = reinterpret_cast<LeafHeader*>(inBuffer);
        inBuffer += (sizeof(LeafHeader)*blockHeader->numberOfLeavesInBlock);
        symbPart            = reinterpret_cast<unsigned char*>(blockHeader)+blockHeader->offsetSymbPart;

        downPart = inAttributes;
    }

    /**
 * @brief FGroupOfParticlesDyn
 * @param inStartingIndex first leaf morton index
 * @param inEndingIndex last leaf morton index + 1
 * @param inNumberOfLeaves total number of leaves in the interval (should be <= inEndingIndex-inEndingIndex)
 */
    FGroupOfParticlesDyn(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfLeaves,
                         const size_t sizePerLeafSymb[], const size_t sizePerLeafDown[])
        : allocatedMemoryInByteSymb(0), memoryBuffer(nullptr), blockHeader(nullptr), blockIndexesTable(nullptr), leafHeader(nullptr),
          symbPart(nullptr), downPart(nullptr), allocatedMemoryInByteDown(0), deleteBuffer(true){

        size_t totalSymb = 0;
        size_t totalDown = 0;
        for(int idxLeaf = 0 ; idxLeaf < inNumberOfLeaves ; ++idxLeaf){
            totalSymb += sizePerLeafSymb[idxLeaf];
            totalDown += sizePerLeafDown[idxLeaf];
        }

        // Find the number of leaf to allocate in the blocks
        const int blockIndexesTableSize = int(inEndingIndex-inStartingIndex);
        FAssertLF(inNumberOfLeaves <= blockIndexesTableSize);
        // Total number of bytes in the block
        const size_t memoryToAllocSymb = sizeof(BlockHeader)
                                    + (blockIndexesTableSize*sizeof(int))
                                    + (inNumberOfLeaves*sizeof(LeafHeader))
                                    + (totalSymb+MemoryAlignementBytes-1);

        // Allocate
        allocatedMemoryInByteSymb = memoryToAllocSymb;
        memoryBuffer = (unsigned char*)FAlignedMemory::AllocateBytes<MemoryAlignementBytes>(memoryToAllocSymb);
        FAssertLF(memoryBuffer);
        memset(memoryBuffer, 0, memoryToAllocSymb);

        allocatedMemoryInByteDown = totalDown;
        downPart = (unsigned char*)FAlignedMemory::AllocateBytes<MemoryAlignementBytes>(allocatedMemoryInByteDown);
        FAssertLF(downPart);
        memset(downPart, 0, allocatedMemoryInByteDown);

        // Move the pointers to the correct position
        unsigned char* bufferPtr = memoryBuffer;
        blockHeader         = reinterpret_cast<BlockHeader*>(bufferPtr);
        bufferPtr += sizeof(BlockHeader);
        blockIndexesTable   = reinterpret_cast<int*>(bufferPtr);
        bufferPtr += (blockIndexesTableSize*sizeof(int));
        leafHeader          = reinterpret_cast<LeafHeader*>(bufferPtr);
        bufferPtr += (inNumberOfLeaves*sizeof(LeafHeader));
        symbPart            = reinterpret_cast<unsigned char*>(size_t(bufferPtr+MemoryAlignementBytes-1) & ~size_t(MemoryAlignementBytes-1));

        FAssertLF(size_t(bufferPtr-memoryBuffer) <= allocatedMemoryInByteSymb);

        // Init header
        blockHeader->startingIndex = inStartingIndex;
        blockHeader->endingIndex   = inEndingIndex;
        blockHeader->numberOfLeavesInBlock  = inNumberOfLeaves;
        blockHeader->blockIndexesTableSize  = blockIndexesTableSize;
        blockHeader->offsetSymbPart = size_t(symbPart)-size_t(blockHeader);

        if(inNumberOfLeaves){
            leafHeader[0].sizeSymb = sizePerLeafSymb[0];
            leafHeader[0].sizeDown = sizePerLeafDown[0];
            leafHeader[0].offSetSymb = 0;
            leafHeader[0].offSetDown = 0;

            for(int idxLeaf = 1 ; idxLeaf < inNumberOfLeaves ; ++idxLeaf){
                leafHeader[idxLeaf].sizeSymb = sizePerLeafSymb[idxLeaf];
                leafHeader[idxLeaf].sizeDown = sizePerLeafDown[idxLeaf];
                leafHeader[idxLeaf].offSetSymb = leafHeader[idxLeaf-1].sizeSymb + leafHeader[idxLeaf-1].offSetSymb;
                leafHeader[idxLeaf].offSetDown = leafHeader[idxLeaf-1].sizeDown + leafHeader[idxLeaf-1].offSetDown;
            }
            FAssertLF(leafHeader[inNumberOfLeaves-1].offSetSymb+leafHeader[inNumberOfLeaves-1].sizeSymb == totalSymb);
            FAssertLF(leafHeader[inNumberOfLeaves-1].offSetDown+leafHeader[inNumberOfLeaves-1].sizeDown == totalDown);
        }

        // Set all index to not used
        for(int idxLeafPtr = 0 ; idxLeafPtr < blockIndexesTableSize ; ++idxLeafPtr){
            blockIndexesTable[idxLeafPtr] = LeafIsEmptyFlag;
        }
    }

    /** Call the destructor of leaves and dealloc block memory */
    ~FGroupOfParticlesDyn(){
        if(deleteBuffer){
            FAlignedMemory::DeallocBytes(memoryBuffer);
            FAlignedMemory::DeallocBytes(downPart);
        }
    }

    /** Give access to the buffer to send the data */
    const unsigned char* getRawBuffer() const{
        return memoryBuffer;
    }

    /** The the size of the allocated buffer */
    size_t getBufferSizeInByte() const {
        return allocatedMemoryInByteSymb;
    }

    /** Give access to the buffer to send the data */
    const unsigned char* getRawAttributesBuffer() const{
        return downPart;
    }

    /** Give access to the buffer to send the data */
    unsigned char* getRawAttributesBuffer(){
        return downPart;
    }

    /** The the size of the allocated buffer */
    size_t getAttributesBufferSizeInByte() const {
        return allocatedMemoryInByteDown;
    }

    /** To know if the object will delete the memory block */
    bool getDeleteMemory() const{
        return deleteBuffer;
    }

    /** The index of the fist leaf (set from the constructor) */
    MortonIndex getStartingIndex() const {
        return blockHeader->startingIndex;
    }

    /** The index of the last leaf + 1 (set from the constructor) */
    MortonIndex getEndingIndex() const {
        return blockHeader->endingIndex;
    }

    /** The number of leaf (set from the constructor) */
    int getNumberOfLeavesInBlock() const {
        return blockHeader->numberOfLeavesInBlock;
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
        return isInside(inIndex) && (blockIndexesTable[inIndex-blockHeader->startingIndex] != LeafIsEmptyFlag);
    }

    /** Allocate a new leaf by calling its constructor */
    void newLeaf(const MortonIndex inIndex, const int id){
        FAssertLF(isInside(inIndex));
        FAssertLF(!exists(inIndex));
        FAssertLF(id < blockHeader->blockIndexesTableSize);
        blockIndexesTable[inIndex-blockHeader->startingIndex] = id;
    }

    /** Iterate on each allocated leaves */
    template<class ParticlesAttachedClass>
    void forEachLeaf(std::function<void(ParticlesAttachedClass*)> function){
        for(int idxLeafPtr = 0 ; idxLeafPtr < blockHeader->blockIndexesTableSize ; ++idxLeafPtr){
            if(blockIndexesTable[idxLeafPtr] != LeafIsEmptyFlag){
                const int id = blockIndexesTable[idxLeafPtr];
                ParticlesAttachedClass leaf( (leafHeader[id].sizeSymb? symbPart + leafHeader[id].offSetSymb : nullptr),
                                             (downPart && leafHeader[id].sizeDown ?downPart + leafHeader[id].offSetDown : nullptr) );
                function(&leaf);
            }
        }
    }


    /** Return the address of the leaf if it exists (or NULL) */
    template<class ParticlesAttachedClass>
    ParticlesAttachedClass getLeaf(const MortonIndex leafIndex){
        if(blockIndexesTable[leafIndex - blockHeader->startingIndex] != LeafIsEmptyFlag){
            const int id = blockIndexesTable[leafIndex - blockHeader->startingIndex];
            return ParticlesAttachedClass((leafHeader[id].sizeSymb? symbPart + leafHeader[id].offSetSymb : nullptr),
                                          (downPart && leafHeader[id].sizeDown ?downPart + leafHeader[id].offSetDown : nullptr) );
        }
        return ParticlesAttachedClass();
    }

    /** Return the buffer for a leaf or null if it does not exist */
    unsigned char* getLeafSymbBuffer(const MortonIndex leafIndex){
        if(blockIndexesTable[leafIndex - blockHeader->startingIndex] != LeafIsEmptyFlag){
            const int id = blockIndexesTable[leafIndex - blockHeader->startingIndex];
            if(leafHeader[id].sizeSymb){
                return (symbPart + leafHeader[id].offSetSymb);
            }
        }
        return nullptr;
    }

    /** Return the buffer for a leaf or null if it does not exist */
    unsigned char* getLeafDownBuffer(const MortonIndex leafIndex){
        if(blockIndexesTable[leafIndex - blockHeader->startingIndex] != LeafIsEmptyFlag){
            const int id = blockIndexesTable[leafIndex - blockHeader->startingIndex];
            if(leafHeader[id].sizeDown){
                return (downPart?downPart + leafHeader[id].offSetDown : nullptr);
            }
        }
        return nullptr;
    }
};

#endif // FGROUPOFPARTICLESDYN_HPP

