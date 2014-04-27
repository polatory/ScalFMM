#ifndef FGROUPOFPARTICLES_HPP
#define FGROUPOFPARTICLES_HPP


#include "../Utils/FAssert.hpp"
#include "../Containers/FTreeCoordinate.hpp"

#include <list>
#include <functional>


/**
* @brief The FGroupOfParticles class manages the leaves in block allocation.
*/
template <unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FGroupOfParticles {
    /** One header is allocated at the beginning of each block */
    struct BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int numberOfLeavesInBlock;
        int blockIndexesTableSize;
    };

    /** Information about a leaf */
    struct LeafHeader {
        int nbParticles;
        size_t offSet;
    };

protected:
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    int*            blockIndexesTable;
    //< Pointer to leaves information
    LeafHeader*     leafHeader;
    //< The total number of particles in the group
    const int nbParticlesInGroup;

    //< Pointers to particle position x, y, z
    FReal* particlePosition[3];
    //< Bytes difference/offset between position
    size_t positionOffset;

    //< Pointers to the particles data inside the block memory
    AttributeClass*      particleAttributes[NbAttributesPerParticle];
    //< Bytes difference/offset between attributes
    size_t attributeOffset;

    //< This value is for not used leaves
    static const int LeafIsEmptyFlag = -1;

public:
    /**
 * @brief FGroupOfParticles
 * @param inStartingIndex first leaf morton index
 * @param inEndingIndex last leaf morton index + 1
 * @param inNumberOfLeaves total number of leaves in the interval (should be <= inEndingIndex-inEndingIndex)
 */
    FGroupOfParticles(const MortonIndex inStartingIndex, const MortonIndex inEndingIndex, const int inNumberOfLeaves, const int inNbParticles)
        : memoryBuffer(0), blockHeader(0), blockIndexesTable(0), leafHeader(0), nbParticlesInGroup(inNbParticles),
          positionOffset(0), attributeOffset(0) {
        memset(particlePosition, 0, sizeof(particlePosition));
        memset(particleAttributes, 0, sizeof(particleAttributes));

        // Find the number of leaf to allocate in the blocks
        const int blockIndexesTableSize = int(inEndingIndex-inStartingIndex);
        FAssertLF(inNumberOfLeaves <= blockIndexesTableSize);
        // Total number of bytes in the block
        const size_t sizeOfOneParticle = (3*sizeof(FReal) + NbAttributesPerParticle*sizeof(AttributeClass));
        const size_t memoryToAlloc = sizeof(BlockHeader) + (blockIndexesTableSize*sizeof(int)) + (inNumberOfLeaves*sizeof(LeafHeader))
                + inNbParticles*sizeOfOneParticle;

        // Allocate
        memoryBuffer = new unsigned char[memoryToAlloc];
        FAssertLF(memoryBuffer);
        memset(memoryBuffer, 0, memoryToAlloc);

        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        leafHeader          = reinterpret_cast<LeafHeader*>(memoryBuffer+sizeof(BlockHeader)+(blockIndexesTableSize*sizeof(int)));

        // Init header
        blockHeader->startingIndex = inStartingIndex;
        blockHeader->endingIndex   = inEndingIndex;
        blockHeader->numberOfLeavesInBlock  = inNumberOfLeaves;
        blockHeader->blockIndexesTableSize = blockIndexesTableSize;

        // Init particle pointers
        positionOffset = (sizeof(FReal) * inNbParticles);
        particlePosition[0] = reinterpret_cast<FReal*>(leafHeader + inNumberOfLeaves);
        particlePosition[1] = (particlePosition[0] + inNbParticles);
        particlePosition[2] = (particlePosition[1] + inNbParticles);

        // Redirect pointer to data
        attributeOffset = (sizeof(AttributeClass) * inNbParticles);
        unsigned char* previousPointer = reinterpret_cast<unsigned char*>(particlePosition[2] + inNbParticles);
        for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            particleAttributes[idxAttribute] = reinterpret_cast<AttributeClass*>(previousPointer);
            previousPointer += sizeof(AttributeClass)*inNbParticles;
        }

        // Set all index to not used
        for(int idxLeafPtr = 0 ; idxLeafPtr < blockIndexesTableSize ; ++idxLeafPtr){
            blockIndexesTable[idxLeafPtr] = LeafIsEmptyFlag;
        }
    }

    /** Call the destructor of leaves and dealloc block memory */
    ~FGroupOfParticles(){
        delete[] memoryBuffer;
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

    /** Get the total number of particles in the group */
    int getNbParticlesInGroup() const {
        return nbParticlesInGroup;
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
    void newLeaf(const MortonIndex inIndex, const int id, const int nbParticles, const size_t offsetInGroup){
        FAssertLF(isInside(inIndex));
        FAssertLF(!exists(inIndex));
        FAssertLF(id < blockHeader->blockIndexesTableSize);
        blockIndexesTable[inIndex-blockHeader->startingIndex] = id;
        leafHeader[id].nbParticles = nbParticles;
        leafHeader[id].offSet = offsetInGroup;
    }

    /** Iterate on each allocated leaves */
    template<class ParticlesAttachedClass>
    void forEachLeaf(std::function<void(ParticlesAttachedClass*)> function){
        for(int idxLeafPtr = 0 ; idxLeafPtr < blockHeader->blockIndexesTableSize ; ++idxLeafPtr){
            if(blockIndexesTable[idxLeafPtr] != LeafIsEmptyFlag){
                const int id = blockIndexesTable[idxLeafPtr];
                ParticlesAttachedClass leaf(leafHeader[id].nbParticles,
                                            particlePosition[0] + leafHeader[id].offSet,
                        positionOffset,
                        particleAttributes[0] + leafHeader[id].offSet,
                        attributeOffset);
                function(&leaf);
            }
        }
    }


    /** Return the address of the leaf if it exists (or NULL) */
    template<class ParticlesAttachedClass>
    ParticlesAttachedClass getLeaf(const MortonIndex leafIndex){
        if(blockIndexesTable[leafIndex - blockHeader->startingIndex] != LeafIsEmptyFlag){
            const int id = blockIndexesTable[leafIndex - blockHeader->startingIndex];
            return ParticlesAttachedClass(leafHeader[id].nbParticles,
                                          particlePosition[0] + leafHeader[id].offSet,
                                            positionOffset,
                                            particleAttributes[0] + leafHeader[id].offSet,
                                            attributeOffset);
        }
        return ParticlesAttachedClass();
    }
};



#endif // FGROUPOFPARTICLES_HPP
