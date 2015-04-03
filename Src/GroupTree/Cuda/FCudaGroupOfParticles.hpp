// @SCALFMM_PRIVATE
#ifndef FCUDAGROUPOFPARTICLES_HPP
#define FCUDAGROUPOFPARTICLES_HPP

#include "FCudaGlobal.hpp"
#include "../../Utils/FGlobal.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

template <class FReal, unsigned NbSymbAttributes, unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FCudaGroupOfParticles {
    /** One header is allocated at the beginning of each block */
    struct alignas(FStarPUDefaultAlign::StructAlign) BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int numberOfLeavesInBlock;
        int blockIndexesTableSize;

        //< The real number of particles allocated
        FSize nbParticlesAllocatedInGroup;
        //< Starting point of position
        size_t offsetPosition;
        //< Bytes difference/offset between position
        size_t positionsLeadingDim;
        //< Bytes difference/offset between attributes
        size_t attributeLeadingDim;
        //< The total number of particles in the group
        FSize nbParticlesInGroup;
    };

    /** Information about a leaf */
    struct alignas(FStarPUDefaultAlign::StructAlign) LeafHeader {
        FSize nbParticles;
        size_t offSet;
    };


protected:
    static const int MemoryAlignementBytes     = FP2PDefaultAlignement;
    static const int MemoryAlignementParticles = MemoryAlignementBytes/sizeof(FReal);

    /** This function return the correct number of particles that should be used to have a correct pack.
     * If alignement is 32 and use double (so 4 particles in pack), then this function returns:
     * RoundToUpperParticles(1) = 1 + 3 = 4
     * RoundToUpperParticles(63) = 63 + 1 = 64
     */
    template <class NumClass>
    __device__ static NumClass RoundToUpperParticles(const NumClass& nbParticles){
        return nbParticles + (MemoryAlignementParticles - (nbParticles%MemoryAlignementParticles));
    }

    //< This value is for not used leaves
    static const int LeafIsEmptyFlag = -1;

    //< The size of memoryBuffer in byte
    int allocatedMemoryInByte;
    //< Pointer to a block memory
    unsigned char* memoryBuffer;

    //< Pointer to the header inside the block memory
    BlockHeader*    blockHeader;
    //< Pointer to the indexes table inside the block memory
    int*            blockIndexesTable;
    //< Pointer to leaves information
    LeafHeader*     leafHeader;
    //< The total number of particles in the group
    const FSize nbParticlesInGroup;

    //< Pointers to particle position x, y, z
    FReal* particlePosition[3];

    //< Pointers to the particles data inside the block memory
    AttributeClass* attributesBuffer;
    AttributeClass* particleAttributes[NbSymbAttributes+NbAttributesPerParticle];

public:
    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    __device__ FCudaGroupOfParticles(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte,
                                     unsigned char* inAttributes)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), leafHeader(nullptr), nbParticlesInGroup(0),
          attributesBuffer(nullptr){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        leafHeader          = reinterpret_cast<LeafHeader*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));

        // Init particle pointers
        // Assert blockHeader->positionsLeadingDim == (sizeof(FReal) * blockHeader->nbParticlesAllocatedInGroup);
        particlePosition[0] = reinterpret_cast<FReal*>(memoryBuffer + blockHeader->offsetPosition);
        particlePosition[1] = (particlePosition[0] + blockHeader->nbParticlesAllocatedInGroup);
        particlePosition[2] = (particlePosition[1] + blockHeader->nbParticlesAllocatedInGroup);

        // Redirect pointer to data
        // Assert blockHeader->attributeLeadingDim == (sizeof(AttributeClass) * blockHeader->nbParticlesAllocatedInGroup);
        AttributeClass* symAttributes = (AttributeClass*)(&particlePosition[2][blockHeader->nbParticlesAllocatedInGroup]);
        for(unsigned idxAttribute = 0 ; idxAttribute < NbSymbAttributes ; ++idxAttribute){
            particleAttributes[idxAttribute] = symAttributes;
            symAttributes += blockHeader->nbParticlesAllocatedInGroup;
        }
        if(inAttributes){
            attributesBuffer = (AttributeClass*)inAttributes;
            for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
                particleAttributes[idxAttribute+NbSymbAttributes] = &attributesBuffer[idxAttribute*blockHeader->nbParticlesAllocatedInGroup];
            }
        }
    }

    /** The index of the fist leaf (set from the constructor) */
    __device__ MortonIndex getStartingIndex() const {
        return blockHeader->startingIndex;
    }

    /** The index of the last leaf + 1 (set from the constructor) */
    __device__ MortonIndex getEndingIndex() const {
        return blockHeader->endingIndex;
    }

    /** The number of leaf (set from the constructor) */
    __device__ int getNumberOfLeavesInBlock() const {
        return blockHeader->numberOfLeavesInBlock;
    }

    /** Get the total number of particles in the group */
    __device__ FSize getNbParticlesInGroup() const {
        return nbParticlesInGroup;
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
        return isInside(inIndex) && (blockIndexesTable[inIndex-blockHeader->startingIndex] != LeafIsEmptyFlag);
    }

    /** Return the address of the leaf if it exists (or NULL) */
    template<class ParticlesAttachedClass>
    __device__ ParticlesAttachedClass getLeaf(const MortonIndex leafIndex){
        if(blockIndexesTable[leafIndex - blockHeader->startingIndex] != LeafIsEmptyFlag){
            const int id = blockIndexesTable[leafIndex - blockHeader->startingIndex];
            return ParticlesAttachedClass(leafHeader[id].nbParticles,
                                          particlePosition[0] + leafHeader[id].offSet,
                                            blockHeader->positionsLeadingDim,
                                            (attributesBuffer?particleAttributes[NbSymbAttributes] + leafHeader[id].offSet:nullptr),
                                            blockHeader->attributeLeadingDim);
        }
        return ParticlesAttachedClass();
    }
};

#endif // FCUDAGROUPOFPARTICLES_HPP

