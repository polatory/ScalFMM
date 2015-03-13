// @SCALFMM_PRIVATE
#ifndef FCUDAGROUPOFPARTICLES_HPP
#define FCUDAGROUPOFPARTICLES_HPP

#include "FCudaGlobal.hpp"


template <unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FCudaGroupOfParticles {
    /** One header is allocated at the beginning of each block */
    struct alignas(1) BlockHeader{
        MortonIndex startingIndex;
        MortonIndex endingIndex;
        int numberOfLeavesInBlock;
        int blockIndexesTableSize;

        //< The real number of particles allocated
        int nbParticlesAllocatedInGroup;
        //< Bytes difference/offset between position
        size_t positionOffset;
        //< Bytes difference/offset between attributes
        size_t attributeOffset;
        //< The total number of particles in the group
        int nbParticlesInGroup;
    };

    /** Information about a leaf */
    struct alignas(1) LeafHeader {
        int nbParticles;
        size_t offSet;
    };


protected:
    static const int MemoryAlignementBytes     = 32;
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
    const int nbParticlesInGroup;

    //< Pointers to particle position x, y, z
    FReal* particlePosition[3];

    //< Pointers to the particles data inside the block memory
    AttributeClass*      particleAttributes[NbAttributesPerParticle];

public:
    /**
     * Init from a given buffer
     * @param inBuffer
     * @param inAllocatedMemoryInByte
     */
    __device__ FCudaGroupOfParticles(unsigned char* inBuffer, const size_t inAllocatedMemoryInByte)
        : allocatedMemoryInByte(inAllocatedMemoryInByte), memoryBuffer(inBuffer),
          blockHeader(nullptr), blockIndexesTable(nullptr), leafHeader(nullptr), nbParticlesInGroup(0){
        // Move the pointers to the correct position
        blockHeader         = reinterpret_cast<BlockHeader*>(memoryBuffer);
        blockIndexesTable   = reinterpret_cast<int*>(memoryBuffer+sizeof(BlockHeader));
        leafHeader          = reinterpret_cast<LeafHeader*>(memoryBuffer+sizeof(BlockHeader)+(blockHeader->blockIndexesTableSize*sizeof(int)));

        // Init particle pointers
        blockHeader->positionOffset = (sizeof(FReal) * blockHeader->nbParticlesAllocatedInGroup);
        particlePosition[0] = reinterpret_cast<FReal*>((reinterpret_cast<size_t>(leafHeader + blockHeader->numberOfLeavesInBlock)
                                                       +MemoryAlignementBytes-1) & ~(MemoryAlignementBytes-1));
        particlePosition[1] = (particlePosition[0] + blockHeader->nbParticlesAllocatedInGroup);
        particlePosition[2] = (particlePosition[1] + blockHeader->nbParticlesAllocatedInGroup);

        // Redirect pointer to data
        blockHeader->attributeOffset = (sizeof(AttributeClass) * blockHeader->nbParticlesAllocatedInGroup);
        unsigned char* previousPointer = reinterpret_cast<unsigned char*>(particlePosition[2] + blockHeader->nbParticlesAllocatedInGroup);
        for(unsigned idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            particleAttributes[idxAttribute] = reinterpret_cast<AttributeClass*>(previousPointer);
            previousPointer += sizeof(AttributeClass)*blockHeader->nbParticlesAllocatedInGroup;
        }
    }

    /** Call the destructor of leaves and dealloc block memory */
    __device__ ~FCudaGroupOfParticles(){
    }

    /** Give access to the buffer to send the data */
    __device__ const unsigned char* getRawBuffer() const{
        return memoryBuffer;
    }

    /** The the size of the allocated buffer */
    __device__ int getBufferSizeInByte() const {
        return allocatedMemoryInByte;
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
    __device__ int getNbParticlesInGroup() const {
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
                                            blockHeader->positionOffset,
                                            particleAttributes[0] + leafHeader[id].offSet,
                                            blockHeader->attributeOffset);
        }
        return ParticlesAttachedClass();
    }
};

#endif // FCUDAGROUPOFPARTICLES_HPP

