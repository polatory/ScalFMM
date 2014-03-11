#ifndef FGROUPATTACHEDLEAF_HPP
#define FGROUPATTACHEDLEAF_HPP

#include "../Utils/FGlobal.hpp"

/**
 * This class is "attached" to a buffer.
 * It is a wrapper/interface to a memory area that stores all the particles info.
 * The idea is to hidde the group allocation from the group tree but
 * to keep the same interface than the FBasicParticlesContainer.
 */
template <unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FGroupAttachedLeaf {
protected:
    //< Nb of particles in the current leaf
    const int nbParticles;
    //< Pointers to the positions of the particles
    FReal* positionsPointers[3];
    //< Pointers to the attributes of the particles
    AttributeClass* attributes[NbAttributesPerParticle];

    // Forbid copy even if there is no real reason to do that
    FGroupAttachedLeaf(const FGroupAttachedLeaf&) = delete;
    FGroupAttachedLeaf& operator=(const FGroupAttachedLeaf&) = delete;

public:
    /**
     * @brief FGroupAttachedLeaf
     * @param inNbParticles the number of particles in the leaf
     * @param inPositionBuffer the memory address of the X array of particls
     * @param inLeadingPosition each position is access by inPositionBuffer + in bytes inLeadingPosition*idx
     * @param inAttributesBuffer the memory address of the first attribute
     * @param inLeadingAttributes each attribute is access by inAttributesBuffer + in bytes inLeadingAttributes*idx
     */
    FGroupAttachedLeaf(const int inNbParticles, FReal* inPositionBuffer, const size_t inLeadingPosition,
                       AttributeClass* inAttributesBuffer, const size_t inLeadingAttributes)
        : nbParticles(inNbParticles){
        // Redirect pointers to position
        positionsPointers[0] = inPositionBuffer;
        positionsPointers[1] = reinterpret_cast<FReal*>(reinterpret_cast<unsigned char*>(inPositionBuffer) + inLeadingPosition);
        positionsPointers[2] = reinterpret_cast<FReal*>(reinterpret_cast<unsigned char*>(inPositionBuffer) + inLeadingPosition*2);

        // Redirect pointers to data
        for(int idxAttribute = 0 ; idxAttribute < NbAttributesPerParticle ; ++idxAttribute){
            particleAttributes[idxAttribute] = reinterpret_cast<AttributeClass*>(reinterpret_cast<unsigned char*>(inAttributesBuffer) + idxAttribute*inLeadingAttributes);
        }
    }

    /**
     * @brief getNbParticles
     * @return the number of particles in the leaf
     */
    int getNbParticles() const{
        return nbParticles;
    }

    /**
     * @brief getPositions
     * @return a FReal*[3] to get access to the positions
     */
    const FReal*const* getPositions() const {
        return positionsPointers;
    }

    /**
     * @brief getWPositions
     * @return get the position in write mode
     */
    FReal* const* getWPositions() {
        return positionsPointers;
    }

    /**
     * @brief getAttribute
     * @param index
     * @return the attribute at index index
     */
    AttributeClass* getAttribute(const int index) {
        return attributes[index];
    }

    /**
     * @brief getAttribute
     * @param index
     * @return
     */
    const AttributeClass* getAttribute(const int index) const {
        return attributes[index];
    }

    /**
     * Get the attribute with a forcing compile optimization
     */
    template <int index>
    AttributeClass* getAttribute() {
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        return attributes[index];
    }

    /**
     * Get the attribute with a forcing compile optimization
     */
    template <int index>
    const AttributeClass* getAttribute() const {
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        return attributes[index];
    }
};

#endif // FGROUPATTACHEDLEAF_HPP
