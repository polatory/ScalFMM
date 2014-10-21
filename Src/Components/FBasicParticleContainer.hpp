// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FBASICPARTICLECONTAINER_HPP
#define FBASICPARTICLECONTAINER_HPP

#include "FAbstractParticleContainer.hpp"
#include "FAbstractSerializable.hpp"

#include "Utils/FAlignedMemory.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FPoint.hpp"
#include "FParticleType.hpp"

#include <array>

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FBasicParticle
 * Please read the license
 *
 * This class defines a container which can holds one type (AttributeClass)
 * for each particle.
 * The memory is allocated for all informations, the positions and the
 * request type.
 * For example is one want to store a struct for each particle:
 * @code
 * @code struct AStruct{
 * @code ...
 * @code };
 * @code FBasicParticleContainer<1, AStruct> container;
 * And then the access is done using:
 * @code AStruct* strucs = container.getAttributes<0>();
 */
template <unsigned NbAttributesPerParticle, class AttributeClass = FReal >
class FBasicParticleContainer : public FAbstractParticleContainer, public FAbstractSerializable {
protected:
    /** The number of particles in the container */
    int nbParticles;
    /** 3 pointers to 3 arrays of real to store the position */
    FReal* positions[3];
    /** The attributes requested by the user */
    AttributeClass* attributes[NbAttributesPerParticle];

    /** The allocated memory */
    int allocatedParticles;

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    /** Ending call for pushing the attributes */
    template<int index>
    void addParticleValue(const int /*insertPosition*/){
    }

    /** Ending call for pushing array of attributes */
    template<int index>
    void addParticleValueS(const int /*insertPosition*/,const int /*nbParticles*/){
    }


    /** Filling call for each attributes values */
    template<int index, typename... Args>
    void addParticleValue(const int insertPosition, const AttributeClass value, Args... args){
        // Compile test to ensure indexing
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        // insert the value
        attributes[index][insertPosition] = value;
        // Continue for reamining values
        addParticleValue<index+1>( insertPosition, args...);
    }

    /** Filling call for each attributes values
     * add multiples attributes from arrays
     */
    template<int index, typename... Args>
    void addParticleValueS(const int insertPosition, const int nbParts, const AttributeClass* value, Args... args){
        // Compile test to ensure indexing
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        for(int idxPart = 0; idxPart<nbParts ; ++idxPart){
            // insert the value
            attributes[index][insertPosition+idxPart] = value[idxPart];
            // Continue for reamining values
        }
        addParticleValueS<index+1>( insertPosition, nbParts, args...);

    }


    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    void increaseSizeIfNeeded(FSize sizeInput = 1){
        if( nbParticles+(sizeInput-1) >= allocatedParticles ){
            // allocate memory
            const int moduloParticlesNumber = (32/sizeof(FReal)); // We want to be rounded to 32B
            allocatedParticles = (FMath::Max(10,int(FReal(nbParticles+sizeInput)*1.5)) + moduloParticlesNumber - 1) & ~(moduloParticlesNumber-1);
            // init with 0
            const size_t allocatedBytes = (sizeof(FReal)*3 + sizeof(AttributeClass)*NbAttributesPerParticle)*allocatedParticles;
            FReal* newData  = reinterpret_cast<FReal*>(FAlignedMemory::Allocate32BAligned(allocatedBytes));
            memset( newData, 0, allocatedBytes);
            // copy memory
            const char*const toDelete  = reinterpret_cast<const char*>(positions[0]);
            for(int idx = 0 ; idx < 3 ; ++idx){
                memcpy(newData + (allocatedParticles * idx), positions[idx], sizeof(FReal) * nbParticles);
                positions[idx] = newData + (allocatedParticles * idx);
            }
            // copy attributes
            AttributeClass* startAddress = reinterpret_cast<AttributeClass*>(positions[2] + allocatedParticles);
            for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
                memcpy(startAddress + (allocatedParticles * idx), attributes[idx], sizeof(AttributeClass) * nbParticles);
                attributes[idx] = startAddress + (idx * allocatedParticles);
            }
            // delete old
            FAlignedMemory::Dealloc32BAligned(toDelete);
        }
    }

public:
    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    FBasicParticleContainer(const FBasicParticleContainer&)            = delete;
    FBasicParticleContainer& operator=(const FBasicParticleContainer&) = delete;

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    /** Basic contructor */
    FBasicParticleContainer() : nbParticles(0), allocatedParticles(0){
        memset(positions, 0, sizeof(positions[0]) * 3);
        memset(attributes, 0, sizeof(attributes[0]) * NbAttributesPerParticle);
    }

    /** Simply dalloc the memory using first pointer
   */
    ~FBasicParticleContainer(){
        FAlignedMemory::Dealloc32BAligned(positions[0]);
    }

    /**
   * @brief getNbParticles
   * @return  the number of particles
   */
    int getNbParticles() const{
        return nbParticles;
    }
    /**
   * @brief reset the number of particles
   * @warning Only the number of particles is set to 0, the particles are still here.
   */
    void resetNumberOfParticles()
    {
        nbParticles = 0 ;
    }
    /**
   * @brief getPositions
   * @return a FReal*[3] to get access to the positions
   */
    const FReal*const* getPositions() const {
        return positions;
    }

    /**
   * @brief getWPositions
   * @return get the position in write mode
   */
    FReal* const* getWPositions() {
        return positions;
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

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////


    /**
   * Push multiple particles
   * Should have a particle position fallowed by attributes
   * @param Array of position, number of parts to insert, followed by array of attribute
   */
    template<typename... Args>
    void pushArray(const FPoint * inParticlePosition, FSize numberOfParts, Args... args){
        const int positionToInsert = nbParticles;
        //Tests if enough space
        increaseSizeIfNeeded(numberOfParts);
        for(int idxPart = 0; idxPart<numberOfParts ; ++idxPart){
            // insert particle data
            positions[0][positionToInsert + idxPart] = inParticlePosition[idxPart].getX();
            positions[1][positionToInsert + idxPart] = inParticlePosition[idxPart].getY();
            positions[2][positionToInsert + idxPart] = inParticlePosition[idxPart].getZ();

        }
        // insert attribute data
        addParticleValueS<0>( nbParticles, numberOfParts ,args...);
        nbParticles += numberOfParts;

    }

    /**
     * Allocate the data to store requiereNbParticles particles
     */
    void reserve(const FSize requiereNbParticles){
        increaseSizeIfNeeded(requiereNbParticles - nbParticles);
    }


    /**
   * Push called bu FSimpleLeaf
   * Should have a particle position fallowed by attributes
   */
    template<typename... Args>
    void push(const FPoint& inParticlePosition, Args... args){
        // enought space?
        increaseSizeIfNeeded();

        // insert particle data
        positions[0][nbParticles] = inParticlePosition.getX();
        positions[1][nbParticles] = inParticlePosition.getY();
        positions[2][nbParticles] = inParticlePosition.getZ();
        // insert attribute data
        addParticleValue<0>( nbParticles, args...);
        nbParticles += 1;
    }


    /**
   * Push called bu FSimpleLeaf
   * Should have a particle position fallowed by attributes
   */
    template<typename... Args>
    void push(const FPoint& inParticlePosition, const std::array<AttributeClass , NbAttributesPerParticle>& values){
        // enought space?
        increaseSizeIfNeeded();

        // insert particle data
        positions[0][nbParticles] = inParticlePosition.getX();
        positions[1][nbParticles] = inParticlePosition.getY();
        positions[2][nbParticles] = inParticlePosition.getZ();
        // insert attribute data
        for(unsigned idxVal = 0 ; idxVal < NbAttributesPerParticle ; ++idxVal){
            attributes[idxVal][nbParticles] = values[idxVal];
        }
        nbParticles += 1;
    }

    /**
   * Push called usually by FTypedLeaf with the isTarget flag in addition
   */
    template<typename... Args>
    void push(const FPoint& inParticlePosition, const FParticleType /*particleType*/, Args... args){
        push(inParticlePosition, args...);
    }

    /** set nb particles to 0 */
    void clear(){
        nbParticles = 0;
    }

    /** to enable rearranging
   * indexesToRemove must be sorted
   * it removes all the particles at position indexesToRemove
   */
    void removeParticles(const int indexesToRemove[], const int nbParticlesToRemove){
        int offset = 1;
        int idxIndexes = 1;
        int idxIns = indexesToRemove[0] + 1;
        for( ; idxIns < nbParticles && idxIndexes < nbParticlesToRemove ; ++idxIns){
            if( idxIns == indexesToRemove[idxIndexes] ){
                idxIndexes += 1;
                offset += 1;
            }
            else{
                for(int idx = 0 ; idx < 3 ; ++idx){
                    positions[idx][idxIns-offset] = positions[idx][idxIns];
                }
                for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
                    attributes[idx][idxIns-offset] = attributes[idx][idxIns];
                }
            }
        }
        for( ; idxIns < nbParticles ; ++idxIns){
            for(int idx = 0 ; idx < 3 ; ++idx){
                positions[idx][idxIns-offset] = positions[idx][idxIns];
            }
            for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
                attributes[idx][idxIns-offset] = attributes[idx][idxIns];
            }
        }
        nbParticles -= nbParticlesToRemove;
    }

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    AttributeClass* getRawData(){
        return reinterpret_cast<AttributeClass*>(positions[2] + allocatedParticles);
    }

    const AttributeClass* getRawData() const {
        return reinterpret_cast<AttributeClass*>(positions[2] + allocatedParticles);
    }

    int getLeadingRawData() const {
        return allocatedParticles;
    }

    /////////////////////////////////////////////////////
    /////////////////////////////////////////////////////

    /** The size to send a leaf */
    int getSavedSize() const{
        return int(sizeof(nbParticles) + nbParticles * (3 * sizeof(FReal) + NbAttributesPerParticle * sizeof(AttributeClass)));
    }

    /** Save the current cell in a buffer */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        buffer << nbParticles;
        for(int idx = 0 ; idx < 3 ; ++idx){
            buffer.write(positions[idx], nbParticles);
        }
        for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
            buffer.write(attributes[idx], nbParticles);
        }
    }
    /** Restore the current cell from a buffer */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        buffer >> nbParticles;
        if( nbParticles >= allocatedParticles ){
            // allocate memory
            const int moduloParticlesNumber = (32/sizeof(FReal)); // We want to be rounded to 32B
            allocatedParticles = (FMath::Max(10,int(FReal(nbParticles+1)*1.5)) + moduloParticlesNumber - 1) & ~(moduloParticlesNumber-1);
            // init with 0
            const size_t allocatedBytes = (sizeof(FReal)*3 + sizeof(AttributeClass)*NbAttributesPerParticle)*allocatedParticles;
            FReal* newData  = reinterpret_cast<FReal*>(FAlignedMemory::Allocate32BAligned(allocatedBytes));
            memset( newData, 0, allocatedBytes);

            FAlignedMemory::Dealloc32BAligned(positions[0]);
            for(int idx = 0 ; idx < 3 ; ++idx){
                positions[idx] = newData + (allocatedParticles * idx);
            }
            AttributeClass* startAddress = reinterpret_cast<AttributeClass*>(positions[2] + allocatedParticles);
            for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
                attributes[idx] = startAddress + (idx * allocatedParticles);
            }
        }
        for(int idx = 0 ; idx < 3 ; ++idx){
            buffer.fillArray(positions[idx], nbParticles);
        }
        for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
            buffer.fillArray(attributes[idx], nbParticles);
        }
    }

    /** Reset the attributes to zeros */
    void resetToInitialState(){
        for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
            memset(attributes[idx], 0, sizeof(AttributeClass) * allocatedParticles);
        }
    }

    /** Reset the attributes to zeros */
    void resetToInitialState(const int idxAttribute){
        memset(attributes[idxAttribute], 0, sizeof(AttributeClass) * allocatedParticles);
    }
};


#endif //FBASICPARTICLECONTAINER_HPP
