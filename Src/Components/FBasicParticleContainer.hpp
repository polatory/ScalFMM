// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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

#include <type_traits>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticle
* Please read the license
*
* This class defines a basic particle used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractParticle without
* inheriting from it.
*/
template <unsigned NbAttributesPerParticle, class AttributeClass = FReal >
class FBasicParticleContainer : public FAbstractParticleContainer {
protected:
    int nbParticles;
    FReal* positions[3];
    AttributeClass* attributes[NbAttributesPerParticle];

    int allocatedParticles;

    template<int index>
    void addParticleValue(const int /*insertPosition*/){
    }

    template<int index, typename... Args>
    void addParticleValue(const int insertPosition, const AttributeClass value, Args... args){
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        attributes[index][insertPosition] = value;
        addParticleValue<index+1>( insertPosition, args...);
    }

public:
    FBasicParticleContainer(const FBasicParticleContainer&) = delete;
    FBasicParticleContainer& operator=(const FBasicParticleContainer&) = delete;

    FBasicParticleContainer() : nbParticles(0), allocatedParticles(0){
        memset(positions, 0, sizeof(positions[0]) * 3);
        memset(attributes, 0, sizeof(attributes[0]) * NbAttributesPerParticle);
    }

    ~FBasicParticleContainer(){
        delete[] reinterpret_cast<char*>(positions[0]);
    }

    int getNbParticles() const{
        return nbParticles;
    }

    const FReal*const*const getPositions() const {
        return positions;
    }

    FReal*const*const getWPositions() {
        return positions;
    }

    AttributeClass* getAttribute(const int index) {
        return attributes[index];
    }

    const AttributeClass* getAttribute(const int index) const {
        return attributes[index];
    }

    template <int index>
    AttributeClass* getAttribute() {
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        return attributes[index];
    }

    template <int index>
    const AttributeClass* getAttribute() const {
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        return attributes[index];
    }

    template<typename... Args>
    void push(const FPoint& inParticlePosition, Args... args){
        if( nbParticles == allocatedParticles ){
            allocatedParticles = FMath::Max(10,int(FReal(nbParticles+1)*1.5));
            FReal* newData  = reinterpret_cast<FReal*>(new char[(sizeof(FReal)*3 + sizeof(AttributeClass)*NbAttributesPerParticle)*allocatedParticles]);
            memset( newData, 0, (sizeof(FReal)*3 + sizeof(AttributeClass)*NbAttributesPerParticle)*allocatedParticles);

            const char*const toDelete  = reinterpret_cast<const char*>(positions[0]);
            for(int idx = 0 ; idx < 3 ; ++idx){
                memcpy(newData + (allocatedParticles * idx), positions[idx], sizeof(FReal) * nbParticles);
                positions[idx] = newData + (allocatedParticles * idx);
            }
            AttributeClass* startAddress = reinterpret_cast<AttributeClass*>(positions[2] + allocatedParticles);
            for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
                memcpy(startAddress + (allocatedParticles * idx), attributes[idx], sizeof(AttributeClass) * nbParticles);
                attributes[idx] = startAddress + (idx * allocatedParticles);
            }
            delete[] toDelete;
        }
        positions[0][nbParticles] = inParticlePosition.getX();
        positions[1][nbParticles] = inParticlePosition.getY();
        positions[2][nbParticles] = inParticlePosition.getZ();
        addParticleValue<0>( nbParticles, args...);
        nbParticles += 1;
    }

    template<typename... Args>
    void push(const FPoint& inParticlePosition, const bool /*isTarget*/, Args... args){
        push(inParticlePosition, args...);
    }

    void clear(){
        nbParticles = 0;
    }

    /** Save the current cell in a buffer */
    void save(FBufferWriter& buffer) const{
        buffer << nbParticles;
        for(int idx = 0 ; idx < 3 ; ++idx){
            buffer.write(positions[idx], nbParticles);
        }
        for(unsigned idx = 0 ; idx < NbAttributesPerParticle ; ++idx){
            buffer.write(attributes[idx], nbParticles);
        }
    }
    /** Restore the current cell from a buffer */
    void restore(FBufferReader& buffer){
        buffer >> nbParticles;
        if( nbParticles >= allocatedParticles ){
            allocatedParticles = FMath::Max(10,int(FReal(nbParticles+1)*1.5));
            FReal* newData  = reinterpret_cast<FReal*>(new char[(sizeof(FReal)*3 + sizeof(AttributeClass)*NbAttributesPerParticle)*allocatedParticles]);
            memset( newData, 0, (sizeof(FReal)*3 + sizeof(AttributeClass)*NbAttributesPerParticle)*allocatedParticles);

            delete[] reinterpret_cast<char*>(positions[0]);
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

    /** to enable rearranging */
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

};


#endif //FBASICPARTICLECONTAINER_HPP


