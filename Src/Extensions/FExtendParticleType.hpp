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
#ifndef FEXTENDPARTICLETYPE_HPP
#define FEXTENDPARTICLETYPE_HPP

#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendParticleType
* Please read the license
* This class is an extenssion.
* It proposes a target/source extenssion for particles.
*/
class FExtendParticleType {
protected:
    /** Particle potential type */
    enum Type {
        Source,
        Target,
        Undef
    };

    /** Current type */
    Type type;

public:
    /** Default constructor */
    FExtendParticleType(const Type inType = Undef) : type(inType) {
    }

    /** Copy constructor */
    FExtendParticleType(const FExtendParticleType& other) : type(other.type) {
    }

    /** Copy operator */
    FExtendParticleType& operator=(const FExtendParticleType& other) {
        this->type = other.type;
        return *this;
    }

    /** To get the type */
    Type getParticleType() const {
        return this->type;
    }

    /** To set the type */
    void setParticleType(const Type inType) {
        this->type = inType;
    }

    /** To know if a particle is a target */
    bool isTarget() const{
        return this->type == Target;
    }

    /** To know if a particle is a source */
    bool isSource() const{
        return this->type == Source;
    }

    /** To know if a particle has an undefined type */
    bool isUndefinedType() const{
        return this->type == Undef;
    }

    /** To know if a particle is a target */
    void setAsTarget() {
        this->type = Target;
    }

    /** To know if a particle is a source */
    void setAsSource() {
        this->type = Source;
    }

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << type;
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> type;
    }
};


#endif //FEXTENDPARTICLETYPE_HPP


