// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FEXTENDPARTICLETYPE_HPP
#define FEXTENDPARTICLETYPE_HPP



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendParticleType
* Please read the license
* This class is an extenssion.
* It proposes a target/source extenssion.
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
    FExtendParticleType() : type(Undef) {
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

};


#endif //FEXTENDPARTICLETYPE_HPP


