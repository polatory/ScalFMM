// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
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


