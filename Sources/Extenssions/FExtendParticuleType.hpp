#ifndef FEXTENDPARTICULETYPE_HPP
#define FEXTENDPARTICULETYPE_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendParticuleType
* Please read the license
* This class is an extenssion.
* It proposes a target/source extenssion.
*/
class FExtendParticuleType {
protected:
    /** Particule potential type */
    enum Type {
        Source,
        Target,
        Undef
    };

    /** Current type */
    Type type;

public:
    /** Default constructor */
    FExtendParticuleType() : type(Undef) {
    }

    /** Copy constructor */
    FExtendParticuleType(const FExtendParticuleType& other) : type(other.type) {
    }

    /** Destructor */
    virtual ~FExtendParticuleType(){
    }

    /** Copy operator */
    FExtendParticuleType& operator=(const FExtendParticuleType& other) {
        this->type = other.type;
        return *this;
    }

    /** To get the type */
    Type getParticuleType() const {
        return this->type;
    }

    /** To set the type */
    void setParticuleType(const Type inType) {
        this->type = inType;
    }

    /** To know if a particule is a target */
    bool isTarget() const{
        return this->type == Target;
    }

    /** To know if a particule is a source */
    bool isSource() const{
        return this->type == Source;
    }

    /** To know if a particule has an undefined type */
    bool isUndefinedType() const{
        return this->type == Undef;
    }

    /** To know if a particule is a target */
    void setAsTarget() {
        this->type = Target;
    }

    /** To know if a particule is a source */
    void setAsSource() {
        this->type = Source;
    }

};


#endif //FEXTENDPARTICULETYPE_HPP

// [--LICENSE--]
