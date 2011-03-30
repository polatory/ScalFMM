#ifndef FEXTENDVALUE_HPP
#define FEXTENDVALUE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendValue
* Please read the license
* This class is an extenssion.
* It proposes a value (FReal).
*/
class FExtendValue {
protected:
    FReal value;   //< A simple value

public:
    /** Default constructor */
    FExtendValue() : value(0) {
    }

    /** Copy constructor */
    FExtendValue(const FExtendValue& other) : value(other.value) {
    }

    /** Destructor */
    virtual ~FExtendValue(){
    }

    /** Copy Constructor */
    FExtendValue& operator=(const FExtendValue& other) {
        this->value = other.value;
        return *this;
    }

    /** To get the value */
    FReal getValue() const {
        return this->value;
    }

    /** To set the value */
    void setValue(const FReal invalue) {
        this->value = invalue;
    }

};


#endif //FEXTENDVALUE_HPP

// [--LICENSE--]
