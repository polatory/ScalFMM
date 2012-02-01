#ifndef FEXTENDPOSITION_HPP
#define FEXTENDPOSITION_HPP
// [--License--]

#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPosition
* Please read the license
* This class is an extenssion.
* It proposes a mortonIndex.
*/
class FExtendPosition {
protected:
    F3DPosition position; //< The position

public:
    /** Default constructor */
    FExtendPosition() {
    }

    /** Copy constructor */
    FExtendPosition(const FExtendPosition& other) : position(other.position) {
    }

    /** Copy operator */
    FExtendPosition& operator=(const FExtendPosition& other) {
        this->position = other.position;
        return *this;
    }

    /** To get the position */
    const F3DPosition& getPosition() const {
        return this->position;
    }

    /** To set the position */
    void setPosition(const F3DPosition& inPosition) {
        this->position = inPosition;
    }

    /** To set the position from 3 FReals */
    void setPosition(const FReal inX, const FReal inY, const FReal inZ) {
        this->position.setX(inX);
        this->position.setY(inY);
        this->position.setZ(inZ);
    }

    /** Set Position */
    void incPosition(const F3DPosition& inPosition) {
        this->position += inPosition;
    }

    /** Set Position with 3 FReals */
    void incPosition(const FReal inPx, const FReal inPy, const FReal inPz) {
        this->position.incX(inPx);
        this->position.incY(inPy);
        this->position.incZ(inPz);
    }

};


#endif //FEXTENDPOSITION_HPP

// [--END--]
