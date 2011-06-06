#ifndef FEXTENDCOORDINATE_HPP
#define FEXTENDCOORDINATE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCoordinate
* Please read the license
* This class is an extenssion.
* It proposes a mortonIndex.
*/
class FExtendCoordinate {
protected:
    FTreeCoordinate coordinate; //< The position

public:
    /** Default constructor */
    FExtendCoordinate() {
    }

    /** Copy constructor */
    FExtendCoordinate(const FExtendCoordinate& other) : coordinate(other.coordinate) {
    }

    /** Destructor */
    virtual ~FExtendCoordinate(){
    }

    /** Copy operator */
    FExtendCoordinate& operator=(const FExtendCoordinate& other) {
        this->coordinate = other.coordinate;
        return *this;
    }

    /** To get the position */
    const FTreeCoordinate& getCoordinate() const {
        return this->coordinate;
    }

    /** To set the position */
    void setCoordinate(const FTreeCoordinate& inCoordinate) {
        this->coordinate = inCoordinate;
    }

    /** To set the position from 3 FReals */
    void setCoordinate(const long inX, const long inY, const long inZ) {
        this->coordinate.setX(inX);
        this->coordinate.setY(inY);
        this->coordinate.setZ(inZ);
    }

};


#endif //FEXTENDCOORDINATE_HPP

// [--LICENSE--]
