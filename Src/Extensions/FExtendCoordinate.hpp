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
#ifndef FEXTENDCOORDINATE_HPP
#define FEXTENDCOORDINATE_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCoordinate
* Please read the license
* This class is an extenssion.
* It proposes a tree coordinate.
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
    void setCoordinate(const int inX, const int inY, const int inZ) {
        this->coordinate.setX(inX);
        this->coordinate.setY(inY);
        this->coordinate.setZ(inZ);
    }


    /** Save current object */
    void save(FBufferWriter& buffer) const {
        coordinate.save(buffer);
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        coordinate.restore(buffer);
    }
};


#endif //FEXTENDCOORDINATE_HPP


