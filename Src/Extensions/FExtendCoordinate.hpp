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
#ifndef FEXTENDCOORDINATE_HPP
#define FEXTENDCOORDINATE_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"
#include "../Containers/FMpiBufferReader.hpp"
#include "../Containers/FMpiBufferWriter.hpp"

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
  
  /////////////////////////////////////////////////////
  ///////////////  Test with FMpiBuffer*  /////////////
  /////////////////////////////////////////////////////
  
  /** Save current object */
  void save(FMpiBufferWriter& buffer) const {
    coordinate.save(buffer);
  }
  /** Retrieve current object */
  void restore(FMpiBufferReader& buffer) {
    coordinate.restore(buffer);
  }
  
};


#endif //FEXTENDCOORDINATE_HPP


