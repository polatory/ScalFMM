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
#ifndef FEXTENDCELLTYPE_HPP
#define FEXTENDCELLTYPE_HPP

#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCellType
* This class is an extenssion.
* It proposes a target/source extenssion for cell.
* Because cells may have child that contains only
* sources or targets (in Tsm system) then it is important
* to not compute for nothing.
*/
class FExtendCellType {
protected:
    /** Particle potential type */
    static const int Neither = 0;
    static const int ContainsSrc = 1;
    static const int ContainsTargets = 2;

    /** Current type */
    int type;

public:
    /** Default constructor */
    FExtendCellType() : type(Neither) {
    }

    /** Copy constructor */
    FExtendCellType(const FExtendCellType& other) : type(other.type) {
    }

    /** Copy operator */
    FExtendCellType& operator=(const FExtendCellType& other) {
        this->type = other.type;
        return *this;
    }

    /** To know if a cell has sources */
    bool hasSrcChild() const {
        return this->type & ContainsSrc;
    }

    /** To know if a cell has targets */
    bool hasTargetsChild() const {
        return this->type & ContainsTargets;
    }

    /** To set cell as sources container */
    void setSrcChildTrue() {
        this->type |= ContainsSrc;
    }

    /** To set cell as targets container */
    void setTargetsChildTrue() {
        this->type |= ContainsTargets;
    }

public:
    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << type;
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> type;
    }
};


#endif //FEXTENDCELLTYPE_HPP


