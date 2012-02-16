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
#ifndef FEXTENDCELLTYPE_HPP
#define FEXTENDCELLTYPE_HPP

#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCellType
* Please read the license
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


