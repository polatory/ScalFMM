// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#ifndef FEXTENDCELLTYPE_HPP
#define FEXTENDCELLTYPE_HPP

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
    /** Current type */
    bool containsTargets;
    bool containsSources;

public:
    /** Default constructor */
    FExtendCellType() : containsTargets(false), containsSources(false) {
    }

    /** Copy constructor */
    FExtendCellType(const FExtendCellType& other) : containsTargets(other.containsTargets),
            containsSources(other.containsSources){
    }

    /** Copy operator */
    FExtendCellType& operator=(const FExtendCellType& other) {
        this->containsTargets = other.containsTargets;
        this->containsSources = other.containsSources;
        return *this;
    }

    /** To know if a cell has sources */
    bool hasSrcChild() const {
        return containsSources;
    }

    /** To know if a cell has targets */
    bool hasTargetsChild() const {
        return containsTargets;
    }

    /** To set cell as sources container */
    void setSrcChildTrue() {
        containsSources = true;
    }

    /** To set cell as targets container */
    void setTargetsChildTrue() {
        containsTargets = true;
    }

public:
    /** Save current object */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const {
        buffer << containsTargets;
        buffer << containsSources;
    }
    /** Retrieve current object */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer) {
        buffer >> containsTargets;
        buffer >> containsSources;
    }
    /** reset to unknown type */
    void resetToInitialState(){
        containsTargets = false;
        containsSources = false;
    }

    FSize getSavedSize() const {
        return FSize(sizeof(containsTargets) + sizeof(containsSources));
    }
};


#endif //FEXTENDCELLTYPE_HPP


