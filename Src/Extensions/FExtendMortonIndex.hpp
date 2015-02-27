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
#ifndef FEXTENDMORTONINDEX_HPP
#define FEXTENDMORTONINDEX_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendMortonIndex
* Please read the license
* This class is an extenssion.
* It proposes a mortonIndex.
*/
class FExtendMortonIndex {
protected:
    MortonIndex mortonIndex;    //< Morton index (need by most elements)

public:
    /** Default constructor */
    FExtendMortonIndex() : mortonIndex(0) {
    }

    /** Copy constructor */
    FExtendMortonIndex(const FExtendMortonIndex& other) : mortonIndex(other.mortonIndex) {
    }

    /** Copy operator */
    FExtendMortonIndex& operator=(const FExtendMortonIndex& other) {
        this->mortonIndex = other.mortonIndex;
        return *this;
    }

    /** To get the morton index */
    MortonIndex getMortonIndex() const {
        return this->mortonIndex;
    }

    /** To set the morton index */
    void setMortonIndex(const MortonIndex inMortonIndex) {
        this->mortonIndex = inMortonIndex;
    }

    /** Save current object */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const {
        buffer << mortonIndex;
    }
    /** Retrieve current object */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer) {
        buffer >> mortonIndex;
    }

    int getSavedSize() const {
        return int(sizeof(mortonIndex));
    }
};


#endif //FEXTENDMORTONINDEX_HPP


