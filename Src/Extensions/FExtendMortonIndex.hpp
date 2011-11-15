#ifndef FEXTENDMORTONINDEX_HPP
#define FEXTENDMORTONINDEX_HPP
// /!\ Please, you must read the license at the bottom of this page

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

};


#endif //FEXTENDMORTONINDEX_HPP

// [--LICENSE--]
