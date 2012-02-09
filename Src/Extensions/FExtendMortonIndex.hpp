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

};


#endif //FEXTENDMORTONINDEX_HPP


