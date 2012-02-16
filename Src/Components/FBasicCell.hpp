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
#ifndef FBASICCELL_HPP
#define FBASICCELL_HPP


#include "../Extensions/FExtendPosition.hpp"
#include "../Extensions/FExtendMortonIndex.hpp"
#include "../Extensions/FExtendCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell
* Please read the license
*
* This class defines a basic cell used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractCell without
* inheriting from it.
*/
class FBasicCell : public FExtendMortonIndex, public FExtendCoordinate {
public:
    /** Default destructor */
    virtual ~FBasicCell(){
    }

    void save(FBufferWriter& buffer) const{
        FExtendMortonIndex::save(buffer);
        FExtendCoordinate::save(buffer);
    }
    void restore(FBufferReader& buffer){
        FExtendMortonIndex::restore(buffer);
        FExtendCoordinate::restore(buffer);
    }
};


#endif //FBASICCELL_HPP


