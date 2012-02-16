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
#ifndef FBASICPARTICLE_HPP
#define FBASICPARTICLE_HPP


#include "../Extensions/FExtendPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticle
* Please read the license
*
* This class defines a basic particle used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractParticle without
* inheriting from it.
*/
class FBasicParticle : public FExtendPosition{
public:
    void save(FBufferWriter& buffer) const{
        FExtendPosition::save(buffer);
    }
    void restore(FBufferReader& buffer){
        FExtendPosition::restore(buffer);
    }
};


#endif //FBASICPARTICLE_HPP


