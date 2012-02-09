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
#ifndef FFmaPARTICLE_HPP
#define FFmaPARTICLE_HPP


#include "FBasicParticle.hpp"
#include "../Extensions/FExtendPhysicalValue.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaParticle
* Please read the license
*
* This class defines a particle for FMA loader.
* As defined in FFmaLoader it needs {FBasicParticle,FExtendPhysicalValue}
*/
class FFmaParticle : public FBasicParticle, public FExtendPhysicalValue {
public:
    /** Default destructor */
    virtual ~FFmaParticle(){
    }
};


#endif //FFmaPARTICLE_HPP


