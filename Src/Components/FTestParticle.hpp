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
#ifndef FTESTPARTICLE_HPP
#define FTESTPARTICLE_HPP


#include "FBasicParticle.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTestParticle
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know whit it is.
*
* Particles just need the data down.
*/
class FTestParticle : public FBasicParticle {
protected:
    // To store data during downard pass
    long long int dataDown;
public:
    FTestParticle(): dataDown(0){
    }

    /** Default destructor */
    virtual ~FTestParticle(){
    }

    /** Get the down data */
    long long int getDataDown() const {
        return this->dataDown;
    }

    /** Set down data */
    void setDataDown(const long long int inData){
        this->dataDown = inData;
    }
};


#endif //FTESTPARTICLE_HPP


