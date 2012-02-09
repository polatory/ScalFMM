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
#ifndef FEXTENDFMBCELL_HPP
#define FEXTENDFMBCELL_HPP


#include <cstring>

#include "../Utils/FComplexe.hpp"

#include "FFmbKernels.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendFmbCell
* Please read the license.
*
* This class is an extenssion.
* It is needed by the Fmb Kernels.
*/
class FExtendFmbCell {
public:
    // FMB_Info_P is declared in FAbstractFmbKernels
    static const int MultipoleSize = int(((FMB_Info_P)+1) * ((FMB_Info_P)+2) * 0.5); //< The size of the multipole

protected:
    FComplexe multipole_exp[MultipoleSize]; //< For multipole extenssion
    FComplexe local_exp[MultipoleSize];     //< For local extenssion

public:
    /** Default constructor */
    FExtendFmbCell() {
    }

    /** Constructor */
    FExtendFmbCell(const FExtendFmbCell& other){
        (*this) = other;
    }

    /** Default destructor */
    virtual ~FExtendFmbCell(){
    }

    /** Copy constructor */
    FExtendFmbCell& operator=(const FExtendFmbCell& other) {
        memcpy(multipole_exp, other.multipole_exp, sizeof(FComplexe)*MultipoleSize);
        memcpy(local_exp, other.local_exp, sizeof(FComplexe)*MultipoleSize);
        return *this;
    }

    /** Get Multipole */
    const FComplexe* getMultipole() const {
        return this->multipole_exp;
    }
    /** Get Local */
    const FComplexe* getLocal() const {
        return this->local_exp;
    }

    /** Get Multipole */
    FComplexe* getMultipole() {
        return this->multipole_exp;
    }
    /** Get Local */
    FComplexe* getLocal() {
        return this->local_exp;
    }

};


#endif //FEXTENDFMBCELL_HPP


