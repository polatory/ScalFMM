#ifndef FEXTENDFMBCELL_HPP
#define FEXTENDFMBCELL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <string.h>

#include "../Utils/FComplexe.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendFmbCell
* Please read the license.
*
* This class is an extenssion.
* It is needed by the Fmb Kernels.
*/
template <int P>
class FExtendFmbCell {
protected:
    static const int FMB_Info_P = P;        //< P >> FMB_Info.P
    static const int MultipoleSize = int(((FMB_Info_P)+1) * ((FMB_Info_P)+2) * 0.5); //< The size of the multipole

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
        memcpy(multipole_exp, other.multipole_exp, sizeof(multipole_exp));
        memcpy(local_exp, other.local_exp, sizeof(local_exp));
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

// [--LICENSE--]
