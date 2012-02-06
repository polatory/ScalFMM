// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
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


