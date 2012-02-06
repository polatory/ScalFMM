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
#ifndef FEXTENDPOTENTIAL_HPP
#define FEXTENDPOTENTIAL_HPP


#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPotential
* Please read the license
* This class is an extenssion.
* It proposes a Potential (FReal).
*/
class FExtendPotential {
protected:
    FReal potential;   //< The potential extended

public:
    /** Default constructor */
    FExtendPotential() : potential(0) {
    }

    /** Copy constructor */
    FExtendPotential(const FExtendPotential& other) : potential(other.potential) {
    }

    /** Copy operator */
    FExtendPotential& operator=(const FExtendPotential& other) {
        this->potential = other.potential;
        return *this;
    }

    /** To get the potential */
    FReal getPotential() const {
        return this->potential;
    }

    /** To set the potential */
    void setPotential(const FReal inPotential) {
        this->potential = inPotential;
    }

    /** To inc the potential */
    void incPotential(const FReal inPotential) {
        this->potential += inPotential;
    }

};


#endif //FEXTENDPOTENTIAL_HPP


