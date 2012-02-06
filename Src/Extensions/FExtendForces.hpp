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
#ifndef FEXTENDFORCES_HPP
#define FEXTENDFORCES_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendForces
* Please read the license
*
* This class is an extenssion.
* It proposes a 3d array as a forces vector.
*/
class FExtendForces {
protected:
    F3DPosition forces; //< 3D vector stored in a position object

public:
    /** Default constructor */
    FExtendForces() {
    }

    /** Copy constructor */
    FExtendForces(const FExtendForces& other) : forces(other.forces) {
    }

    /** Copy operator */
    FExtendForces& operator=(const FExtendForces& other) {
        this->forces = other.forces;
        return *this;
    }

    /** Return the forces */
    const F3DPosition& getForces() const {
        return this->forces;
    }

    /** Set Forces */
    void incForces(const F3DPosition& inForces) {
        this->forces += inForces;
    }

    /** Set Forces with 3 FReals */
    void incForces(const FReal inFx, const FReal inFy, const FReal inFz) {
        this->forces.incX(inFx);
        this->forces.incY(inFy);
        this->forces.incZ(inFz);
    }

    /** set the forces from 3 variables */
    void setForces(const FReal inFx, const FReal inFy, const FReal inFz) {
        this->forces.setPosition(inFx , inFy, inFz);
    }
};


#endif //FEXTENDFORCES_HPP


