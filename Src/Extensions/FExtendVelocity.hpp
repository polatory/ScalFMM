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
#ifndef FEXTENDVELOCITY_HPP
#define FEXTENDVELOCITY_HPP



#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendVelocity
* Please read the license
*
* This class is an extenssion.
* It proposes a 3d array as a velocity vector.
*/
class FExtendVelocity {
protected:
    F3DPosition velocity; //< 3D vector stored in a position object

public:
    /** Default constructor */
    FExtendVelocity() {
    }

    /** Copy constructor */
    FExtendVelocity(const FExtendVelocity& other) : velocity(other.velocity) {
    }

    /** Copy operator */
    FExtendVelocity& operator=(const FExtendVelocity& other) {
        this->velocity = other.velocity;
        return *this;
    }

    /** Return the velocity */
    const F3DPosition& getVelocity() const {
        return this->velocity;
    }

    /** Set Velocity */
    void incVelocity(const F3DPosition& inVelocity) {
        this->velocity += inVelocity;
    }

    /** Set Velocity with 3 FReals */
    void incVelocity(const FReal inVx, const FReal inVy, const FReal inVz) {
        this->velocity.incX(inVx);
        this->velocity.incY(inVy);
        this->velocity.incZ(inVz);
    }

    /** set the velocity from 3 variables */
    void setVelocity(const FReal inVx, const FReal inVy, const FReal inVz) {
        this->velocity.setPosition(inVx , inVy, inVz);
    }
};


#endif //FEXTENDVELOCITY_HPP


