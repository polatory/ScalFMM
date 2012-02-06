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
#ifndef FEXTENDPOSITION_HPP
#define FEXTENDPOSITION_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPosition
* Please read the license
* This class is an extenssion.
* It proposes a mortonIndex.
*/
class FExtendPosition {
protected:
    F3DPosition position; //< The position

public:
    /** Default constructor */
    FExtendPosition() {
    }

    /** Copy constructor */
    FExtendPosition(const FExtendPosition& other) : position(other.position) {
    }

    /** Copy operator */
    FExtendPosition& operator=(const FExtendPosition& other) {
        this->position = other.position;
        return *this;
    }

    /** To get the position */
    const F3DPosition& getPosition() const {
        return this->position;
    }

    /** To set the position */
    void setPosition(const F3DPosition& inPosition) {
        this->position = inPosition;
    }

    /** To set the position from 3 FReals */
    void setPosition(const FReal inX, const FReal inY, const FReal inZ) {
        this->position.setX(inX);
        this->position.setY(inY);
        this->position.setZ(inZ);
    }

    /** Set Position */
    void incPosition(const F3DPosition& inPosition) {
        this->position += inPosition;
    }

    /** Set Position with 3 FReals */
    void incPosition(const FReal inPx, const FReal inPy, const FReal inPz) {
        this->position.incX(inPx);
        this->position.incY(inPy);
        this->position.incZ(inPz);
    }

};


#endif //FEXTENDPOSITION_HPP


