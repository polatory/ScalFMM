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
#ifndef FEXTENDCOORDINATE_HPP
#define FEXTENDCOORDINATE_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCoordinate
* Please read the license
* This class is an extenssion.
* It proposes a mortonIndex.
*/
class FExtendCoordinate {
protected:
    FTreeCoordinate coordinate; //< The position

public:
    /** Default constructor */
    FExtendCoordinate() {
    }

    /** Copy constructor */
    FExtendCoordinate(const FExtendCoordinate& other) : coordinate(other.coordinate) {
    }

    /** Copy operator */
    FExtendCoordinate& operator=(const FExtendCoordinate& other) {
        this->coordinate = other.coordinate;
        return *this;
    }

    /** To get the position */
    const FTreeCoordinate& getCoordinate() const {
        return this->coordinate;
    }

    /** To set the position */
    void setCoordinate(const FTreeCoordinate& inCoordinate) {
        this->coordinate = inCoordinate;
    }

    /** To set the position from 3 FReals */
    void setCoordinate(const int inX, const int inY, const int inZ) {
        this->coordinate.setX(inX);
        this->coordinate.setY(inY);
        this->coordinate.setZ(inZ);
    }

};


#endif //FEXTENDCOORDINATE_HPP


