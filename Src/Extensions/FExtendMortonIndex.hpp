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
#ifndef FEXTENDMORTONINDEX_HPP
#define FEXTENDMORTONINDEX_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendMortonIndex
* Please read the license
* This class is an extenssion.
* It proposes a mortonIndex.
*/
class FExtendMortonIndex {
protected:
    MortonIndex mortonIndex;    //< Morton index (need by most elements)

public:
    /** Default constructor */
    FExtendMortonIndex() : mortonIndex(0) {
    }

    /** Copy constructor */
    FExtendMortonIndex(const FExtendMortonIndex& other) : mortonIndex(other.mortonIndex) {
    }

    /** Copy operator */
    FExtendMortonIndex& operator=(const FExtendMortonIndex& other) {
        this->mortonIndex = other.mortonIndex;
        return *this;
    }

    /** To get the morton index */
    MortonIndex getMortonIndex() const {
        return this->mortonIndex;
    }

    /** To set the morton index */
    void setMortonIndex(const MortonIndex inMortonIndex) {
        this->mortonIndex = inMortonIndex;
    }

};


#endif //FEXTENDMORTONINDEX_HPP


