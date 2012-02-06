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
#ifndef FExtendPhysicalValue_HPP
#define FExtendPhysicalValue_HPP


#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPhysicalValue
* Please read the license
* This class is an extenssion.
* It proposes a physicalValue (FReal).
*/
class FExtendPhysicalValue {
protected:
    FReal physicalValue;   //< A simple physicalValue

public:
    /** Default constructor */
    FExtendPhysicalValue() : physicalValue(0) {
    }

    /** Copy constructor */
    FExtendPhysicalValue(const FExtendPhysicalValue& other) : physicalValue(other.physicalValue) {
    }

    /** Copy Constructor */
    FExtendPhysicalValue& operator=(const FExtendPhysicalValue& other) {
        this->physicalValue = other.physicalValue;
        return *this;
    }

    /** To get the physicalValue */
    FReal getPhysicalValue() const {
        return this->physicalValue;
    }

    /** To set the physicalValue */
    void setPhysicalValue(const FReal inphysicalValue) {
        this->physicalValue = inphysicalValue;
    }

};


#endif //FExtendPhysicalValue_HPP


