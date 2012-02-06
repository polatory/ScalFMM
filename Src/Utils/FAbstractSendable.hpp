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
#ifndef FABSTRACTSENDABLE_HPP
#define FABSTRACTSENDABLE_HPP



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractSendable
* Please read the license
*/
class FAbstractSendable {
protected:
    /** Empty Destructor */
    virtual ~FAbstractSendable(){}

    //static const int SerializedSizeUp = sizeof(?);
    virtual void serializeUp(void* const buffer) const  = 0;
    virtual void deserializeUp(const void* const buffer) = 0;

    //static const int SerializedSizeDown = sizeof(?);
    virtual void serializeDown(void* const buffer) const = 0;
    virtual void deserializeDown(const void* const buffer) = 0;
};

#endif //FABSTRACTSENDABLE_HPP


