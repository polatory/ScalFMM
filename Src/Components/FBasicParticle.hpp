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
#ifndef FBASICPARTICLE_HPP
#define FBASICPARTICLE_HPP


#include "../Extensions/FExtendPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticle
* Please read the license
*
* This class defines a basic particle used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractParticle without
* inheriting from it.
*/
class FBasicParticle : public FExtendPosition{
public:
    /** Default destructor */
    virtual ~FBasicParticle(){
    }
};


#endif //FBASICPARTICLE_HPP


