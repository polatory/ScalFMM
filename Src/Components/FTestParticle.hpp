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
#ifndef FTESTPARTICLE_HPP
#define FTESTPARTICLE_HPP


#include "FBasicParticle.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTestParticle
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know whit it is.
*
* Particles just need the data down.
*/
class FTestParticle : public FBasicParticle {
protected:
    // To store data during downard pass
    long long int dataDown;
public:
    FTestParticle(): dataDown(0){
    }

    /** Default destructor */
    virtual ~FTestParticle(){
    }

    /** Get the down data */
    long long int getDataDown() const {
        return this->dataDown;
    }

    /** Set down data */
    void setDataDown(const long long int inData){
        this->dataDown = inData;
    }
};


#endif //FTESTPARTICLE_HPP


