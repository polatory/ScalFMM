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
#ifndef FHLOADER_HPP
#define FHLOADER_HPP


#include <cstdlib>
#include <time.h>


#include "../Utils/FGlobal.hpp"

#include "FAbstractLoader.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FRandomLoader
* Please read the license
*/
template <class ParticleClass>
class FRandomLoader : public FAbstractLoader<ParticleClass> {
protected:
    const int nbParticles;            //< the number of particles
    const FReal boxWidth;             //< the box width
    const F3DPosition centerOfBox;    //< The center of box

public:
    /**
    * The constructor need the simulation data
    */
    FRandomLoader(const int inNbParticles, const FReal inBoxWidth = 1.0,
                  const F3DPosition& inCenterOfBox = F3DPosition(0,0,0), const unsigned int inSeed = static_cast<unsigned int>(time(NULL)))
        : nbParticles(inNbParticles), boxWidth(inBoxWidth), centerOfBox(inCenterOfBox) {
        srand(inSeed);
    }

    /**
    * Default destructor
    */
    virtual ~FRandomLoader(){
    }

    /**
      * @return true
      */
    bool isOpen() const{
        return true;
    }

    /**
      * To get the number of particles from this loader
      * @param the number of particles the loader can fill
      */
    FSize getNumberOfParticles() const{
        return FSize(this->nbParticles);
    }

    /**
      * The center of the box
      * @return box center
      */
    F3DPosition getCenterOfBox() const{
        return this->centerOfBox;
    }

    /**
      * The box width
      * @return box width
      */
    FReal getBoxWidth() const{
        return this->boxWidth;
    }

    /**
      * Fill a particle
      * @warning to work with the loader, particles has to expose a setPosition method
      * @param the particle to fill
      */
    void fillParticle(ParticleClass& inParticle){
        inParticle.setPosition(
                    (getRandom() * boxWidth) + centerOfBox.getX() - boxWidth/2,
                    (getRandom() * boxWidth) + centerOfBox.getY() - boxWidth/2,
                    (getRandom() * boxWidth) + centerOfBox.getZ() - boxWidth/2);
    }

    /** Get a random number between 0 & 1 */
    FReal getRandom() const{
        return FReal(rand())/FReal(RAND_MAX);
    }
};



template <class ParticleClass>
class FRandomLoaderTsm : public FRandomLoader<ParticleClass> {
public:
    FRandomLoaderTsm(const int inNbParticles, const FReal inBoxWidth = 1.0,
                  const F3DPosition& inCenterOfBox = F3DPosition(0,0,0), const unsigned int inSeed = static_cast<unsigned int>(time(NULL)))
        : FRandomLoader<ParticleClass>(inNbParticles,inBoxWidth,inCenterOfBox,inSeed) {
    }

    void fillParticle(ParticleClass& inParticle){
        FRandomLoader<ParticleClass>::fillParticle(inParticle);
        if(FRandomLoader<ParticleClass>::getRandom() > 0.5 ) inParticle.setAsTarget();
        else inParticle.setAsSource();
    }
};


#endif //FHLOADER_HPP


