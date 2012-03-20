// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FHLOADER_HPP
#define FHLOADER_HPP


#include <cstdlib>
#include <time.h>


#include "../Utils/FGlobal.hpp"

#include "FAbstractLoader.hpp"
#include "../Utils/FPoint.hpp"

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
    const FPoint centerOfBox;    //< The center of box

public:
    /**
    * The constructor need the simulation data
    */
    FRandomLoader(const int inNbParticles, const FReal inBoxWidth = 1.0,
                  const FPoint& inCenterOfBox = FPoint(0,0,0), const unsigned int inSeed = static_cast<unsigned int>(time(NULL)))
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
    FPoint getCenterOfBox() const{
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
                  const FPoint& inCenterOfBox = FPoint(0,0,0), const unsigned int inSeed = static_cast<unsigned int>(time(NULL)))
        : FRandomLoader<ParticleClass>(inNbParticles,inBoxWidth,inCenterOfBox,inSeed) {
    }

    void fillParticle(ParticleClass& inParticle){
        FRandomLoader<ParticleClass>::fillParticle(inParticle);
        if(FRandomLoader<ParticleClass>::getRandom() > 0.5 ) inParticle.setAsTarget();
        else inParticle.setAsSource();
    }
};


#endif //FHLOADER_HPP


