// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
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
    const size_t nbParticles;            //< the number of particles
    const FReal boxWidth;             //< the box width
    const FPoint centerOfBox;    //< The center of box

public:
    /**
    * The constructor need the simulation data
    */
    FRandomLoader(const size_t inNbParticles, const FReal inBoxWidth = 1.0,
                  const FPoint& inCenterOfBox = FPoint(0,0,0), const unsigned int inSeed = static_cast<unsigned int>(0))
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


/** This class is a random loader but it also generate
  * randomly the particles type (target or source)
  */
template <class ParticleClass>
class FRandomLoaderTsm : public FRandomLoader<ParticleClass> {
public:
    FRandomLoaderTsm(const size_t inNbParticles, const FReal inBoxWidth = 1.0,
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


