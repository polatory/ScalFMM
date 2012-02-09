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
#ifndef FABSTRACTLOADER_HPP
#define FABSTRACTLOADER_HPP


#include "../Utils/FGlobal.hpp"
class F3DPosition;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractLoader
* Please read the license
*
* This class defined the FMB usual loader. A loader is the component
* that fills an octree.
*
* If you want to use a specific file format you then need to inherite from this loader
* and implemente several methods.
*
* Please look at FBasicLoader or FFmaLoader to see an example.
*
* @warning Inherite from this class when defining a loader class
*/
template <class ParticleClass>
class FAbstractLoader {
public:	
	/** Default destructor */
	virtual ~FAbstractLoader(){
	}

        /**
        * Get the number of particles for this simulation
        * @return number of particles that the loader can fill
        */
        virtual FSize getNumberOfParticles() const = 0;

        /**
        * Get the center of the simulation box
        * @return box center needed by the octree
        */
        virtual F3DPosition getCenterOfBox() const = 0;

        /**
        * Get the simulation box width
        * @return box width needed by the octree
        */
        virtual FReal getBoxWidth() const = 0;

        /**
        * To know if the loader is valide (file opened, etc.)
        * @return true if file is open
        */
        virtual bool isOpen() const = 0;

        /**
        * Fill the next particle
        * @param inParticle the particle to fill
        */
        virtual void fillParticle(ParticleClass& inParticle) = 0;
};


#endif //FABSTRACTLOADER_HPP


