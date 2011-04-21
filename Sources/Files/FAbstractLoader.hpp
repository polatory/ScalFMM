#ifndef FABSTRACTLOADER_HPP
#define FABSTRACTLOADER_HPP
// /!\ Please, you must read the license at the bottom of this page

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
        virtual long getNumberOfParticles() const = 0;

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
        virtual bool isValide() const = 0;

        /**
        * Fill the next particle
        * @param inParticle the particle to fill
        */
        virtual void fillParticle(ParticleClass* const inParticle) = 0;
};


#endif //FABSTRACTLOADER_HPP

// [--LICENSE--]
