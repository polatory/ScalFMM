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


