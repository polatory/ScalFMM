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
#ifndef FABSTRACTPARTICLE_HPP
#define FABSTRACTPARTICLE_HPP


/* forward declaration to avoid include */
class F3DPosition;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractBody
* @brief
* Please read the license
*
* This class define the method that every particle class
* has to implement.
*
* In fact FOctree & FFmmAlgorithm need this function to be implemented.
* But you cannot use this interface with the extension (as an example :
* because the compiler will faill to know if getPosition is coming
* from this interface or from the extension)
*
*
* @warning Inherite from this class when implement a specific particle type
*/
class FAbstractParticle{
public:	
	/** Default destructor */
	virtual ~FAbstractParticle(){
	}

	/**
	* Must be implemented by each user Particle class
	* @return the position of the current cell
	*/
        virtual const F3DPosition& getPosition() const = 0;
};


#endif //FABSTRACTPARTICLE_HPP


