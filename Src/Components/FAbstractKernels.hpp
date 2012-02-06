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
#ifndef FABSTRACTKERNELS_HPP
#define FABSTRACTKERNELS_HPP


#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractKernels
* @brief
* Please read the license
*
* If you want to create you own kernels you have to inherit from this class.
*/
template< class ParticleClass, class CellClass, class ContainerClass >
class FAbstractKernels{
public:
    /** Default destructor */
    virtual ~FAbstractKernels(){
    }

    /**
        * P2M
        * particles to multipole
        * @param pole the multipole to fill using the particles
        * @param particles the particles from the same spacial boxe
        */
    virtual void P2M(CellClass* const pole, const ContainerClass* const particles) = 0;

    /**
        * M2M
        * Multipole to multipole
        * @param pole the father (the boxe that contains other ones)
        * @param child the boxe to take values from
        * @param the current computation level
        * the child array has a size of 8 elements (address if exists or 0 otherwise).
        * You must test if a pointer is 0 to know if an element exists inside this array
        */
    virtual void M2M(CellClass* const FRestrict pole, const CellClass*const FRestrict *const FRestrict child, const int inLevel) = 0;

    /**
        * M2L
        * Multipole to local
        * @param local the element to fill using distant neighbors
        * @param distantNeighbors is an array containing fathers's direct neighbors's child - direct neigbors
        * @param size the number of neighbors
        * @param inLevel the current level of the computation
        */
    virtual void M2L(CellClass* const FRestrict local, const CellClass* distantNeighbors[189],
                     const int size, const int inLevel) = 0;

    /**
        * L2L
        * Local to local
        * @param the father to take value from
        * @param the child to downward values (child may have already been impacted by M2L)
        * @param level the current level of computation
        * the child array has a size of 8 elements (address if exists or 0 otherwise).
        * You must test if a pointer is 0 to know if an element exists inside this array
        */
    virtual void L2L(const CellClass* const FRestrict local, CellClass* FRestrict * const FRestrict child, const int inLevel) = 0;

    /**
        * L2P
        * Local to particles
        * @param local the leaf element (smaller boxe local element)
        * @param particles the list of particles inside this boxe
        */
    virtual void L2P(const CellClass* const local, ContainerClass* const particles) = 0;

    /**
        * P2P
        * Particles to particles
        * @param targets current boxe targets particles
        * @param sources current boxe sources particles
        * @param directNeighborsParticles the particles from direct neighbors (this is an array of list)
        * @param size the number of direct neighbors (the size of the array directNeighborsParticles)
        */
    virtual void P2P(ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                     const ContainerClass* const directNeighborsParticles[26], const int size) = 0;

    /**
        * P2P
        * Particles to particles
        * @param inCurrentLeafIndex
        * @param targets current boxe targets particles
        * @param sources current boxe sources particles
        * @param directNeighborsParticles the particles from direct neighbors (this is an array of list)
        * @param inNeighborsIndex the indexes of neighbors
        * @param size the number of direct neighbors (the size of the array directNeighborsParticles)
        */
    virtual void P2P(const MortonIndex inCurrentLeafIndex,
             ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
             ContainerClass* const directNeighborsParticles[26], const MortonIndex inNeighborsIndex[26], const int size) = 0;


    //////////////////////////////////////////////////////////////////////////////
    // Periodic methods
    //////////////////////////////////////////////////////////////////////////////

    /**
        * M2L
        * Multipole to local
        * @param local the element to fill using distant neighbors
        * @param distantNeighbors is an array containing fathers's direct neighbors's child - direct neigbors
        * @param neighborsRelativePositions the relative position of the neighbors (can be -2,-2,-2)
        * @param size the number of neighbors
        * @param inLevel the current level of the computation
        */
    virtual void M2L(CellClass* const FRestrict local, const CellClass* distantNeighbors[189],
                     const FTreeCoordinate neighborsRelativePositions[189], const int size, const int level) = 0;


    /**
        * P2P
        * Particles to particles
        * @param inCurrentLeafIndex the leaf index
        * @param targets current boxe targets particles
        * @param sources current boxe sources particles
        * @param directNeighborsParticles the particles from direct neighbors (this is an array of list)
        * @param neighborsRelativePositions the relative position of neighbors
        * @param size the number of direct neighbors (the size of the array directNeighborsParticles)
        */
    virtual void P2P(const MortonIndex inCurrentLeafIndex,
                    ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                    ContainerClass* const directNeighborsParticles[26], const FTreeCoordinate neighborsRelativeOffset[26],
                    const int size) = 0;
};


#endif //FABSTRACTKERNELS_HPP


