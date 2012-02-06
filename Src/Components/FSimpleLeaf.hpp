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
#ifndef FSIMPLELEAF_HPP
#define FSIMPLELEAF_HPP


#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSimpleLeaf
* @brief
* Please read the license
* This class is used as a leaf in simple system (source AND target)
* here there is only one list that store all particles.
*/
template< class ParticleClass , class ContainerClass>
class FSimpleLeaf : public FAbstractLeaf<ParticleClass,ContainerClass> {
    ContainerClass particles;

public:
    /** Default destructor */
    virtual ~FSimpleLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        */
    void push(const ParticleClass& particle){
        this->particles.push(particle);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    ContainerClass* getSrc() {
        return &this->particles;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    ContainerClass* getTargets() {
        return &this->particles;
    }

};


#endif //FSIMPLELEAF_HPP


