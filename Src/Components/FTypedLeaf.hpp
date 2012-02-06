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
#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP


#include "../Utils/FAssertable.hpp"
#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTypedLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particles
* (source XOR target) or simple system (source AND target)
*/
template< class ParticleClass , class ContainerClass>
class FTypedLeaf  : public FAbstractLeaf<ParticleClass,ContainerClass>, public FAssertable {
    ContainerClass sources;
    ContainerClass targets;

public:
    /** Default destructor */
    virtual ~FTypedLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        */
    void push(const ParticleClass& particle){
        if(particle.isTarget()) this->targets.push(particle);
        else if(particle.isSource()) this->sources.push(particle);
        else fassert(false, "Error particle has undefined type.", __LINE__, __FILE__);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    ContainerClass* getSrc() {
        return &this->sources;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    ContainerClass* getTargets() {
        return &this->targets;
    }

};


#endif //FTYPEDLEAF_HPP


