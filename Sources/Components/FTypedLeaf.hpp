#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Containers/FList.hpp"
#include "../Utils/FAssertable.hpp"
#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTypedLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particules
* (source XOR target) or simple system (source AND target)
*/
template< class ParticuleClass >
class FTypedLeaf  : public FAbstractLeaf<ParticuleClass>, public FAssertable {
    FList<ParticuleClass*> sources;
    FList<ParticuleClass*> targets;

public:
    /** Default destructor */
    virtual ~FTypedLeaf(){
    }

    /**
        * To add a new particule in the leaf
        * @param particule the new particule
        */
    void push(ParticuleClass* const particule){
        if(particule->isTarget()) this->targets.pushFront(particule);
        else if(particule->isSource()) this->sources.pushFront(particule);
        else assert(false, "Error particule has undefined type.", __LINE__, __FILE__);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particules that are sources
        */
    FList<ParticuleClass*>* getSources() {
        return &this->sources;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particules that are targets
        */
    FList<ParticuleClass*>* getTargets() {
        return &this->targets;
    }

};


#endif //FTYPEDLEAF_HPP

// [--LICENSE--]
