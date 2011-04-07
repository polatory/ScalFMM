#ifndef FSIMPLELEAF_HPP
#define FSIMPLELEAF_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Containers/FList.hpp"
#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSimpleLeaf
* @brief
* Please read the license
* This class is used as a leaf in simple system (source AND target)
* here there is only one list that store all particules.
*/
template< class ParticuleClass >
class FSimpleLeaf : public FAbstractLeaf<ParticuleClass> {
    FList<ParticuleClass*> particules;

public:
    /** Default destructor */
    virtual ~FSimpleLeaf(){
    }

    /**
        * To add a new particule in the leaf
        * @param particule the new particule
        */
    void push(ParticuleClass* const particule){
        this->particules.pushFront(particule);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particules that are sources
        */
    FList<ParticuleClass*>* getSources() {
        return &this->particules;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particules that are targets
        */
    FList<ParticuleClass*>* getTargets() {
        return &this->particules;
    }

};


#endif //FSIMPLELEAF_HPP

// [--LICENSE--]
