
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPTESTPARTICLECONTAINER_HPP
#define FGROUPTESTPARTICLECONTAINER_HPP

#include "FGroupAttachedLeaf.hpp"

class FGroupTestParticleContainer : public FGroupAttachedLeaf<2, long long int> {
    typedef FGroupAttachedLeaf<2, long long int> Parent;

public:
    FGroupTestParticleContainer(){}
    FGroupTestParticleContainer(const int inNbParticles, FReal* inPositionBuffer, const size_t inLeadingPosition,
                                long long int* inAttributesBuffer, const size_t inLeadingAttributes)
        : Parent(inNbParticles, inPositionBuffer, inLeadingPosition, inAttributesBuffer, inLeadingAttributes) {

    }

    /**
     * @brief getDataDown
     * @return
     */
    long long int* getDataDown(){
        return Parent::getAttribute<0>();
    }

    /**
     * @brief getDataDown
     * @return
     */
    const long long int* getDataDown() const {
        return Parent::getAttribute<0>();
    }
};

#endif // FGROUPTESTPARTICLECONTAINER_HPP
