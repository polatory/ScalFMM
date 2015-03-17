
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FP2PGROUPPARTICLECONTAINER_HPP
#define FP2PGROUPPARTICLECONTAINER_HPP

#include "FGroupAttachedLeaf.hpp"

template<int NRHS = 1, int NLHS = 1>
class FP2PGroupParticleContainer : public FGroupAttachedLeaf<NRHS+4*NLHS, FReal> {
    typedef FGroupAttachedLeaf<NRHS+4*NLHS, FReal> Parent;

public:
    FP2PGroupParticleContainer(){}
    FP2PGroupParticleContainer(const int inNbParticles, FReal* inPositionBuffer, const size_t inLeadingPosition,
                       FReal* inAttributesBuffer, const size_t inLeadingAttributes)
        : Parent(inNbParticles, inPositionBuffer, inLeadingPosition, inAttributesBuffer, inLeadingAttributes) {

    }

    FReal* getPhysicalValues(const int idxRhs = 0){
        return Parent::getAttribute(0+idxRhs);
    }

    const FReal* getPhysicalValues(const int idxRhs = 0) const {
        return Parent::getAttribute(0+idxRhs);
    }

    FReal* getPotentials(const int idxLhs = 0){
        return Parent::getAttribute(NRHS+idxLhs);
    }

    const FReal* getPotentials(const int idxLhs = 0) const {
        return Parent::getAttribute(NRHS+idxLhs);
    }

    FReal* getForcesX(const int idxLhs = 0){
        return Parent::getAttribute(NRHS+NLHS+idxLhs);
    }

    const FReal* getForcesX(const int idxLhs = 0) const {
        return Parent::getAttribute(NRHS+NLHS+idxLhs);
    }

    FReal* getForcesY(const int idxLhs = 0){
        return Parent::getAttribute(NRHS+2*NLHS+idxLhs);
    }

    const FReal* getForcesY(const int idxLhs = 0) const {
        return Parent::getAttribute(NRHS+2*NLHS+idxLhs);
    }

    FReal* getForcesZ(const int idxLhs = 0){
        return Parent::getAttribute(NRHS+3*NLHS+idxLhs);
    }

    const FReal* getForcesZ(const int idxLhs = 0) const {
        return Parent::getAttribute(NRHS+3*NLHS+idxLhs);
    }
};

#endif // FP2PGROUPPARTICLECONTAINER_HPP
