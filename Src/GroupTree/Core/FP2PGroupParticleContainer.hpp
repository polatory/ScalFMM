
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FP2PGROUPPARTICLECONTAINER_HPP
#define FP2PGROUPPARTICLECONTAINER_HPP

#include "FGroupAttachedLeaf.hpp"

template<int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FP2PGroupParticleContainer : public FGroupAttachedLeaf<NVALS*(NRHS+4*NLHS), FReal> {
    typedef FGroupAttachedLeaf<NVALS*(NRHS+4*NLHS), FReal> Parent;

public:
    FP2PGroupParticleContainer(){}
    FP2PGroupParticleContainer(const int inNbParticles, FReal* inPositionBuffer, const size_t inLeadingPosition,
                       FReal* inAttributesBuffer, const size_t inLeadingAttributes)
        : Parent(inNbParticles, inPositionBuffer, inLeadingPosition, inAttributesBuffer, inLeadingAttributes) {

    }


    FReal* getPhysicalValues(const int idxVals = 0, const int idxRhs = 0){
      return Parent::getAttribute((0+idxRhs)*NVALS+idxVals);
    }

    const FReal* getPhysicalValues(const int idxVals = 0, const int idxRhs = 0) const {
        return Parent::getAttribute((0+idxRhs)*NVALS+idxVals);
    }

    FReal* getPhysicalValuesArray(const int idxVals = 0, const int idxRhs = 0){
        return Parent::getRawData() + ((0+idxRhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    const FReal* getPhysicalValuesArray(const int idxVals = 0, const int idxRhs = 0) const {
        return Parent::getRawData() + ((0+idxRhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    int getLeadingDimension(){
        return Parent::getLeadingRawData();
    }

    FReal* getPotentials(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);
    }

    const FReal* getPotentials(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);
    }

    FReal* getPotentialsArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    const FReal* getPotentialsArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    FReal* getForcesX(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+NLHS+idxLhs)*NVALS+idxVals);
    }

    const FReal* getForcesX(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getAttribute((NRHS+NLHS+idxLhs)*NVALS+idxVals);
    }

    FReal* getForcesXArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    const FReal* getForcesXArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    FReal* getForcesY(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+2*NLHS+idxLhs)*NVALS+idxVals);
    }

    const FReal* getForcesY(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getAttribute((NRHS+2*NLHS+idxLhs)*NVALS+idxVals);
    }

    FReal* getForcesYArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+2*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    const FReal* getForcesYArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+2*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    FReal* getForcesZ(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+3*NLHS+idxLhs)*NVALS+idxVals);
    }

    const FReal* getForcesZ(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getAttribute((NRHS+3*NLHS+idxLhs)*NVALS+idxVals);
    }

    FReal* getForcesZArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+3*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    const FReal* getForcesZArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+3*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();
    }

    void resetForcesAndPotential(){
        for(int idx = 0 ; idx < 4*NLHS*NVALS ; ++idx){
            Parent::resetToInitialState(idx + NRHS*NVALS);
        }
    }

    int getNVALS() const {
        return NVALS;
    }

};

#endif // FP2PGROUPPARTICLECONTAINER_HPP
