// ===================================================================================
// Copyright ScalFmm 2011 INRIA
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FP2PPARTICLECONTAINER_HPP
#define FP2PPARTICLECONTAINER_HPP

#include "Components/FBasicParticleContainer.hpp"

template<int NRHS = 1, int NLHS = 1>
class FP2PParticleContainer : public FBasicParticleContainer<NRHS+4*NLHS> {
    typedef FBasicParticleContainer<NRHS+4*NLHS> Parent;

public:
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

#endif // FP2PPARTICLECONTAINER_HPP
