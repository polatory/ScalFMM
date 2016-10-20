// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#ifndef FROTATIONLPARTICLECONTAINER_HPP
#define FROTATIONLPARTICLECONTAINER_HPP


#include "../../Components/FBasicParticleContainer.hpp"

template <class FReal>
class FRotationParticleContainer : public FBasicParticleContainer<FReal, 5> {
    typedef FBasicParticleContainer<FReal, 5> Parent;
public:
    FReal* getPhysicalValues(){
        return Parent::getAttribute(0);
    }

    const FReal* getPhysicalValues() const {
        return Parent::getAttribute(0);
    }

    FReal* getPotentials(){
        return Parent::getAttribute(1);
    }

    const FReal* getPotentials() const {
        return Parent::getAttribute(1);
    }

    FReal* getForcesX(){
        return Parent::getAttribute(2);
    }

    const FReal* getForcesX() const {
        return Parent::getAttribute(2);
    }

    FReal* getForcesY(){
        return Parent::getAttribute(3);
    }

    const FReal* getForcesY() const {
        return Parent::getAttribute(3);
    }

    FReal* getForcesZ(){
        return Parent::getAttribute(4);
    }

    const FReal* getForcesZ() const {
        return Parent::getAttribute(4);
    }
};

#endif // FROTATIONLPARTICLECONTAINER_HPP
