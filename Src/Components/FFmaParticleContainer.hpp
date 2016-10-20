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
#ifndef FFMAPARTICLECONTAINER_HPP
#define FFMAPARTICLECONTAINER_HPP

#include "FBasicParticleContainer.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* Please read the license
*
* This class defines a particle container for FMA loader.
* position + 1 real for the physical value
*/
template <class FReal>
class FFmaParticleContainer : public FBasicParticleContainer<FReal, 1> {
    typedef FBasicParticleContainer<1> Parent;

public:
    /**
     * @brief getPhysicalValues to get the array of physical values
     * @return
     */
    FReal* getPhysicalValues(){
        return Parent::getAttribute<0>();
    }

    /**
     * @brief getPhysicalValues to get the array of physical values
     * @return
     */
    const FReal* getPhysicalValues() const {
        return Parent::getAttribute<0>();
    }
};


#endif //FFMAPARTICLECONTAINER_HPP


