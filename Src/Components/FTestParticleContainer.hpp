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
#ifndef FTESTPARTICLECONTAINER_HPP
#define FTESTPARTICLECONTAINER_HPP


#include "FBasicParticleContainer.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTestParticle
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know whit it is.
* We store the positions + 1 long long int
*/
template <class FReal>
class FTestParticleContainer : public FBasicParticleContainer<FReal, 1, long long int> {
    typedef FBasicParticleContainer<FReal, 1, long long int> Parent;

public:
    /**
     * @brief getDataDown
     * @return
     */
    long long int* getDataDown(){
        return Parent::template getAttribute<0>();
    }

    /**
     * @brief getDataDown
     * @return
     */
    const long long int* getDataDown() const {
        return Parent::template getAttribute<0>();
    }
};


#endif //FTESTPARTICLECONTAINER_HPP


