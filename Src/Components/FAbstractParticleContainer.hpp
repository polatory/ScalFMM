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
#ifndef FABSTRACTPARTICLECONTAINER_HPP
#define FABSTRACTPARTICLECONTAINER_HPP

#include "../Utils/FGlobal.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FPoint.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @brief
* Please read the license
*
* This class define the method that every particle container
* has to implement.
*
* @warning Inherit from this class when implement a specific particle type
*/
template <class FReal>
class FAbstractParticleContainer {
public:
    /** Default destructor */
    virtual ~FAbstractParticleContainer(){
    }

    /**
     * This method should be inherited (or your leaf will do nothing)
     * the point is coming from the tree and is followed by what let the leaf
     * pass through its push method.
     */
    template<typename... Args>
    void push(const FPoint<FReal>& /*inParticlePosition*/, Args ... /*args*/){
        FLOG( FLog::Controller.write("Warning, push is not implemented!").write(FLog::Flush) );
    }
};


#endif //FABSTRACTPARTICLECONTAINER_HPP


