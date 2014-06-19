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
#ifndef FP2PPARTICLECONTAINERINDEXED_HPP
#define FP2PPARTICLECONTAINERINDEXED_HPP

#include "../../Containers/FVector.hpp"

#include "FP2PParticleContainer.hpp"
#include "Components/FParticleType.hpp"

template<int NRHS = 1, int NLHS = 1>
class FP2PParticleContainerIndexed : public FP2PParticleContainer<NRHS,NLHS> {
    typedef FP2PParticleContainer<NRHS,NLHS> Parent;

    FVector<int> indexes;

public:
    template<typename... Args>
    void push(const FPoint& inParticlePosition, const int index, Args... args){
        Parent::push(inParticlePosition, args... );
        indexes.push(index);
    }

    template<typename... Args>
    void push(const FPoint& inParticlePosition, const FParticleType particleType, const int index, Args... args){
        Parent::push(inParticlePosition, particleType, args... );
        indexes.push(index);
    }

    const FVector<int>& getIndexes() const{
        return indexes;
    }
};

#endif // FP2PPARTICLECONTAINERINDEXED_HPP
