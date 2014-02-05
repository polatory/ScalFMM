#ifndef FP2PPARTICLECONTAINERINDEXED_HPP
#define FP2PPARTICLECONTAINERINDEXED_HPP

#include "../../Containers/FVector.hpp"

#include "FP2PParticleContainer.hpp"
#include "../../Components/FParticleType.hpp"

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
