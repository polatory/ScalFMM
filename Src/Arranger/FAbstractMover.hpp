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
#ifndef FABSTRACTLEAFINTERFACE_HPP
#define FABSTRACTLEAFINTERFACE_HPP

template<class FReal,class OctreeClass,class ParticleClass>
class FAbstractMover{
public:
    virtual void getParticlePosition(ParticleClass* lf, const FSize idxPart, FPoint<FReal>* particlePos) = 0;
    virtual void removeFromLeafAndKeep(ParticleClass* lf, const FPoint<FReal>& particlePos, const FSize idxPart, FParticleType type) = 0;
    virtual void insertAllParticles(OctreeClass* tree) = 0;
};





#endif //FABSTRACTLEAFINTERFACE_HPP
