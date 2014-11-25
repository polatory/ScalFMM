// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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
#ifndef FABSTRACTLEAFINTERFACE_HPP
#define FABSTRACTLEAFINTERFACE_HPP

template<class OctreeClass,class ParticleClass>
class FAbstractLeafInterface{
public:
    //    virtual FPoint getParticlePosition(ParticleClass* lf, const int idxPart) = 0;
    virtual void getParticlePosition(ParticleClass* lf, const int idxPart, FPoint& newPos) = 0;
    virtual void removeFromLeafAndKeep(ParticleClass* lf, const FPoint& newPos,const int idxPart) = 0;
    virtual void insertAllParticles(OctreeClass* tree) = 0;
};





#endif //FABSTRACTLEAFINTERFACE_HPP
