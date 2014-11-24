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
#ifndef FP2PLEAFINTERFACE_HPP
#define FP2PLEAFINTERFACE_HPP


#include "Arranger/FAbstractLeafInterface.hpp"
#include "Components/FBasicParticleContainer.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

template<class OctreeClass>
class FP2PLeafInterface : public FAbstractLeafInterface<OctreeClass,FP2PParticleContainerIndexed<> >{
private:
    typedef FP2PParticleContainerIndexed<> ContainerClass;
    ContainerClass * toStoreRemovedParts;
    FVector<int> indexes;

public:
    FP2PLeafInterface() : toStoreRemovedParts(),indexes(0){
        toStoreRemovedParts = new ContainerClass();
    }

    virtual ~FP2PLeafInterface(){
        indexes.clear();
        delete toStoreRemovedParts;
    }

    void getParticlePosition(ContainerClass* lf,const int idxPart, FPoint& newPos){
        newPos = FPoint(lf->getPositions()[0][idxPart],lf->getPositions()[1][idxPart],lf->getPositions()[2][idxPart]);
    }

    void removeFromLeafAndKeep(ContainerClass* lf, const int idxPart){
        //Store index
        indexes.push(lf->getIndexes()[idxPart]);
        //Store particles attributes
        toStoreRemovedParts->push(FPoint(lf->getPositions()[0][idxPart],lf->getPositions()[1][idxPart],lf->getPositions()[2][idxPart]),
                                  idxPart,
                                  lf->getPhysicalValues()[idxPart],
                                  lf->getForcesX()[idxPart],lf->getForcesY()[idxPart],lf->getForcesZ()[idxPart]);
        lf->removeParticles(&idxPart,1);
    }
    //Easy one use insert(std::array)
    void insertAllParticles(OctreeClass* tree){
        for(int idxToInsert = 0; idxToInsert<toStoreRemovedParts->getNbParticles() ; ++idxToInsert){
            tree->insert(FPoint(toStoreRemovedParts->getPositions()[0][idxToInsert],
                                toStoreRemovedParts->getPositions()[1][idxToInsert],
                                toStoreRemovedParts->getPositions()[2][idxToInsert]),
                         indexes[idxToInsert],
                         toStoreRemovedParts->getPhysicalValues()[idxToInsert],
                         toStoreRemovedParts->getForcesX()[idxToInsert],
                         toStoreRemovedParts->getForcesY()[idxToInsert],
                         toStoreRemovedParts->getForcesZ()[idxToInsert]);
        }
        toStoreRemovedParts->clear();
    }
};

#endif //FP2PLEAFINTERFACE_HPP
