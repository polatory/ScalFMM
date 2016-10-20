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

#ifndef FARRANGERPERIODIC_HPP
#define FARRANGERPERIODIC_HPP

#include "Utils/FGlobalPeriodic.hpp"
#include "Arranger/FOctreeArranger.hpp"

template <class FReal, class OctreeClass, class ContainerClass, class LeafInterface >
class FArrangerPeriodic : public FOctreeArranger<FReal,OctreeClass,ContainerClass,LeafInterface>{

    FReal getPeriodicPos(FReal pos, PeriodicCondition periodicPlus, PeriodicCondition periodicMinus,FReal maxDir,FReal minDir,const int dir){
        FReal res = pos;
        if( TestPeriodicCondition(dir, periodicPlus) ){
            while(res >= maxDir){
                res += (-(this->boxWidth));
            }
        }

        if( TestPeriodicCondition(dir, periodicMinus) ){
            while(res < minDir){
                res += (this->boxWidth);
            }
        }
        return res;
    }

public:

    FArrangerPeriodic(OctreeClass * octree) : FOctreeArranger<FReal,OctreeClass,ContainerClass,LeafInterface>(octree){
    }

    // To put in inhereed class
    void checkPosition(FPoint<FReal>& particlePos) override {
        particlePos.setX( getPeriodicPos( particlePos.getX(), DirMinusX, DirPlusX, (this->MaxBox).getX(),(this->MinBox).getX(),DirX));
        particlePos.setY( getPeriodicPos( particlePos.getY(), DirMinusY, DirPlusY, (this->MaxBox).getY(),(this->MinBox).getY(),DirY));
        particlePos.setZ( getPeriodicPos( particlePos.getZ(), DirMinusZ, DirPlusZ, (this->MaxBox).getZ(),(this->MinBox).getZ(),DirZ));
    }
};


#endif // FCONVERTERPERIODIC_HPP
