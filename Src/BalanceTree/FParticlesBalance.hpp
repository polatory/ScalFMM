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

#ifndef FPARTICLESBALANCE_H
#define FPARTICLESBALANCE_H

#include "./FAbstractBalanceAlgorithm.hpp"
#include "../Utils/FMath.hpp"

/**
 * @author Cyrille Piacibello
 * @class FLeafBalance
 *
 * @brief This class inherits from FAbstractBalanceAlgorithm. It
 * provides balancing methods based on particles distribution.
 */
class FParticlesBalance : public FAbstractBalanceAlgorithm{

public:
  /**
   * getRight interval based on particles distribution
   */
  FSize getRight(const FSize numberOfLeaves,
                 const int numberOfProc, const int idxOfProc){
    int acc = 0;
    FSize i = 0;
    const double step = (double(numberOfPart) / double(numberOfProc));
    FSize aimRight = FSize(FMath::Ceil(step * double(idxOfProc+1)));
    if(aimRight > numberOfPart) aimRight = numberOfPart;
    while(acc < aimRight){
      acc+=numberOfPartPerLeaf[i];
      ++i;
    }
    if(FMath::Abs(aimRight-acc) < FMath::Abs(aimRight-(acc-numberOfPartPerLeaf[i]))) return i;
    else
      return i-1;
  }

  /**
   * get left interval based on particles distribution
   */
  FSize getLeft(const FSize numberOfLeaves,
                const int numberOfProc, const int idxOfProc){
    int acc = 0;
    FSize i = 0;
    const double step = (double(numberOfPart) / double(numberOfProc));
    const FSize aimLeft = FSize(FMath::Ceil(step * double(idxOfProc)));
    while (acc < aimLeft){
      acc+=numberOfPartPerLeaf[i];
      ++i;
    }
    if(FMath::Abs(aimLeft-acc) < FMath::Abs(aimLeft-(acc-numberOfPartPerLeaf[i]))) return i;
    else
      return i-1;
  }

};


#endif // FPARTICLESBALANCE_H
