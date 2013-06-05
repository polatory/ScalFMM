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
#ifndef FTAYLORCELL_HPP
#define FTAYLORCELL_HPP

#include "../../Components/FBasicCell.hpp"
#include "../../Containers/FVector.hpp"

/**
 *@author Cyrille Piacibello
 *@class FTaylorCell
 *
 *This class is a cell used for the Taylor Expansion Kernel.
 *
 *
 */
template <int P>
class FTaylorCell : public FBasicCell {
protected:
  //Size of Multipole Vector
  static const int MultipoleSize = ((P+1)*(P+2)*(P+3))/6;
  //Size of Local Vector
  static const int LocalSize = ((P+1)*(P+2)*(P+3))/6;
  
  //Multipole vector
  FVector<FReal> multipole_exp = FVector(MultipoleSize);
  //Local vector
  FVector<FReal> local_exp =  FVector(LocalSize);

public:
  /**
   *Default Constructor
   */
  FTaylorCell(){
  }

  //Get multipole Vector
  FVector * getMultipole(void)
  {
    return multipole_exp;
  }

  //Get local Vector
  FVector * getLocal(void)
  {
    return local_exp;
  }

};

#endif
