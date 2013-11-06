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
#include "../../Utils/FMemUtils.hpp"
#include "../../Extensions/FExtendCellType.hpp"

/**
 *@author Cyrille Piacibello
 *@class FTaylorCell
 *
 *This class is a cell used for the Taylor Expansion Kernel.
 *
 *
 */
template <int P, int order>
class FTaylorCell : public FBasicCell {
protected:
  //Size of Multipole Vector
  static const int MultipoleSize = ((P+1)*(P+2)*(P+3))*order/6;
  //Size of Local Vector
  static const int LocalSize = ((P+1)*(P+2)*(P+3))*order/6;
  
  //Multipole vector
  FReal multipole_exp[MultipoleSize];
  //Local vector
  FReal local_exp[LocalSize];

public:
  /**
   *Default Constructor
   */
  FTaylorCell(){
    FMemUtils::memset(multipole_exp,0,MultipoleSize*sizeof(FReal(0)));
    FMemUtils::memset(local_exp,0,LocalSize*sizeof(FReal(0)));
  }

  //Get multipole Vector for setting
  FReal * getMultipole(void)
  {
    return multipole_exp;
  }

  //Get multipole Vector for reading
  const FReal * getMultipole(void) const
  {
    return multipole_exp;
  }

  //Get local Vector
  FReal * getLocal(void)
  {
    return local_exp;
  }

  //Get local Vector for reading
  const FReal * getLocal(void) const
  {
    return local_exp;
  }

  /** Make it like the begining */
  void resetToInitialState(){
      FMemUtils::memset(multipole_exp,0,MultipoleSize*sizeof(FReal(0)));
      FMemUtils::memset(local_exp,0,LocalSize*sizeof(FReal(0)));
  }

};

template <int P, int order>
class FTypedTaylorCell : public FTaylorCell<P,order>, public FExtendCellType {
public:
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FTaylorCell<P,order>::save(buffer);
        FExtendCellType::save(buffer);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FTaylorCell<P,order>::restore(buffer);
        FExtendCellType::restore(buffer);
    }
    void resetToInitialState(){
        FTaylorCell<P,order>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }
};

#endif
