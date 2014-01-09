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
#ifndef FUNIFCELL_HPP
#define FUNIFCELL_HPP


#include "../../Extensions/FExtendMortonIndex.hpp"
#include "../../Extensions/FExtendCoordinate.hpp"

#include "./FUnifTensor.hpp"

#include "../../Extensions/FExtendCellType.hpp"

#include "../../Utils/FComplexe.hpp"

/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FUnifCell
 * Please read the license
 *
 * This class defines a cell used in the Lagrange based FMM.
 * 
 * PB: !!! exactly the same as FChebCell except the TensorTraits used is 
 * the one from FUnifTensor. This is a quick work around to avoid multiple 
 * definition of TensorTraits until it is factorized for all interpolations.
 *
 * PB: This class also contains the storage and accessors for the tranformed 
 * expansion  (in Fourier space, i.e. complex valued).
 *
 * @param NVALS is the number of right hand side.
 */
template <int ORDER, int NMUL = 1, int NLOC = 1>
class FUnifCell : public FExtendMortonIndex, public FExtendCoordinate
{
  static const int VectorSize = TensorTraits<ORDER>::nnodes;
  static const int TransformedVectorSize = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);

  FReal multipole_exp[NMUL * VectorSize]; //< Multipole expansion
  FReal     local_exp[NLOC * VectorSize]; //< Local expansion

  // PB: Store multipole and local expansion in Fourier space
  FComplexe transformed_multipole_exp[NMUL * TransformedVectorSize];
  FComplexe     transformed_local_exp[NLOC * TransformedVectorSize];

public:
  FUnifCell(){
    memset(multipole_exp, 0, sizeof(FReal) * NMUL * VectorSize);
    memset(local_exp, 0, sizeof(FReal) * NLOC * VectorSize);
    memset(transformed_multipole_exp, 0, 
           sizeof(FComplexe) * NMUL * TransformedVectorSize);
    memset(transformed_local_exp, 0, 
           sizeof(FComplexe) * NLOC * TransformedVectorSize);
  }

  ~FUnifCell() {}

  /** Get Multipole */
  const FReal* getMultipole(const int inRhs) const
  {	return this->multipole_exp + inRhs*VectorSize;
  }
  /** Get Local */
  const FReal* getLocal(const int inRhs) const{
    return this->local_exp + inRhs*VectorSize;
  }

  /** Get Multipole */
  FReal* getMultipole(const int inRhs){
    return this->multipole_exp + inRhs*VectorSize;
  }
  /** Get Local */
  FReal* getLocal(const int inRhs){
    return this->local_exp + inRhs*VectorSize;
  }

  /** To get the leading dim of a vec */
  int getVectorSize() const{
    return VectorSize;
  }

  /** Make it like the begining */
  void resetToInitialState(){
    memset(multipole_exp, 0, sizeof(FReal) * NMUL * VectorSize);
    memset(local_exp, 0, sizeof(FReal) * NLOC * VectorSize);
    memset(transformed_multipole_exp, 0, 
           sizeof(FComplexe) * NMUL * TransformedVectorSize);
    memset(transformed_local_exp, 0, 
           sizeof(FComplexe) * NLOC * TransformedVectorSize);
  }

  /** Get Transformed Multipole */
  const FComplexe* getTransformedMultipole(const int inRhs) const
  {	return this->transformed_multipole_exp + inRhs*TransformedVectorSize;
  }
  /** Get Transformed Local */
  const FComplexe* getTransformedLocal(const int inRhs) const{
    return this->transformed_local_exp + inRhs*TransformedVectorSize;
  }

  /** Get Transformed Multipole */
  FComplexe* getTransformedMultipole(const int inRhs){
    return this->transformed_multipole_exp + inRhs*TransformedVectorSize;
  }
  /** Get Transformed Local */
  FComplexe* getTransformedLocal(const int inRhs){
    return this->transformed_local_exp + inRhs*TransformedVectorSize;
  }

};

template <int ORDER, int NMUL = 1, int NLOC = 1>
class FTypedUnifCell : public FUnifCell<ORDER,NMUL,NLOC>, public FExtendCellType {
public:
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FUnifCell<ORDER,NMUL,NLOC>::save(buffer);
    FExtendCellType::save(buffer);
  }
  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FUnifCell<ORDER,NMUL,NLOC>::restore(buffer);
    FExtendCellType::restore(buffer);
  }
  void resetToInitialState(){
    FUnifCell<ORDER,NMUL,NLOC>::resetToInitialState();
    FExtendCellType::resetToInitialState();
  }
};

#endif //FUNIFCELL_HPP
