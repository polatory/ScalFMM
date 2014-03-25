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
#include "../../Components/FBasicCell.hpp"
#include "../../Extensions/FExtendCellType.hpp"

#include "../../Utils/FComplexe.hpp"

/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FUnifCell
 * Please read the license
 *
 * This class defines a cell used in the Lagrange based FMM.
 * 
 * PB: This class also contains the storage and accessors for the transformed 
 * expansion (in Fourier space, i.e. complex valued).
 *
 * @param NVALS is the number of right hand side.
 */
template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FUnifCell : public FBasicCell
{
  static const int VectorSize = TensorTraits<ORDER>::nnodes;
  static const int TransformedVectorSize = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);

  FReal multipole_exp[NRHS * NVALS * VectorSize]; //< Multipole expansion
  FReal     local_exp[NLHS * NVALS * VectorSize]; //< Local expansion
  // PB: Store multipole and local expansion in Fourier space
  FComplexe transformed_multipole_exp[NRHS * NVALS * TransformedVectorSize];
  FComplexe     transformed_local_exp[NLHS * NVALS * TransformedVectorSize];

public:
  FUnifCell(){
    memset(multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
    memset(local_exp, 0, sizeof(FReal) * NLHS * NVALS * VectorSize);
    memset(transformed_multipole_exp, 0, 
           sizeof(FComplexe) * NRHS * NVALS * TransformedVectorSize);
    memset(transformed_local_exp, 0, 
           sizeof(FComplexe) * NLHS * NVALS * TransformedVectorSize);
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

  /** Get Transformed Multipole */
  const FComplexe* getTransformedMultipole(const int inRhs) const{	
    return this->transformed_multipole_exp + inRhs*TransformedVectorSize;
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

  /** To get the leading dim of a vec */
  int getTransformedVectorSize() const{
    return TransformedVectorSize;
  }

  /** Make it like the begining */
  void resetToInitialState(){
    memset(multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
    memset(local_exp, 0, sizeof(FReal) * NLHS * NVALS * VectorSize);
    memset(transformed_multipole_exp, 0, 
           sizeof(FComplexe) * NRHS * NVALS * TransformedVectorSize);
    memset(transformed_local_exp, 0, 
           sizeof(FComplexe) * NLHS * NVALS * TransformedVectorSize);
  }

  ///////////////////////////////////////////////////////
  // to extend FAbstractSendable
  ///////////////////////////////////////////////////////
  template <class BufferWriterClass>
  void serializeUp(BufferWriterClass& buffer) const{
    buffer.write(multipole_exp, VectorSize*NVALS*NRHS);
    buffer.write(transformed_multipole_exp, TransformedVectorSize*NVALS*NRHS);
  }

  template <class BufferReaderClass>
  void deserializeUp(BufferReaderClass& buffer){
    buffer.fillArray(multipole_exp, VectorSize*NVALS*NRHS);
    buffer.fillArray(transformed_multipole_exp, TransformedVectorSize*NVALS*NRHS);
  }
  
  template <class BufferWriterClass>
  void serializeDown(BufferWriterClass& buffer) const{
    buffer.write(local_exp, VectorSize*NVALS*NLHS);
    buffer.write(transformed_local_exp, TransformedVectorSize*NVALS*NLHS);
  }
  
  template <class BufferReaderClass>
  void deserializeDown(BufferReaderClass& buffer){
    buffer.fillArray(local_exp, VectorSize*NVALS*NLHS);
    buffer.fillArray(transformed_local_exp, TransformedVectorSize*NVALS*NLHS);
  }
  
  ///////////////////////////////////////////////////////
  // to extend Serializable
  ///////////////////////////////////////////////////////
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FBasicCell::save(buffer);
    buffer.write(multipole_exp, VectorSize*NVALS*NRHS);
    buffer.write(transformed_multipole_exp, TransformedVectorSize*NVALS*NRHS);
    buffer.write(local_exp, VectorSize*NVALS*NLHS);
    buffer.write(transformed_local_exp, TransformedVectorSize*NVALS*NLHS);
  }
  
  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FBasicCell::restore(buffer);
    buffer.fillArray(multipole_exp, VectorSize*NVALS*NRHS);
    buffer.fillArray(transformed_multipole_exp, TransformedVectorSize*NVALS*NRHS);
    buffer.fillArray(local_exp, VectorSize*NVALS*NLHS);
    buffer.fillArray(transformed_local_exp, TransformedVectorSize*NVALS*NLHS);
  }
  
  static constexpr int GetSize(){
    return (NRHS+NLHS)*NVALS*VectorSize * (int) sizeof(FReal) + (NRHS+NLHS)*NVALS*TransformedVectorSize * (int) sizeof(FComplexe);
  }


};

template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FTypedUnifCell : public FUnifCell<ORDER,NRHS,NLHS,NVALS>, public FExtendCellType {
public:
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FUnifCell<ORDER,NRHS,NLHS,NVALS>::save(buffer);
    FExtendCellType::save(buffer);
  }
  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FUnifCell<ORDER,NRHS,NLHS,NVALS>::restore(buffer);
    FExtendCellType::restore(buffer);
  }
  void resetToInitialState(){
    FUnifCell<ORDER,NRHS,NLHS,NVALS>::resetToInitialState();
    FExtendCellType::resetToInitialState();
  }
};

#endif //FUNIFCELL_HPP
