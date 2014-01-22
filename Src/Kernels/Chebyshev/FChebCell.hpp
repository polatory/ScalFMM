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
#ifndef FCHEBCELL_HPP
#define FCHEBCELL_HPP

#include "../../Components/FBasicCell.hpp"
#include "../../Extensions/FExtendMortonIndex.hpp"
#include "../../Extensions/FExtendCoordinate.hpp"

#include "./FChebTensor.hpp"
#include "../../Extensions/FExtendCellType.hpp"

/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @class FChebCell
 * Please read the license
 *
 * This class defines a cell used in the Chebyshev based FMM.
 * @param NVALS is the number of right hand side.
 */
template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FChebCell : public FBasicCell
{
  static const int VectorSize = TensorTraits<ORDER>::nnodes * 2;

  FReal multipole_exp[NRHS * NVALS * VectorSize]; //< Multipole expansion
  FReal     local_exp[NLHS * NVALS * VectorSize]; //< Local expansion

public:
  FChebCell(){
    memset(multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
    memset(local_exp, 0, sizeof(FReal) * NLHS * NVALS * VectorSize);
  }

  ~FChebCell() {}
	
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
    memset(multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
    memset(local_exp, 0, sizeof(FReal) * NLHS * NVALS * VectorSize);
  }
  
  ///////////////////////////////////////////////////////
  // to extend FAbstractSendable
  ///////////////////////////////////////////////////////
  template <class BufferWriterClass>
  void serializeUp(BufferWriterClass& buffer) const{
    buffer.write(multipole_exp, VectorSize*NVALS*NRHS);
  }
  template <class BufferReaderClass>
  void deserializeUp(BufferReaderClass& buffer){
    buffer.fillArray(multipole_exp, VectorSize*NVALS*NRHS);
  }
  
  template <class BufferWriterClass>
  void serializeDown(BufferWriterClass& buffer) const{
    buffer.write(local_exp, VectorSize*NVALS*NLHS);
  }
  template <class BufferReaderClass>
  void deserializeDown(BufferReaderClass& buffer){
    buffer.fillArray(local_exp, VectorSize*NVALS*NLHS);
  }
  
  ///////////////////////////////////////////////////////
  // to extend Serializable
  ///////////////////////////////////////////////////////
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FBasicCell::save(buffer);
    buffer.write(multipole_exp, VectorSize*NVALS*NRHS);
    buffer.write(local_exp, VectorSize*NVALS*NLHS);
  }
  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FBasicCell::restore(buffer);
    buffer.fillArray(multipole_exp, VectorSize*NVALS*NRHS);
    buffer.fillArray(local_exp, VectorSize*NVALS*NLHS);
  }
  
  static int GetSize(){
    return sizeof(FReal) * VectorSize*(NRHS+NLHS)*NVALS;
  }

};

template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FTypedChebCell : public FChebCell<ORDER,NRHS,NLHS,NVALS>, public FExtendCellType {
public:
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FChebCell<ORDER,NRHS,NLHS,NVALS>::save(buffer);
    FExtendCellType::save(buffer);
  }
  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FChebCell<ORDER,NRHS,NLHS,NVALS>::restore(buffer);
    FExtendCellType::restore(buffer);
  }
  void resetToInitialState(){
    FChebCell<ORDER,NRHS,NLHS,NVALS>::resetToInitialState();
    FExtendCellType::resetToInitialState();
  }
};

#endif //FCHEBCELL_HPP


