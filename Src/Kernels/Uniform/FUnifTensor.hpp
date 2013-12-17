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
#ifndef FUNIFTENSOR_HPP
#define FUNIFTENSOR_HPP

#include "../../Utils/FMath.hpp"

#include "./FUnifRoots.hpp"
#include "./../Interpolation/FInterpTensor.hpp"


/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * Please read the license
 */



/**
 * @class FUnifTensor
 *
 * The class FUnifTensor provides function considering the tensor product
 * interpolation.
 *
 * @tparam ORDER interpolation order \f$\ell\f$
 */
template <int ORDER>
class FUnifTensor : public FInterpTensor<ORDER,FUnifRoots<ORDER>>
{
  enum {nnodes = TensorTraits<ORDER>::nnodes};
  typedef FUnifRoots<ORDER> BasisType;
  typedef FInterpTensor<ORDER,BasisType> ParentTensor;

 public:

  /**
   * Sets the diff of ids of the coordinates of all \f$\ell^6\f$ interpolation
   * nodes duet
   *
   * @param[out] NodeIdsDiff diff of ids of coordinates of interpolation nodes
   */
  static
    void setNodeIdsDiff(unsigned int NodeIdsDiff[nnodes*nnodes])
  {
    unsigned int node_ids[nnodes][3];
    ParentTensor::setNodeIds(node_ids);

    for (unsigned int i=0; i<nnodes; ++i) {
      for (unsigned int j=0; j<nnodes; ++j) {
        // 0 <= id < 2*ORDER-1
        unsigned int idl = node_ids[i][0]-node_ids[j][0]+ORDER-1;
        unsigned int idm = node_ids[i][1]-node_ids[j][1]+ORDER-1;
        unsigned int idn = node_ids[i][2]-node_ids[j][2]+ORDER-1;
        NodeIdsDiff[i*nnodes+j]
          = idn*(2*ORDER-1)*(2*ORDER-1) + idm*(2*ORDER-1) + idl;
      }
    }
  }

};


#endif /*FUNIFTENSOR_HPP*/
