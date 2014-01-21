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
#ifndef FINTERPMATRIXKERNEL_HPP
#define FINTERPMATRIXKERNEL_HPP

#include <stdexcept>

#include "../../Utils/FPoint.hpp"
#include "../../Utils/FNoCopyable.hpp"
#include "../../Utils/FMath.hpp"
#include "../../Utils/FBlas.hpp"

// extendable
enum KERNEL_FUNCTION_IDENTIFIER {ONE_OVER_R,
                                 ONE_OVER_R_SQUARED,
                                 LENNARD_JONES_POTENTIAL,
                                 RX,
                                 ID_OVER_R,
                                 R_IJ,
                                 R_IJK};

// probably not extendable :)
enum KERNEL_FUNCTION_TYPE {HOMOGENEOUS, NON_HOMOGENEOUS};


/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FInterpMatrixKernels
 * Please read the license
 *
 * This class provides the evaluators and scaling functions of the matrix 
 * kernels. A matrix kernel should be understood in the sense of a kernel 
 * of interaction (or the fundamental solution of a given equation). 
 * It can either be scalar (DIM=1) or tensorial (DIM>1) depending on the 
 * dimension of the equation considered. DIM denotes the number of components 
 * that are actually stored (e.g. 6 for a \f$3\times3\f$ symmetric tensor). 
 * Notes on application scheme: 
 * Let there be a kernel \f$K\f$ such that \f$X_i=K_{ij}Y_j\f$
 * with \f$X\f$ the lhs of size NLHS and \f$Y\f$ the rhs of size NRHS. 
 * The table applyTab provides the indices in the reduced storage table 
 * corresponding to the application scheme depicted ealier.
 *
 */
struct FInterpAbstractMatrixKernel : FNoCopyable
{
  virtual ~FInterpAbstractMatrixKernel(){} // to remove warning
  virtual FReal evaluate(const FPoint&, const FPoint&) const = 0;
  // I need both functions because required arguments are not always given
  virtual FReal getScaleFactor(const FReal, const int) const = 0;
  virtual FReal getScaleFactor(const FReal) const = 0;
};



/// One over r
struct FInterpMatrixKernelR : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = ONE_OVER_R;
  static const  unsigned int DIM = 1; //PB: dimension of kernel
  static const unsigned int NRHS = 1;
  static const unsigned int NLHS = 1;

  FInterpMatrixKernelR(const unsigned int = 0) {}

  // returns position in reduced storage
  int getPosition(const unsigned int) const
  {return 0;} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    return FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
                                   xy.getY()*xy.getY() +
                                   xy.getZ()*xy.getZ());
  }

  FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
  {
    const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
    return getScaleFactor(CellWidth);
  }

  FReal getScaleFactor(const FReal CellWidth) const
  {
    return FReal(2.) / CellWidth;
  }
};



/// One over r^2
struct FInterpMatrixKernelRR : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = ONE_OVER_R_SQUARED;
  static const  unsigned int DIM = 1; //PB: dimension of kernel
  static const unsigned int NRHS = 1;
  static const unsigned int NLHS = 1;

  FInterpMatrixKernelRR(const unsigned int) {}

  // returns position in reduced storage
  int getPosition(const unsigned int) const
  {return 0;} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    return FReal(1.) / FReal(xy.getX()*xy.getX() +
                             xy.getY()*xy.getY() +
                             xy.getZ()*xy.getZ());
  }

  FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
  {
    const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
    return getScaleFactor(CellWidth);
  }

  FReal getScaleFactor(const FReal CellWidth) const
  {
    return FReal(4.) / CellWidth;
  }
};



/// One over r^12 - One over r^6
struct FInterpMatrixKernelLJ : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = LENNARD_JONES_POTENTIAL;
  static const  unsigned int DIM = 1; //PB: dimension of kernel
  static const unsigned int NRHS = 1;
  static const unsigned int NLHS = 1;

  FInterpMatrixKernelLJ(const unsigned int = 0) {}

  // returns position in reduced storage
  int getPosition(const unsigned int) const
  {return 0;} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    const FReal r = xy.norm();
    const FReal r3 = r*r*r;
    const FReal one_over_r6 = FReal(1.) / (r3*r3);
    //return one_over_r6 * one_over_r6;
    //return one_over_r6;
    return one_over_r6 * one_over_r6 - one_over_r6;
  }

  FReal getScaleFactor(const FReal, const int) const
  {
    // return 1 because non homogeneous kernel functions cannot be scaled!!!
    return FReal(1.);
  }

  FReal getScaleFactor(const FReal) const
  {
    // return 1 because non homogeneous kernel functions cannot be scaled!!!
    return FReal(1.);
  }

};


/// r_0 / R^2
struct FInterpMatrixKernelRX : FInterpAbstractMatrixKernel
{
  // PB: leave NON_HOMOGENEOUS while it is used as test matrix kernel
  static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = RX;
  static const  unsigned int DIM = 1; //PB: dimension of kernel
  static const unsigned int NRHS = 1;
  static const unsigned int NLHS = 1;

  FInterpMatrixKernelRX(const unsigned int = 0) {}

  // returns position in reduced storage
  int getPosition(const unsigned int) const
  {return 0;} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    const FReal r = xy.norm();
    const FReal r2 = r*r;
    return xy.getX()/**xy.getX()*//r2;
  }

  FReal getScaleFactor(const FReal, const int) const
  {
    // return 1 because non homogeneous kernel functions cannot be scaled!!!
    return FReal(1.);
  }

  FReal getScaleFactor(const FReal) const
  {
    // return 1 because non homogeneous kernel functions cannot be scaled!!!
    return FReal(1.);
  }

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// Tensorial Matrix Kernels (DIM>1)
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// Test Tensorial kernel 1/R*Id_3
struct FInterpMatrixKernel_IOR : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = ID_OVER_R;
  static const unsigned int DIM = 6; //PB: dimension of kernel
  const unsigned int indexTab[12]={0,0,0,1,1,2,
                                   0,1,2,1,2,2};
  static const unsigned int NRHS = 3;
  static const unsigned int NLHS = 3;

  // store positions in sym tensor 
  const unsigned int applyTab[9]={0,1,2,
                                  1,3,4,
                                  2,4,5};

  const unsigned int _i,_j;

  FInterpMatrixKernel_IOR(const unsigned int d = 0)
  : _i(indexTab[d]), _j(indexTab[d+DIM])
  {}

  // returns position in reduced storage from position in full 3x3 matrix
  int getPosition(const unsigned int n) const
  {return applyTab[n];} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);

    // low rank approx does not support nul kernels
//    if(_i==_j)
      return FReal(1.)/xy.norm();
//    else
//      return FReal(0.);

  }

  FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
  {
    const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
    return getScaleFactor(CellWidth);
  }

  FReal getScaleFactor(const FReal CellWidth) const
  {
        return FReal(2.) / CellWidth;
  }

};

/// R_{,ij}
// PB: IMPORTANT! This matrix kernel does not present the symmetries 
// required by ChebSym kernel => only suited for Unif interpolation
struct FInterpMatrixKernel_R_IJ : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = R_IJ;
  static const unsigned int DIM = 6; //PB: dimension of kernel
  const unsigned int indexTab[12]={0,0,0,1,1,2,
                                   0,1,2,1,2,2};
  static const unsigned int NRHS = 3;
  static const unsigned int NLHS = 3;

  // store positions in sym tensor 
  const unsigned int applyTab[9]={0,1,2,
                                  1,3,4,
                                  2,4,5};

  const unsigned int _i,_j;

  FInterpMatrixKernel_R_IJ(const unsigned int d = 0) 
    : _i(indexTab[d]), _j(indexTab[d+DIM])
  {}

  // returns position in reduced storage from position in full 3x3 matrix
  int getPosition(const unsigned int n) const
  {return applyTab[n];} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    const FReal one_over_r = FReal(1.)/xy.norm();
    const FReal one_over_r3 = one_over_r*one_over_r*one_over_r;
    double ri,rj;
    
    if(_i==0) ri=xy.getX();
    else if(_i==1) ri=xy.getY();
    else if(_i==2) ri=xy.getZ();
    else throw std::runtime_error("Update i!");

    if(_j==0) rj=xy.getX();
    else if(_j==1) rj=xy.getY();
    else if(_j==2) rj=xy.getZ();
    else throw std::runtime_error("Update j!");

//    return xy.getX() * xy.getX() * one_over_r3; //PB: test r0^2/R^3

    if(_i==_j)
      return one_over_r - ri * ri * one_over_r3;
    else
      return - ri * rj * one_over_r3;

  }

  FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
  {
    const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
    return getScaleFactor(CellWidth);
  }

  // R_{,ij} is homogeneous to [L]/[L]^{-2}=[L]^{-1}
  // => scales like ONE_OVER_R
  FReal getScaleFactor(const FReal CellWidth) const
  {
        return FReal(2.) / CellWidth;
  }

};


/// R_{,ijk}
struct FInterpMatrixKernel_R_IJK : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = R_IJK;
  static const  unsigned int DIM = 10; //PB: dimension of kernel
  const unsigned int indexTab[30]={0,0,0,1,1,1,2,2,2,0,
                                   0,1,2,0,1,2,0,1,2,1,
                                   0,1,2,0,1,2,0,1,2,2};
  static const unsigned int NRHS = 3;
  static const unsigned int NLHS = 3;

  // store positions in sym tensor 
  const unsigned int applyTab[27]={0,3,6,3,1,9,6,9,2,
                                   3,1,9,1,4,7,9,7,5,
                                   6,9,2,9,7,5,2,5,8};

  const unsigned int _i,_j,_k;

  FInterpMatrixKernel_R_IJK(const unsigned int d = 0) 
  : _i(indexTab[d]), _j(indexTab[d+DIM]), _k(indexTab[d+2*DIM])
  {}

  // returns position in reduced storage from position in full 3x3x3 matrix
  int getPosition(const unsigned int n) const
  {return applyTab[n];} 

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    const FReal one_over_r = FReal(1.)/xy.norm();
    const FReal one_over_r2 = one_over_r*one_over_r;
    const FReal one_over_r3 = one_over_r2*one_over_r;
    double ri,rj,rk;

    if(_i==0) ri=xy.getX();
    else if(_i==1) ri=xy.getY();
    else if(_i==2) ri=xy.getZ();
    else throw std::runtime_error("Update i!");

    if(_j==0) rj=xy.getX();
    else if(_j==1) rj=xy.getY();
    else if(_j==2) rj=xy.getZ();
    else throw std::runtime_error("Update j!");

    if(_k==0) rk=xy.getX();
    else if(_k==1) rk=xy.getY();
    else if(_k==2) rk=xy.getZ();
    else throw std::runtime_error("Update k!");

    const FReal ri2=ri*ri;
    const FReal rj2=rj*rj;

    if(_i==_j){
      if(_j==_k) //i=j=k
        return FReal(3.) * ( FReal(-1.) + ri2 * one_over_r2 ) * ri * one_over_r3;
      else //i=j!=k
        return ( FReal(-1.) + FReal(3.) * ri2 * one_over_r2 ) * rk * one_over_r3;
    }
    else //(_i!=j)
      if(_i==_k) //i=k!=j
        return ( FReal(-1.) + FReal(3.) * ri2 * one_over_r2 ) * rj * one_over_r3;
      else if(_j==_k) //i!=k=j
        return ( FReal(-1.) + FReal(3.) * rj2 * one_over_r2 ) * ri * one_over_r3;
      else //i!=k!=j
        return FReal(3.) * ri * rj * rk * one_over_r2 * one_over_r3;

  }

  FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
  {
    const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
    return getScaleFactor(CellWidth);
  }

  // R_{,ijk} is homogeneous to [L]/[L]^{-3}=[L]^{-2}
  // => scales like ONE_OVER_RR
  FReal getScaleFactor(const FReal CellWidth) const
  {
    return FReal(4.) / (CellWidth*CellWidth);
  }
};



/*!  Functor which provides the interface to assemble a matrix based on the
  number of rows and cols and on the coordinates x and y and the type of the
  generating matrix-kernel function.
*/
template <typename MatrixKernel>
class EntryComputer
{
  const MatrixKernel Kernel;

  const unsigned int nx, ny;
  const FPoint *const px, *const py;

  const FReal *const weights;

public:
  explicit EntryComputer(const unsigned int _nx, const FPoint *const _px,
                         const unsigned int _ny, const FPoint *const _py,
                         const FReal *const _weights = NULL,
                         const unsigned int idxK = 0)
    : Kernel(idxK),	nx(_nx), ny(_ny), px(_px), py(_py), weights(_weights) {}

  void operator()(const unsigned int xbeg, const unsigned int xend,
                  const unsigned int ybeg, const unsigned int yend,
                  FReal *const data) const
  {
    unsigned int idx = 0;
    if (weights) {
      for (unsigned int j=ybeg; j<yend; ++j)
        for (unsigned int i=xbeg; i<xend; ++i)
          data[idx++] = weights[i] * weights[j] * Kernel.evaluate(px[i], py[j]);
    } else {
      for (unsigned int j=ybeg; j<yend; ++j)
        for (unsigned int i=xbeg; i<xend; ++i)
          data[idx++] = Kernel.evaluate(px[i], py[j]);
    }

    /*
    // apply weighting matrices
    if (weights) {
    if ((xend-xbeg) == (yend-ybeg) && (xend-xbeg) == nx)
    for (unsigned int n=0; n<nx; ++n) {
    FBlas::scal(nx, weights[n], data + n,  nx); // scale rows
    FBlas::scal(nx, weights[n], data + n * nx); // scale cols
    }
    else if ((xend-xbeg) == 1 && (yend-ybeg) == ny)
    for (unsigned int j=0; j<ny; ++j)	data[j] *= weights[j];
    else if ((yend-ybeg) == 1 && (xend-xbeg) == nx)
    for (unsigned int i=0; i<nx; ++i)	data[i] *= weights[i];
    }
    */

  }
};





#endif // FINTERPMATRIXKERNEL_HPP

// [--END--]
