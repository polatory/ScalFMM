// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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

#include "Utils/FPoint.hpp"
#include "Utils/FNoCopyable.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FBlas.hpp"

// extendable
enum KERNEL_FUNCTION_IDENTIFIER {ONE_OVER_R,
                                 ONE_OVER_R_SQUARED,
                                 LENNARD_JONES_POTENTIAL,
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
 * It can either be scalar (NCMP=1) or tensorial (NCMP>1) depending on the 
 * dimension of the equation considered. NCMP denotes the number of components 
 * that are actually stored (e.g. 6 for a \f$3\times3\f$ symmetric tensor). 
 * Notes on application scheme: 
 * Let there be a kernel \f$K\f$ such that \f$X_i=K_{ij}Y_j\f$
 * with \f$X\f$ the lhs of size NLHS and \f$Y\f$ the rhs of size NRHS. 
 * The table applyTab provides the indices in the reduced storage table 
 * corresponding to the application scheme depicted ealier.
 *
 * PB: BEWARE! Homogeneous matrix kernels do not support cell width extension
 * yet. Is it possible to find a reference width and a scale factor such that
 * only 1 set of M2L ops can be used for all levels?? 
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
  static const unsigned int NCMP = 1; //< number of components
  static const unsigned int NPV  = 1; //< dim of physical values
  static const unsigned int NPOT = 1; //< dim of potentials
  static const unsigned int NRHS = 1; //< dim of mult exp
  static const unsigned int NLHS = 1; //< dim of loc exp

  FInterpMatrixKernelR(const FReal = 0.0, const unsigned int = 0) {}

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

  void evaluateBlock(const FPoint& x, const FPoint& y, FReal* block) const
  {
    block[0]=this->evaluate(x,y);
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

/// One over r when the box size is rescaled to 1
struct FInterpMatrixKernelRH :FInterpMatrixKernelR{
	FReal LX,LY,LZ ;
	FInterpMatrixKernelRH(const FReal = 0.0, const unsigned int = 0) : FInterpMatrixKernelR(),
			LX(1.0),LY(1.0),LZ(1.0)
	{}
	  FReal evaluate(const FPoint& x, const FPoint& y) const
	  {
	    const FPoint xy(x-y);
	    return FReal(1.) / FMath::Sqrt(LX*xy.getX()*xy.getX() +
	    		LY*xy.getY()*xy.getY() +
	    		LZ*xy.getZ()*xy.getZ());
	  }
	  	 void setCoeff(const FReal& a,  const FReal& b, const FReal& c)
	  	 {LX= a ; LY = b ; LZ = c ;}
};


/// One over r^2
struct FInterpMatrixKernelRR : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = ONE_OVER_R_SQUARED;
  static const unsigned int NCMP = 1; //< number of components
  static const unsigned int NPV  = 1; //< dim of physical values
  static const unsigned int NPOT = 1; //< dim of potentials
  static const unsigned int NRHS = 1; //< dim of mult exp
  static const unsigned int NLHS = 1; //< dim of loc exp

  FInterpMatrixKernelRR(const FReal = 0.0, const unsigned int = 0) {}

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

  void evaluateBlock(const FPoint& x, const FPoint& y, FReal* block) const
  {
    block[0]=this->evaluate(x,y);
  }

  FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
  {
    const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
    return getScaleFactor(CellWidth);
  }

  FReal getScaleFactor(const FReal CellWidth) const
  {
    return FReal(4.) / (CellWidth*CellWidth);
  }
};



/// One over r^12 - One over r^6
struct FInterpMatrixKernelLJ : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = LENNARD_JONES_POTENTIAL;
  static const unsigned int NCMP = 1; //< number of components
  static const unsigned int NPV  = 1; //< dim of physical values
  static const unsigned int NPOT = 1; //< dim of potentials
  static const unsigned int NRHS = 1; //< dim of mult exp
  static const unsigned int NLHS = 1; //< dim of loc exp

  FInterpMatrixKernelLJ(const FReal = 0.0, const unsigned int = 0) {}

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

  void evaluateBlock(const FPoint& x, const FPoint& y, FReal* block) const
  {
    block[0]=this->evaluate(x,y);
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
// Tensorial Matrix Kernels (NCMP>1)
//
// The definition of the potential p and force f are extended to the case
// of tensorial interaction kernels:
// p_i(x) = K_{ip}(x,y)w_p(y), \forall i=1..NPOT, p=1..NPV
// f_{ik}= w_p(x)K_{ip,k}(x,y)w_p(y) "
//
// Since the interpolation scheme is such that
// p_i(x) \approx S^m(x) L^{m}_{ip}
// f_{ik}= w_p(x) \nabla_k S^m(x) L^{m}_{ip}
// with
// L^{m}_{ip} = K^{mn}_{ip} S^n(y) w_p(y) (local expansion)
// M^{m}_{p} = S^n(y) w_p(y) (multipole expansion)
// then the multipole exp have NPV components and the local exp NPOT*NPV.
//
// NB1: Only the computation of forces requires that the sum over p is 
// performed at L2P step. It could be done at M2L step for the potential.
//
// NB2: An efficient application of the matrix kernel is highly kernel 
// dependent, we recommand overriding the P2M/M2L/L2P function of the kernel 
// you are using in order to have optimal performances + set your own NRHS/NLHS.
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/// R_{,ij}
// PB: IMPORTANT! This matrix kernel does not present the symmetries 
// required by ChebSym kernel => only suited for Unif interpolation
struct FInterpMatrixKernel_R_IJ : FInterpAbstractMatrixKernel
{
  static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = R_IJ;
  static const unsigned int NK   = 3*3; //< total number of components
  static const unsigned int NCMP = 6;   //< number of components
  static const unsigned int NPV  = 3;   //< dim of physical values
  static const unsigned int NPOT = 3;   //< dim of potentials
  static const unsigned int NRHS = NPV; //< dim of mult exp
  static const unsigned int NLHS = NPOT*NRHS; //< dim of loc exp

  // store indices (i,j) corresponding to sym idx
  static const unsigned int indexTab[/*2*NCMP=12*/];

  // store positions in sym tensor (when looping over NRHSxNLHS)
  static const unsigned int applyTab[/*NK=9*/];

  // indices to be set at construction if component-wise evaluation is performed
  const unsigned int _i,_j;

  // Material Parameters
  const FReal _CoreWidth2; // if >0 then kernel is NON homogeneous

  FInterpMatrixKernel_R_IJ(const FReal CoreWidth = 0.0, const unsigned int d = 0)
    : _i(indexTab[d]), _j(indexTab[d+NCMP]), _CoreWidth2(CoreWidth*CoreWidth)
  {}

  // returns position in reduced storage from position in full 3x3 matrix
 unsigned  int getPosition(const unsigned int n) const
  {return applyTab[n];} 

  // returns Core Width squared
  FReal getCoreWidth2() const
  {return _CoreWidth2;}

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    const FPoint xy(x-y);
    const FReal one_over_r = FReal(1.)/sqrt(xy.getX()*xy.getX() + 
                                            xy.getY()*xy.getY() + 
                                            xy.getZ()*xy.getZ() + _CoreWidth2);
    const FReal one_over_r3 = one_over_r*one_over_r*one_over_r;
    FReal ri,rj;
    
    if(_i==0) ri=xy.getX();
    else if(_i==1) ri=xy.getY();
    else if(_i==2) ri=xy.getZ();
    else throw std::runtime_error("Update i!");

    if(_j==0) rj=xy.getX();
    else if(_j==1) rj=xy.getY();
    else if(_j==2) rj=xy.getZ();
    else throw std::runtime_error("Update j!");

    if(_i==_j)
      return one_over_r - ri * ri * one_over_r3;
    else
      return - ri * rj * one_over_r3;

  }

  void evaluateBlock(const FPoint& x, const FPoint& y, FReal* block) const
  {
    const FPoint xy(x-y);
    const FReal one_over_r = FReal(1.)/sqrt(xy.getX()*xy.getX() + 
                                            xy.getY()*xy.getY() + 
                                            xy.getZ()*xy.getZ() + _CoreWidth2);
    const FReal one_over_r3 = one_over_r*one_over_r*one_over_r;
    const FReal r[3] = {xy.getX(),xy.getY(),xy.getZ()};

    for(unsigned int d=0;d<NCMP;++d){
      unsigned int i = indexTab[d];
      unsigned int j = indexTab[d+NCMP];

      if(i==j)
        block[d] = one_over_r - r[i] * r[i] * one_over_r3;
      else
        block[d] = - r[i] * r[j] * one_over_r3;
    }        
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
  static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
  static const KERNEL_FUNCTION_IDENTIFIER Identifier = R_IJK;
  static const unsigned int NK = 3*3*3; //< total number of components
  static const unsigned int NCMP = 10;  //< number of components after symmetry
  static const unsigned int NPOT = 3;   //< dim of potentials
  static const unsigned int NPV  = 3*3; //< dim of physical values
  static const unsigned int NRHS = NPV; //< dim of mult exp
  static const unsigned int NLHS = NPOT*NRHS; //< dim of loc exp

  // store indices (i,j,k) corresponding to sym idx
  static const unsigned int indexTab[/*3*NCMP=30*/];

  // store positions in sym tensor wr to full
  static const unsigned int applyTab[/*NK=27*/];

  // indices to be set at construction if component-wise evaluation is performed
  const unsigned int _i,_j,_k;

  // Material Parameters
  const FReal _CoreWidth2; // if >0 then kernel is NON homogeneous

  FInterpMatrixKernel_R_IJK(const FReal CoreWidth = 0.0, const unsigned int d = 0)
  : _i(indexTab[d]), _j(indexTab[d+NCMP]), _k(indexTab[d+2*NCMP]), _CoreWidth2(CoreWidth*CoreWidth)
  {}

  // returns position in reduced storage from position in full 3x3x3 matrix
  unsigned int getPosition(const unsigned int n) const
  {return applyTab[n];} 

  // returns Core Width squared
  FReal getCoreWidth2() const
  {return _CoreWidth2;}

  FReal evaluate(const FPoint& x, const FPoint& y) const
  {
    // Convention for anti-symmetric kernels xy=y-x instead of xy=x-y
    // This makes the P2P mutual interactions implementation more explicit 
    const FPoint xy(y-x);
    const FReal one_over_r = FReal(1.)/sqrt(xy.getX()*xy.getX() + 
                                            xy.getY()*xy.getY() + 
                                            xy.getZ()*xy.getZ() + _CoreWidth2);
    const FReal one_over_r2 = one_over_r*one_over_r;
    const FReal one_over_r3 = one_over_r2*one_over_r;
    FReal ri,rj,rk;

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

  void evaluateBlock(const FPoint& x, const FPoint& y, FReal* block) const
  {
    // Convention for anti-symmetric kernels xy=y-x instead of xy=x-y
    // This makes the P2P mutual interactions implementation more explicit 
    const FPoint xy(y-x);
    const FReal one_over_r = FReal(1.)/sqrt(xy.getX()*xy.getX() + 
                                            xy.getY()*xy.getY() + 
                                            xy.getZ()*xy.getZ() + _CoreWidth2);
    const FReal one_over_r2 = one_over_r*one_over_r;
    const FReal one_over_r3 = one_over_r2*one_over_r;

    const FReal r[3] = {xy.getX(),xy.getY(),xy.getZ()};

    for(unsigned int d=0;d<NCMP;++d){
      unsigned int i = indexTab[d];
      unsigned int j = indexTab[d+NCMP];
      unsigned int k = indexTab[d+2*NCMP];

      if(i==j){
        if(j==k) //i=j=k
          block[d] = FReal(3.) * ( FReal(-1.) + r[i]*r[i] * one_over_r2 ) * r[i] * one_over_r3;
        else //i=j!=k
          block[d] =  ( FReal(-1.) + FReal(3.) * r[i]*r[i] * one_over_r2 ) * r[k] * one_over_r3;
      }
      else //(i!=j)
        if(i==k) //i=k!=j
          block[d] =  ( FReal(-1.) + FReal(3.) * r[i]*r[i] * one_over_r2 ) * r[j] * one_over_r3;
        else if(j==k) //i!=k=j
          block[d] =  ( FReal(-1.) + FReal(3.) * r[j]*r[j] * one_over_r2 ) * r[i] * one_over_r3;
        else //i!=k!=j
          block[d] =  FReal(3.) * r[i]*r[j]*r[k] * one_over_r2 * one_over_r3;
    }      
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
                         const unsigned int idxK = 0,
                         const FReal matparam = 0.0)
    : Kernel(matparam,idxK),	nx(_nx), ny(_ny), px(_px), py(_py), weights(_weights) {}

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
