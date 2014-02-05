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
#ifndef FUNIFTENSORIALKERNEL_HPP
#define FUNIFTENSORIALKERNEL_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FSmartPointer.hpp"

#include "./FAbstractUnifKernel.hpp"
#include "./FUnifM2LHandler.hpp"
#include "./FUnifTensorialM2LHandler.hpp" //PB: temporary version

class FTreeCoordinate;

/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FUnifTensorialKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Lagrange interpolation based FMM operators. It
 * implements all interfaces (P2P,P2M,M2M,M2L,L2L,L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * PB: 3 IMPORTANT remarks !!!
 *
 * 1) Handling tensorial kernels (DIM,NRHS,NLHS) and having multiple rhs (NVALS) 
 * are considered 2 separate features and are currently combined.
 *
 * 2) The present tensorial version is the most naive one. All tensorial aspects 
 * are handled in the kernel. A more optimal version would be to consider looping 
 * over nRhs/nLhs inside Interpolator::applyP2M/L2P in order to avoid extra 
 * evaluation of the interpolating polynomials. When it comes to applying M2L it is 
 * NOT much faster to loop over NRHSxNLHS inside applyM2L (at least for the Lagrange case)
 * 2-bis) The evaluation of the kernel matrix (see M2LHandler) should be done at once 
 * instead of compo-by-compo (TODO). On the other hand, the ChebyshevSym tensorial kernel 
 * requires the matrix kernel to be evaluated compo-by-compo since we currently use a scalar ACA.
 *
 * 3) We currently use multiple 1D FFT instead of multidim FFT. 
 * TODO switch to multidim if relevant in considered range of size 
 * (see testFFTW and testFFTWMultidim).
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Lagrange interpolation order
 */
template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
class FUnifTensorialKernel
  : public FAbstractUnifKernel< CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
  enum {nRhs = MatrixKernelClass::NRHS,
        nLhs = MatrixKernelClass::NLHS};

protected://PB: for OptiDis

  // private types
  typedef FUnifTensorialM2LHandler<ORDER,MatrixKernelClass> M2LHandlerClass;

  // using from
  typedef FAbstractUnifKernel< CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
  AbstractBaseClass;

  /// Needed for M2L operator
  FSmartPointer<  M2LHandlerClass,FSmartPointerMemory> M2LHandler;

public:
  /**
   * The constructor initializes all constant attributes and it reads the
   * precomputed and compressed M2L operators from a binary file (an
   * runtime_error is thrown if the required file is not valid).
   */
  FUnifTensorialKernel(const int inTreeHeight,
                       const FReal inBoxWidth,
                       const FPoint& inBoxCenter)
    : FAbstractUnifKernel< CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>(inTreeHeight,
                                                                                       inBoxWidth,
                                                                                       inBoxCenter),
      M2LHandler(new M2LHandlerClass())
  {
    // read precomputed compressed m2l operators from binary file
    //M2LHandler->ReadFromBinaryFileAndSet(); // PB: TODO?
    M2LHandler->ComputeAndSet();
  }


  void P2M(CellClass* const LeafCell,
           const ContainerClass* const SourceParticles)
  {
    const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate())); 
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){
        // update multipole index
        int idxMul = idxV*nRhs + idxRhs;

        // 1) apply Sy
        AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                  LeafCell->getMultipole(idxMul), SourceParticles);

        // 2) apply Discrete Fourier Transform
        M2LHandler->applyZeroPaddingAndDFT(LeafCell->getMultipole(idxMul), 
                                           LeafCell->getTransformedMultipole(idxMul));

      }
    }// NVALS
  }


  void M2M(CellClass* const FRestrict ParentCell,
           const CellClass*const FRestrict *const FRestrict ChildCells,
           const int /*TreeLevel*/)
  {
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){
        // update multipole index
        int idxMul = idxV*nRhs + idxRhs;

        // 1) apply Sy
        FBlas::scal(AbstractBaseClass::nnodes, FReal(0.), ParentCell->getMultipole(idxMul));
        for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
          if (ChildCells[ChildIndex]){
            AbstractBaseClass::Interpolator->applyM2M(ChildIndex, ChildCells[ChildIndex]->getMultipole(idxMul),
                                                      ParentCell->getMultipole(idxMul));
          }
        }
        // 2) Apply Discete Fourier Transform
        M2LHandler->applyZeroPaddingAndDFT(ParentCell->getMultipole(idxMul), 
                                           ParentCell->getTransformedMultipole(idxMul));
      }
    }// NVALS
  }


  void M2L(CellClass* const FRestrict TargetCell,
           const CellClass* SourceCells[343],
           const int /*NumSourceCells*/,
           const int TreeLevel)
  {
    const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(FMath::pow(2, TreeLevel)));

    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for (int idxLhs=0; idxLhs < nLhs; ++idxLhs){
        // update local index
        int idxLoc = idxV*nLhs + idxLhs;
        // load transformed local expansion
        FComplexe *const TransformedLocalExpansion = TargetCell->getTransformedLocal(idxLoc);

        for (int idxRhs=0; idxRhs < nRhs; ++idxRhs){
          // update multipole index
          int idxMul = idxV*nRhs + idxRhs;
          // update kernel index such that: x_i = K_{ij}y_j 
          int idxK = idxLhs*nRhs + idxRhs;

          for (int idx=0; idx<343; ++idx){
            if (SourceCells[idx]){
              M2LHandler->applyFC(idx, CellWidth, 
                                  SourceCells[idx]->getTransformedMultipole(idxMul),
                                  TransformedLocalExpansion,idxK);

            }
          }
        }// NRHS
      }// NLHS
    }// NVALS
  }


  void L2L(const CellClass* const FRestrict ParentCell,
           CellClass* FRestrict *const FRestrict ChildCells,
           const int /*TreeLevel*/)
  {
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
        int idxLoc = idxV*nLhs + idxLhs;
        // 1) Apply Inverse Discete Fourier Transform
        M2LHandler->unapplyZeroPaddingAndDFT(ParentCell->getTransformedLocal(idxLoc),
                                             const_cast<CellClass*>(ParentCell)->getLocal(idxLoc));
        // 2) apply Sx
        for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
          if (ChildCells[ChildIndex]){
            AbstractBaseClass::Interpolator->applyL2L(ChildIndex, ParentCell->getLocal(idxLoc), ChildCells[ChildIndex]->getLocal(idxLoc));
          }
        }
      }
    }// NVALS
  }

  void L2P(const CellClass* const LeafCell,
           ContainerClass* const TargetParticles)
  {
    const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));

    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
        int idxLoc = idxV*nLhs + idxLhs;
        // 1)  Apply Inverse Discete Fourier Transform
        M2LHandler->unapplyZeroPaddingAndDFT(LeafCell->getTransformedLocal(idxLoc), 
                                             const_cast<CellClass*>(LeafCell)->getLocal(idxLoc));

        // 2.a) apply Sx
        AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                  LeafCell->getLocal(idxLoc), TargetParticles);

        // 2.b) apply Px (grad Sx)
        AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                          LeafCell->getLocal(idxLoc), TargetParticles);

      }
    }// NVALS
  }

  void P2P(const FTreeCoordinate& /* LeafCellCoordinate */, // needed for periodic boundary conditions
           ContainerClass* const FRestrict TargetParticles,
           const ContainerClass* const FRestrict /*SourceParticles*/,
           ContainerClass* const NeighborSourceParticles[27],
           const int /* size */)
  {
    DirectInteractionComputer<MatrixKernelClass::Identifier, NVALS>::P2P(TargetParticles,NeighborSourceParticles);
  }


  void P2PRemote(const FTreeCoordinate& /*inPosition*/,
                 ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
                 ContainerClass* const inNeighbors[27], const int /*inSize*/){
    DirectInteractionComputer<MatrixKernelClass::Identifier, NVALS>::P2PRemote(inTargets,inNeighbors,27);
  }

};


#endif //FUNIFTENSORIALKERNEL_HPP

// [--END--]