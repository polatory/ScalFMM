#ifndef FCHEBSYMKERNEL_HPP
#define FCHEBSYMKERNEL_HPP
// [--License--]

#include "../../Utils/FGlobal.hpp"

#include "../../Utils/FSmartPointer.hpp"

#include "./FAbstractChebKernel.hpp"
#include "./FChebInterpolator.hpp"
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
//#include "./FChebSymM2LHandler.hpp"
#include "./FChebSymTensorialM2LHandler.hpp" //PB: temporary version

class FTreeCoordinate;


// for verbosity only!!!
//#define COUNT_BLOCKED_INTERACTIONS

// if timings should be logged
//#define LOG_TIMINGS

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebSymTensorialKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Chebyshev interpolation based FMM operators
 * exploiting the symmetries in the far-field. It implements all interfaces
 * (P2P, P2M, M2M, M2L, L2L, L2P) which are required by the FFmmAlgorithm and
 * FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
class FChebSymTensorialKernel
  : public FAbstractChebKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
  typedef FAbstractChebKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	AbstractBaseClass;
  typedef SymmetryHandler<ORDER, MatrixKernelClass::NCMP, MatrixKernelClass::Type> SymmetryHandlerClass;
  enum {nnodes = AbstractBaseClass::nnodes,
        ncmp = MatrixKernelClass::NCMP,
        nRhs = MatrixKernelClass::NRHS,
        nLhs = MatrixKernelClass::NLHS};

  /// Needed for handling all symmetries
  const FSmartPointer<SymmetryHandlerClass,FSmartPointerMemory> SymHandler;

  // permuted local and multipole expansions
  FReal** Loc;
  FReal** Mul;
  unsigned int* countExp;



  /**
   * Allocate memory for storing locally permuted mulipole and local expansions
   */
  void allocateMemoryForPermutedExpansions()
  {
    assert(Loc==NULL && Mul==NULL && countExp==NULL);
    Loc = new FReal* [343];
    Mul = new FReal* [343];
    countExp = new unsigned int [343];

    // set all 343 to NULL
    for (unsigned int idx=0; idx<343; ++idx) {
      Mul[idx] = Loc[idx] = NULL;
    }

    // init only 16 of 343 possible translations due to symmetries
    for (int i=2; i<=3; ++i)
      for (int j=0; j<=i; ++j)
        for (int k=0; k<=j; ++k) {
          const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
          assert(Mul[idx]==NULL || Loc[idx]==NULL);
          Mul[idx] = new FReal [24 * nnodes];
          Loc[idx] = new FReal [24 * nnodes];
        }
  }


#ifdef LOG_TIMINGS
  FTic time;
  FReal t_m2l_1, t_m2l_2, t_m2l_3;
#endif

public:
  /**
   * The constructor initializes all constant attributes and it reads the
   * precomputed and compressed M2L operators from a binary file (an
   * runtime_error is thrown if the required file is not valid).
   */
  FChebSymTensorialKernel(const int inTreeHeight,
                 const FReal inBoxWidth,
                 const FPoint& inBoxCenter,
                 const FReal Epsilon)
    : AbstractBaseClass(inTreeHeight, inBoxWidth, inBoxCenter),
      SymHandler(new SymmetryHandlerClass(AbstractBaseClass::MatrixKernel.getPtr(), Epsilon, inBoxWidth, inTreeHeight)),
      Loc(NULL), Mul(NULL), countExp(NULL)
  {
    this->allocateMemoryForPermutedExpansions();

#ifdef LOG_TIMINGS
    t_m2l_1 = FReal(0.);
    t_m2l_2 = FReal(0.);
    t_m2l_3 = FReal(0.);
#endif
  }



  /** Copy constructor */
  FChebSymTensorialKernel(const FChebSymTensorialKernel& other)
  : AbstractBaseClass(other),
    SymHandler(other.SymHandler),
    Loc(NULL), Mul(NULL), countExp(NULL)
  {
    this->allocateMemoryForPermutedExpansions();
  }



  /** Destructor */
  ~FChebSymTensorialKernel()
  {
    for (unsigned int t=0; t<343; ++t) {
      if (Loc[t]!=NULL) delete [] Loc[t];
      if (Mul[t]!=NULL) delete [] Mul[t];
    }
    if (Loc!=NULL)      delete [] Loc;
    if (Mul!=NULL)      delete [] Mul;
    if (countExp!=NULL) delete [] countExp;

#ifdef LOG_TIMINGS
    std::cout << "- Permutation took " << t_m2l_1 << "s"
              << "\n- GEMMT and GEMM took " << t_m2l_2 << "s"
              << "\n- Unpermutation took " << t_m2l_3 << "s"
              << std::endl;
#endif
  }

  // PB: Only used in testChebAlgorithm
  const SymmetryHandlerClass *const getPtrToSymHandler() const
  {	return SymHandler.getPtr();	}



  void P2M(CellClass* const LeafCell,
           const ContainerClass* const SourceParticles)
  {
    const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){
        // update multipole index
        int idxMul = idxV*nRhs + idxRhs;
        // apply Sy
        AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                  LeafCell->getMultipole(idxMul), SourceParticles);
      }// NRHS
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
        // apply Sy
        FBlas::scal(nnodes*2, FReal(0.), ParentCell->getMultipole(idxMul));
        for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
          if (ChildCells[ChildIndex]){
            AbstractBaseClass::Interpolator->applyM2M(ChildIndex, ChildCells[ChildIndex]->getMultipole(idxMul), ParentCell->getMultipole(idxMul));
          }
        }
      }// NRHS
    }// NVALS
  }


  void M2L(CellClass* const FRestrict TargetCell,
           const CellClass* SourceCells[343],
           const int /*NumSourceCells*/,
           const int TreeLevel)
  {
#ifdef LOG_TIMINGS
    time.tic();
#endif
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for (int idxLhs=0; idxLhs < nLhs; ++idxLhs){
        // update local index
        int idxLoc = idxV*nLhs + idxLhs;

        for (int idxRhs=0; idxRhs < nRhs; ++idxRhs){
          // update multipole index
          int idxMul = idxV*nRhs + idxRhs;
          // update kernel index such that: x_i = K_{ij}y_j 
          int idxK = idxLhs*nRhs + idxRhs;
          unsigned int d = AbstractBaseClass::MatrixKernel->getPosition(idxK);

          // permute and copy multipole expansion
          memset(countExp, 0, sizeof(int) * 343);
          for (unsigned int idx=0; idx<343; ++idx) {
            if (SourceCells[idx]) {
              const unsigned int pidx = SymHandler->pindices[idx];
              const unsigned int count = (countExp[pidx])++;
              FReal *const mul = Mul[pidx] + count*nnodes;
              const unsigned int *const pvec = SymHandler->pvectors[idx];
              const FReal *const MultiExp = SourceCells[idx]->getMultipole(idxMul);

              /*
              // no loop unrolling
              for (unsigned int n=0; n<nnodes; ++n)
              mul[pvec[n]] = MultiExp[n];
              */

              // explicit loop unrolling
              for (unsigned int n=0; n<nnodes/4 * 4; n+=4) {
                mul[pvec[n  ]] = MultiExp[n];
                mul[pvec[n+1]] = MultiExp[n+1];
                mul[pvec[n+2]] = MultiExp[n+2];
                mul[pvec[n+3]] = MultiExp[n+3];
              }
              for (unsigned int n=nnodes/4 * 4; n<nnodes; ++n){
                mul[pvec[n]] = MultiExp[n];
              }
            }
          }

#ifdef LOG_TIMINGS
          t_m2l_1 += time.tacAndElapsed();
#endif


#ifdef COUNT_BLOCKED_INTERACTIONS ////////////////////////////////////
          unsigned int count_lidx = 0;
          unsigned int count_interactions = 0;
          for (unsigned int idx=0; idx<343; ++idx)
            count_interactions += countExp[idx];
          if (count_interactions==189) {
            for (unsigned int idx=0; idx<343; ++idx) {
              if (countExp[idx])
                std::cout << "gidx = " << idx << " gives lidx = " << count_lidx++ << " and has "
                          << countExp[idx] << " interactions" << std::endl;
            }
            std::cout << std::endl;
          }
#endif ///////////////////////////////////////////////////////////////


#ifdef LOG_TIMINGS
          time.tic();
#endif

          // multiply (mat-mat-mul)
          FReal Compressed [nnodes * 24];
          const FReal scale = AbstractBaseClass::MatrixKernel->getScaleFactor(AbstractBaseClass::BoxWidth, TreeLevel);
          for (unsigned int pidx=0; pidx<343; ++pidx) {
            const unsigned int count = countExp[pidx];
            if (count) {
              const unsigned int rank = SymHandler->getLowRank(TreeLevel, pidx, d);
              // rank * count * (2*nnodes-1) flops
              FBlas::gemtm(nnodes, rank, count, FReal(1.),
                           const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx, d))+rank*nnodes,
                           nnodes, Mul[pidx], nnodes, Compressed, rank);
              // nnodes *count * (2*rank-1) flops
              FBlas::gemm( nnodes, rank, count, scale,
                           const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx, d)),
                           nnodes, Compressed, rank, Loc[pidx], nnodes);
            }
          }

#ifdef LOG_TIMINGS
          t_m2l_2 += time.tacAndElapsed();
#endif

#ifdef LOG_TIMINGS
          time.tic();
#endif

          // permute and add contribution to local expansions
          FReal *const LocalExpansion = TargetCell->getLocal(idxLoc);
          memset(countExp, 0, sizeof(int) * 343);
          for (unsigned int idx=0; idx<343; ++idx) {
            if (SourceCells[idx]) {
              const unsigned int pidx = SymHandler->pindices[idx];
              const unsigned int count = (countExp[pidx])++;
              const FReal *const loc = Loc[pidx] + count*nnodes;
              const unsigned int *const pvec = SymHandler->pvectors[idx];

              /*
              // no loop unrolling
              for (unsigned int n=0; n<nnodes; ++n)
              LocalExpansion[n] += loc[pvec[n]];
              */

              // explicit loop unrolling
              for (unsigned int n=0; n<nnodes/4 * 4; n+=4) {
                LocalExpansion[n  ] += loc[pvec[n  ]];
                LocalExpansion[n+1] += loc[pvec[n+1]];
                LocalExpansion[n+2] += loc[pvec[n+2]];
                LocalExpansion[n+3] += loc[pvec[n+3]];
              }
              for (unsigned int n=nnodes/4 * 4; n<nnodes; ++n)
                LocalExpansion[n] += loc[pvec[n]];

            }
          }

#ifdef LOG_TIMINGS
          t_m2l_3 += time.tacAndElapsed();
#endif
        }// NRHS
      }// NLHS
    }// NVALS
  }

  



  /*
   * M2L: More basic version
   */
  /*
  void M2L(CellClass* const FRestrict TargetCell,
           const CellClass* SourceCells[343],
           const int NumSourceCells,
           const int TreeLevel)
  {

    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for (int idxLhs=0; idxLhs < nLhs; ++idxLhs){
        // update local index
        int idxLoc = idxV*nLhs + idxLhs;

        for (int idxRhs=0; idxRhs < nRhs; ++idxRhs){
          // update multipole index
          int idxMul = idxV*nRhs + idxRhs;
        
          // update kernel index such that: x_i = K_{ij}y_j 
          int idxK = idxLhs*nRhs + idxRhs;
          unsigned int d = AbstractBaseClass::MatrixKernel->getPosition(idxK);

          // permute and copy multipole expansion
          memset(countExp, 0, sizeof(int) * 343);
          for (unsigned int idx=0; idx<343; ++idx) {
            if (SourceCells[idx]) {
              const unsigned int pidx = SymHandler->pindices[idx];
              const unsigned int count = (countExp[pidx])++;
              const FReal *const MultiExp = SourceCells[idx]->getMultipole(idxMul);
              for (unsigned int n=0; n<nnodes; ++n)
                Mul[pidx][count*nnodes + SymHandler->pvectors[idx][n]] = MultiExp[n];
            }
          }

          // multiply (mat-mat-mul)
          FReal Compressed [nnodes * 30];
          const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(FMath::pow(2, TreeLevel)));
          const FReal scale(AbstractBaseClass::MatrixKernel->getScaleFactor(CellWidth));
          for (unsigned int pidx=0; pidx<343; ++pidx) {
            const unsigned int count = countExp[pidx];
            if (count) {
              const unsigned int rank = SymHandler->getLowRank(TreeLevel, pidx, d);
              FBlas::gemtm(nnodes, rank, count, FReal(1.),
                           const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx, d)+rank*nnodes),
                           nnodes, Mul[pidx], nnodes, Compressed, rank);
              FBlas::gemm( nnodes, rank, count, scale,
                           const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx, d)),
                           nnodes, Compressed, rank, Loc[pidx], nnodes);
            }
          }

          // permute and add contribution to local expansions
          FReal *const LocalExpansion = TargetCell->getLocal(idxLoc);
          memset(countExp, 0, sizeof(int) * 343);
          for (unsigned int idx=0; idx<343; ++idx) {
            if (SourceCells[idx]) {
              const unsigned int pidx = SymHandler->pindices[idx];
              const unsigned int count = (countExp[pidx])++;
              for (unsigned int n=0; n<nnodes; ++n)
                LocalExpansion[n] += Loc[pidx][count*nnodes + SymHandler->pvectors[idx][n]];
            }
          }

        }// NRHS
      }// NLHS
    }// NVALS

  }
  */

  void L2L(const CellClass* const FRestrict ParentCell,
           CellClass* FRestrict *const FRestrict ChildCells,
           const int /*TreeLevel*/)
  {
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
        int idxLoc = idxV*nLhs + idxLhs;
        // apply Sx
        for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
          if (ChildCells[ChildIndex]){
            AbstractBaseClass::Interpolator->applyL2L(ChildIndex, ParentCell->getLocal(idxLoc), ChildCells[ChildIndex]->getLocal(idxLoc));
          }
        }
      }// NLHS
    }// NVALS
  }


  void L2P(const CellClass* const LeafCell,
           ContainerClass* const TargetParticles)
  {
    const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));
    for(int idxV = 0 ; idxV < NVALS ; ++idxV){
      for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
        int idxLoc = idxV*nLhs + idxLhs;
//        // a) apply Sx
//        AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter,
//                                                  AbstractBaseClass::BoxWidthLeaf,
//                                                  LeafCell->getLocal(idxLoc),
//                                                  TargetParticles);
//        // b) apply Px (grad Sx)
//        AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter,
//                                                          AbstractBaseClass::BoxWidthLeaf,
//                                                          LeafCell->getLocal(idxLoc),
//                                                          TargetParticles);

        // c) apply Sx and Px (grad Sx)
        AbstractBaseClass::Interpolator->applyL2PTotal(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                       LeafCell->getLocal(idxLoc), TargetParticles);
      }// NLHS
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








#endif //FCHEBSYMTENSORIALKERNELS_HPP

// [--END--]
