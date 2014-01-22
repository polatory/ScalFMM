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
#ifndef FUNIFSYMKERNEL_HPP
#define FUNIFSYMKERNEL_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FSmartPointer.hpp"

// Originally in M2LHandler but transferred to the kernel for the symmetric version
#include "../../Utils/FDft.hpp" // PB: for FFT
#include "../../Utils/FComplexe.hpp"
#include "./FUnifTensor.hpp" // PB: for node_diff
//

#include "./FAbstractUnifKernel.hpp"
#include "./FUnifInterpolator.hpp"

#include "./FUnifSymM2LHandler.hpp"

class FTreeCoordinate;


// for verbosity only!!!
//#define COUNT_BLOCKED_INTERACTIONS

// if timings should be logged
#define LOG_TIMINGS

/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FUnifSymKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Lagrange interpolation based FMM operators
 * exploiting the symmetries in the far-field. It implements all interfaces
 * (P2P, P2M, M2M, M2L, L2L, L2P) which are required by the FFmmAlgorithm and
 * FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER interpolation order
 */
template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
class FUnifSymKernel
  : public FAbstractUnifKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
  typedef FAbstractUnifKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	AbstractBaseClass;
  typedef FUnifSymM2LHandler<ORDER, MatrixKernelClass::Type> SymM2LHandlerClass;
  enum {nnodes = AbstractBaseClass::nnodes};

  /// Needed for handling all symmetries
  const FSmartPointer<SymM2LHandlerClass,FSmartPointerMemory> SymM2LHandler;

  // permuted local and multipole expansions
  FReal** Loc;
  FReal** Mul;
  unsigned int* countExp;

  // transformed expansions
  FComplexe** TLoc;
  FComplexe** TMul;

  static const unsigned int rc = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);
  static const unsigned int opt_rc = rc/2+1;

  typedef FUnifTensor<ORDER> TensorType;
  unsigned int node_diff[nnodes*nnodes];

  //  FDft Dft; // Direct Discrete Fourier Transformator
  FFft Dft; // Fast Discrete Fourier Transformator

  /**
   * Allocate memory for storing locally permuted mulipole and local expansions
   */
  void allocateMemoryForPermutedExpansions()
  {
    assert(Loc==NULL && Mul==NULL && countExp==NULL);
    Loc = new FReal* [343];
    Mul = new FReal* [343];
    TLoc = new FComplexe* [343];
    TMul = new FComplexe* [343];
    countExp = new unsigned int [343];

    // set all 343 to NULL
    for (unsigned int idx=0; idx<343; ++idx) {
      Mul[idx] = Loc[idx] = NULL;
      TMul[idx] = TLoc[idx] = NULL;
    }

    // init only 16 of 343 possible translations due to symmetries
    for (int i=2; i<=3; ++i)
      for (int j=0; j<=i; ++j)
        for (int k=0; k<=j; ++k) {
          const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
          assert(Mul[idx]==NULL || Loc[idx]==NULL);
          Mul[idx] = new FReal [24 * nnodes];
          Loc[idx] = new FReal [24 * nnodes];
          TMul[idx] = new FComplexe [24 * rc];
          TLoc[idx] = new FComplexe [24 * rc];
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
  FUnifSymKernel(const int inTreeHeight,
                 const FPoint& inBoxCenter,
                 const FReal inBoxWidth)
    : AbstractBaseClass(inTreeHeight, inBoxCenter, inBoxWidth),
      SymM2LHandler(new SymM2LHandlerClass(AbstractBaseClass::MatrixKernel.getPtr(), inBoxWidth, inTreeHeight)),
      Loc(NULL), Mul(NULL), countExp(NULL),
      TLoc(NULL), TMul(NULL),
      Dft(rc) // initialize Discrete Fourier Transformator,
  {
    this->allocateMemoryForPermutedExpansions();

    // initialize root node ids
    TensorType::setNodeIdsDiff(node_diff);

#ifdef LOG_TIMINGS
    t_m2l_1 = FReal(0.);
    t_m2l_2 = FReal(0.);
    t_m2l_3 = FReal(0.);
#endif
  }



  /** Copy constructor */
  FUnifSymKernel(const FUnifSymKernel& other)
  : AbstractBaseClass(other),
    SymM2LHandler(other.SymM2LHandler),
    Loc(NULL), Mul(NULL), countExp(NULL),
    TLoc(NULL), TMul(NULL)
  {
    this->allocateMemoryForPermutedExpansions();
  }



  /** Destructor */
  ~FUnifSymKernel()
  {
    for (unsigned int t=0; t<343; ++t) {
      if (Loc[t]!=NULL) delete [] Loc[t];
      if (Mul[t]!=NULL) delete [] Mul[t];
      if (TLoc[t]!=NULL) delete [] TLoc[t];
      if (TMul[t]!=NULL) delete [] TMul[t];
    }
    if (Loc!=NULL)      delete [] Loc;
    if (Mul!=NULL)      delete [] Mul;
    if (countExp!=NULL) delete [] countExp;
    if (TLoc!=NULL)      delete [] TLoc;
    if (TMul!=NULL)      delete [] TMul;

#ifdef LOG_TIMINGS
    std::cout << "- Permutation+Pad+Transfo took " << t_m2l_1 << "s"
              << "\n- Apply M2L in Fourier space took " << t_m2l_2 << "s"
              << "\n- Unpermutation+Untransfo+Unpad took " << t_m2l_3 << "s"
              << std::endl;
#endif
  }


  const SymM2LHandlerClass *const getPtrToSymM2LHandler() const
  {	return SymM2LHandler.getPtr();	}



  void P2M(CellClass* const LeafCell,
           const ContainerClass* const SourceParticles)
  {
    const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      // 1) apply Sy
      AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                LeafCell->getMultipole(idxRhs), SourceParticles);
//      // 2) apply Discrete Fourier Transform
//      SymM2LHandler->applyZeroPaddingAndDFT(LeafCell->getMultipole(idxRhs),
//                                         LeafCell->getTransformedMultipole(idxRhs));
    }
  }



  void M2M(CellClass* const FRestrict ParentCell,
           const CellClass*const FRestrict *const FRestrict ChildCells,
           const int /*TreeLevel*/)
  {
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
      // apply Sy
      FBlas::scal(nnodes, FReal(0.), ParentCell->getMultipole(idxRhs));
      for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
        if (ChildCells[ChildIndex]){
          AbstractBaseClass::Interpolator->applyM2M(ChildIndex, ChildCells[ChildIndex]->getMultipole(idxRhs), ParentCell->getMultipole(idxRhs));
        }
      }
//      // 2) Apply Discete Fourier Transform
//      SymM2LHandler->applyZeroPaddingAndDFT(ParentCell->getMultipole(idxRhs),
//                                         ParentCell->getTransformedMultipole(idxRhs));
    }
  }



  void M2L(CellClass* const FRestrict TargetCell,
           const CellClass* SourceCells[343],
           const int /*NumSourceCells*/,
           const int TreeLevel)
  {
#ifdef LOG_TIMINGS
    time.tic();
#endif
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

      // permute and copy multipole expansion
      memset(countExp, 0, sizeof(int) * 343);
      for (unsigned int idx=0; idx<343; ++idx) {
        if (SourceCells[idx]) {
          const unsigned int pidx = SymM2LHandler->pindices[idx];
          const unsigned int count = (countExp[pidx])++;
          FReal *const mul = Mul[pidx] + count*nnodes;
          const unsigned int *const pvec = SymM2LHandler->pvectors[idx];
          const FReal *const MultiExp = SourceCells[idx]->getMultipole(idxRhs);

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

          // transform permuted expansion
          FComplexe *const tmul = TMul[pidx] + count*rc;

          ///////////////////////////////////////////
          FReal pmul[rc];
          FBlas::setzero(rc,pmul);

          // Apply Zero Padding
          for (unsigned int i=0; i<nnodes; ++i)
            pmul[node_diff[i*nnodes]]=mul[i];

          // Apply forward Discrete Fourier Transform
          Dft.applyDFT(pmul,tmul);
          ///////////////////////////////////////////


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
      const FReal scale = AbstractBaseClass::MatrixKernel->getScaleFactor(AbstractBaseClass::BoxWidth, TreeLevel);

      FComplexe tmpTLoc; 

      for (unsigned int pidx=0; pidx<343; ++pidx) {
        const unsigned int count = countExp[pidx];
        if (count) {

          FComplexe *const tmpK=const_cast<FComplexe*>(SymM2LHandler->getK(TreeLevel, pidx));

          for (unsigned int i=0; i<count; ++i){
          //unsigned int i=count;
            FComplexe *const ttmul = TMul[pidx] + i*rc;
            FComplexe *const ttloc = TLoc[pidx] + i*rc;

            // Perform entrywise product manually
            for (unsigned int j=0; j<opt_rc; ++j){
              tmpTLoc=tmpK[j];
              tmpTLoc*=ttmul[j];
              tmpTLoc*=scale;
              ttloc[j]=tmpTLoc;
            }
          }

//          // rank * count * (2*nnodes-1) flops
//          FBlas::gemtm(nnodes, rank, count, FReal(1.),
//                       const_cast<FReal*>(SymM2LHandler->getK(TreeLevel, pidx))+rank*nnodes,
//                       nnodes, Mul[pidx], nnodes, Compressed, rank);
//          // nnodes *count * (2*rank-1) flops
//          FBlas::gemm( nnodes, rank, count, scale,
//                       const_cast<FReal*>(SymM2LHandler->getK(TreeLevel, pidx)),
//                       nnodes, Compressed, rank, Loc[pidx], nnodes);
        }
      }



#ifdef LOG_TIMINGS
      t_m2l_2 += time.tacAndElapsed();
#endif

#ifdef LOG_TIMINGS
      time.tic();
#endif

      // permute and add contribution to local expansions
      FReal *const LocalExpansion = TargetCell->getLocal(idxRhs);
      memset(countExp, 0, sizeof(int) * 343);
      for (unsigned int idx=0; idx<343; ++idx) {
        if (SourceCells[idx]) {
          const unsigned int pidx = SymM2LHandler->pindices[idx];
          const unsigned int count = (countExp[pidx])++;

          FReal *const loc = Loc[pidx] + count*nnodes;
          const FComplexe *const tloc = TLoc[pidx] + count*rc;

          ///////////////////////////////////////////
          FReal ploc[rc];
          FBlas::setzero(rc,ploc);
          // Apply forward Discrete Fourier Transform
          Dft.applyIDFT(tloc,ploc);

          // Unapply Zero Padding
          for (unsigned int j=0; j<nnodes; ++j)
            loc[j]=ploc[node_diff[nnodes-j-1]];
          //loc[j]+=ploc[node_diff[nnodes-j-1]];
          ///////////////////////////////////////////

          const unsigned int *const pvec = SymM2LHandler->pvectors[idx];

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
    }

  }



  void L2L(const CellClass* const FRestrict ParentCell,
           CellClass* FRestrict *const FRestrict ChildCells,
           const int /*TreeLevel*/)
  {
    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
//      // 1) Apply Inverse Discete Fourier Transform
//      SymM2LHandler->unapplyZeroPaddingAndDFT(ParentCell->getTransformedLocal(idxRhs),
//                                           const_cast<CellClass*>(ParentCell)->getLocal(idxRhs));
      // 2) apply Sx
      for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
        if (ChildCells[ChildIndex]){
          AbstractBaseClass::Interpolator->applyL2L(ChildIndex, ParentCell->getLocal(idxRhs), ChildCells[ChildIndex]->getLocal(idxRhs));
        }
      }
    }
  }



  void L2P(const CellClass* const LeafCell,
           ContainerClass* const TargetParticles)
  {
    const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));

    for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

//      // 1)  Apply Inverse Discete Fourier Transform
//      SymM2LHandler->unapplyZeroPaddingAndDFT(LeafCell->getTransformedLocal(idxRhs),
//                                           const_cast<CellClass*>(LeafCell)->getLocal(idxRhs));
//
      // 2.a) apply Sx
      AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                LeafCell->getLocal(idxRhs), TargetParticles);

      // 2.b) apply Px (grad Sx)
      AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                        LeafCell->getLocal(idxRhs), TargetParticles);
    }
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








#endif //FUNIFSYMKERNELS_HPP

// [--END--]
