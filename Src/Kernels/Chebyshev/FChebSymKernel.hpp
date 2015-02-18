#ifndef FCHEBSYMKERNEL_HPP
#define FCHEBSYMKERNEL_HPP
// ===================================================================================
// Copyright ScalFmm 2011 INRIA,
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

#include "Utils/FGlobal.hpp"

#include "Utils/FSmartPointer.hpp"

#include "FAbstractChebKernel.hpp"
#include "FChebInterpolator.hpp"
#include "FChebSymM2LHandler.hpp"

class FTreeCoordinate;


// for verbosity only!!!
//#define COUNT_BLOCKED_INTERACTIONS

// if timings should be logged
//#define LOG_TIMINGS

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebSymKernel
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
template < class CellClass, class ContainerClass,   class MatrixKernelClass, int ORDER, int NVALS = 1>
class FChebSymKernel
    : public FAbstractChebKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
    typedef FAbstractChebKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS> AbstractBaseClass;
    typedef SymmetryHandler<ORDER, MatrixKernelClass::Type> SymmetryHandlerClass;
    enum {nnodes = AbstractBaseClass::nnodes};

    /// Needed for P2P and M2L operators
    const MatrixKernelClass *const MatrixKernel;

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
        assert(Loc==nullptr && Mul==nullptr && countExp==nullptr);
        Loc = new FReal* [343];
        Mul = new FReal* [343];
        countExp = new unsigned int [343];

        // set all 343 to nullptr
        for (unsigned int idx=0; idx<343; ++idx) {
            Mul[idx] = Loc[idx] = nullptr;
        }

        // init only 16 of 343 possible translations due to symmetries
        for (int i=2; i<=3; ++i)
            for (int j=0; j<=i; ++j)
                for (int k=0; k<=j; ++k) {
                    const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
                    assert(Mul[idx]==nullptr || Loc[idx]==nullptr);
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
    FChebSymKernel(const int inTreeHeight,
                   const FReal inBoxWidth,
                   const FPoint& inBoxCenter,
                   const MatrixKernelClass *const inMatrixKernel,
                   const FReal Epsilon)
      : AbstractBaseClass(inTreeHeight, inBoxWidth, inBoxCenter),
        MatrixKernel(inMatrixKernel),
        SymHandler(new SymmetryHandlerClass(MatrixKernel, Epsilon, inBoxWidth, inTreeHeight)),
        Loc(nullptr), Mul(nullptr), countExp(nullptr)
    {
        this->allocateMemoryForPermutedExpansions();
#ifdef LOG_TIMINGS
        t_m2l_1 = FReal(0.);
        t_m2l_2 = FReal(0.);
        t_m2l_3 = FReal(0.);
#endif
    }
    
    /**
     * The constructor initializes all constant attributes and it reads the
     * precomputed and compressed M2L operators from a binary file (an
     * runtime_error is thrown if the required file is not valid).
     */
    FChebSymKernel(const int inTreeHeight,
                   const FReal inBoxWidth,
                   const FPoint& inBoxCenter,
                   const MatrixKernelClass *const inMatrixKernel) :FChebSymKernel(inTreeHeight,inBoxWidth, inBoxCenter,inMatrixKernel,FMath::pow(10.0,static_cast<FReal>(-ORDER)))
    {}

    

    /** Copy constructor */
    FChebSymKernel(const FChebSymKernel& other)
        : AbstractBaseClass(other),
          MatrixKernel(other.MatrixKernel),
          SymHandler(other.SymHandler),
          Loc(nullptr), Mul(nullptr), countExp(nullptr)
    {
        this->allocateMemoryForPermutedExpansions();
    }



    /** Destructor */
    ~FChebSymKernel()
    {
        for (unsigned int t=0; t<343; ++t) {
            if (Loc[t]!=nullptr) delete [] Loc[t];
            if (Mul[t]!=nullptr) delete [] Mul[t];
        }
        if (Loc!=nullptr)      delete [] Loc;
        if (Mul!=nullptr)      delete [] Mul;
        if (countExp!=nullptr) delete [] countExp;

#ifdef LOG_TIMINGS
        std::cout << "- Permutation took " << t_m2l_1 << "s"
                  << "\n- GEMMT and GEMM took " << t_m2l_2 << "s"
                  << "\n- Unpermutation took " << t_m2l_3 << "s"
                  << std::endl;
#endif
    }


    const SymmetryHandlerClass * getPtrToSymHandler() const
    {   return SymHandler.getPtr(); }
    


    void P2M(CellClass* const LeafCell,
                     const ContainerClass* const SourceParticles/*, const int level = AbstractBaseClass::TreeHeight*/)
    {
        // apply Sy
        const FPoint LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafCell->getCoordinate()));
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                      LeafCell->getMultipole(idxRhs), SourceParticles);
        }
    }



    void M2M(CellClass* const FRestrict ParentCell,
             const CellClass*const FRestrict *const FRestrict ChildCells,
             const int /*TreeLevel*/)
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // Reset the Parent expansion to zero
            //   FBlas::scal(nnodes*2, FReal(0.), ParentCell->getMultipole(idxRhs));
            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                // apply Sy
               if (ChildCells[ChildIndex]){
                    AbstractBaseClass::Interpolator->applyM2M(ChildIndex, ChildCells[ChildIndex]->getMultipole(idxRhs), ParentCell->getMultipole(idxRhs));
                }
            }
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
                    const unsigned int pidx = SymHandler->pindices[idx];
                    const unsigned int count = (countExp[pidx])++;
                    FReal *const mul = Mul[pidx] + count*nnodes;
                    const unsigned int *const pvec = SymHandler->pvectors[idx];
                    const FReal *const MultiExp = SourceCells[idx]->getMultipole(idxRhs);

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
            const FReal scale = MatrixKernel->getScaleFactor(AbstractBaseClass::BoxWidth, TreeLevel);
            for (unsigned int pidx=0; pidx<343; ++pidx) {
                const unsigned int count = countExp[pidx];
                if (count) {
                    const unsigned int rank = SymHandler->getLowRank(TreeLevel, pidx);
                    // rank * count * (2*nnodes-1) flops
                    FBlas::gemtm(nnodes, rank, count, FReal(1.),
                                 const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx))+rank*nnodes,
                                 nnodes, Mul[pidx], nnodes, Compressed, rank);
                    // nnodes *count * (2*rank-1) flops
                    FBlas::gemm( nnodes, rank, count, scale,
                                 const_cast<FReal*>(SymHandler->getK(TreeLevel, pidx)),
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
            FReal *const LocalExpansion = TargetCell->getLocal(idxRhs);
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
        }

    }



    /*
    void M2L(CellClass* const FRestrict TargetCell,
                     const CellClass* SourceCells[343],
                     const int NumSourceCells,
                     const int TreeLevel)
    {
        // permute and copy multipole expansion
        memset(countExp, 0, sizeof(int) * 343);
        for (unsigned int idx=0; idx<343; ++idx) {
            if (SourceCells[idx]) {
                const unsigned int pidx = SymHandler->pindices[idx];
                const unsigned int count = (countExp[pidx])++;
                const FReal *const MultiExp = SourceCells[idx]->getMultipole();
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
                const unsigned int rank = SymHandler->LowRank[pidx];
                FBlas::gemtm(nnodes, rank, count, FReal(1.),
                                         SymHandler->K[pidx]+rank*nnodes, nnodes, Mul[pidx], nnodes, Compressed, rank);
                FBlas::gemm( nnodes, rank, count, scale,
                                         SymHandler->K[pidx], nnodes, Compressed, rank, Loc[pidx], nnodes);
            }
        }

        // permute and add contribution to local expansions
        FReal *const LocalExpansion = TargetCell->getLocal();
        memset(countExp, 0, sizeof(int) * 343);
        for (unsigned int idx=0; idx<343; ++idx) {
            if (SourceCells[idx]) {
                const unsigned int pidx = SymHandler->pindices[idx];
                const unsigned int count = (countExp[pidx])++;
                for (unsigned int n=0; n<nnodes; ++n)
                    LocalExpansion[n] += Loc[pidx][count*nnodes + SymHandler->pvectors[idx][n]];
            }
        }

    }
    */

    /*
        void M2L(CellClass* const FRestrict TargetCell,
        const CellClass* SourceCells[343],
        const int NumSourceCells,
        const int TreeLevel)
        {
        FReal *const LocalExpansion = TargetCell->getLocal();
        const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
        const FReal scale(MatrixKernel->getScaleFactor(CellWidth));

        FReal PermLocalExp[nnodes];
        FReal PermMultiExp[nnodes];
        FReal Compressed[nnodes];
        for (int i=-3; i<=3; ++i) {
        for (int j=-3; j<=3; ++j) {
        for (int k=-3; k<=3; ++k) {
                    
        const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
                    
        if (SourceCells[idx]) {
        const FReal *const MultiExp = SourceCells[idx]->getMultipole();

        // permute
        for (unsigned int n=0; n<nnodes; ++n) PermMultiExp[pvectors[idx][n]] = MultiExp[n];

        // mat-vec-mult (svd)
        assert(K[pindices[idx]]!=nullptr);
        const int rank = LowRank[pindices[idx]];
        FBlas::gemtv(nnodes, rank, FReal(1.), K[pindices[idx]]+rank*nnodes, PermMultiExp, Compressed);
        FBlas::gemv( nnodes, rank, scale, K[pindices[idx]], Compressed, PermLocalExp);

        // permute
        for (unsigned int n=0; n<nnodes; ++n) LocalExpansion[n] += PermLocalExp[pvectors[idx][n]];
        }

        }
        }
        }

        }
    */



    void L2L(const CellClass* const FRestrict ParentCell,
             CellClass* FRestrict *const FRestrict ChildCells,
             const int /*TreeLevel*/)
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // apply Sx
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
//      // a) apply Sx
//      AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter,
//                                                AbstractBaseClass::BoxWidthLeaf,
//                                                LeafCell->getLocal(),
//                                                TargetParticles);
//      // b) apply Px (grad Sx)
//      AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter,
//                                                        AbstractBaseClass::BoxWidthLeaf,
//                                                        LeafCell->getLocal(),
//                                                        TargetParticles);

            // c) apply Sx and Px (grad Sx)
            AbstractBaseClass::Interpolator->applyL2PTotal(LeafCellCenter, AbstractBaseClass::BoxWidthLeaf,
                                                           LeafCell->getLocal(idxRhs), TargetParticles);
        }
    }

    void P2P(const FTreeCoordinate& /* LeafCellCoordinate */, // needed for periodic boundary conditions
             ContainerClass* const FRestrict TargetParticles,
             const ContainerClass* const FRestrict /*SourceParticles*/,
             ContainerClass* const NeighborSourceParticles[27],
             const int /* size */)
    {
        DirectInteractionComputer<MatrixKernelClass::NCMP, NVALS>::P2P(TargetParticles,NeighborSourceParticles,MatrixKernel);
    }


    void P2PRemote(const FTreeCoordinate& /*inPosition*/,
                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
                   ContainerClass* const inNeighbors[27], const int /*inSize*/)
    {
        DirectInteractionComputer<MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,27,MatrixKernel);
    }

};








#endif //FCHEBSYMKERNELS_HPP

// [--END--]
