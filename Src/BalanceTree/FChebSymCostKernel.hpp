// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE


#ifndef FCHEBFLOPSSYMKERNEL_HPP
#define FCHEBFLOPSSYMKERNEL_HPP

#include <stdexcept>

#include "Utils/FGlobal.hpp"

#include "Utils/FSmartPointer.hpp"

#include "Components/FAbstractKernels.hpp"

#include "Kernels/Chebyshev/FChebInterpolator.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpSymmetries.hpp"

class FTreeCoordinate;

/**
 * \author Quentin Khan, Matthias Messner (original file: FChebFlopsSymKernel)
 * \brief Cost computation of Chebyshev interpolation based FMM.
 *
 * Please read the license
 *
 * This kernel implements the cost computation of the Chebyshev interpolation
 * based FMM operators exploiting the symmetries in the far-field. It implements
 * all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by the
 * FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * \tparam Freal Type of real number representation
 * \tparam CellClass Type of cell
 * \tparam ContainerClass Type of container to store particles
 * \tparam MatrixKernelClass Type of matrix kernel function
 * \tparam ORDER Chebyshev interpolation order
 * \tparam OctreeClass Class of the octree to work on.
 */
template < typename FReal, class CellClass, class ContainerClass, class MatrixKernelClass, int ORDER, class OctreeClass>
class FChebSymCostKernel : public FAbstractKernels<CellClass, ContainerClass>
{
    enum {nnodes = TensorTraits<ORDER>::nnodes};
public:
    // Class types available to the rest of the world
    /// Type of real number representation
    using _FReal = FReal;
    /// Type of cell
    using _CellClass = CellClass;
    /// Type of container to store particles
    using _ContainerClass = ContainerClass;
    /// Type of matrix kernel function
    using _MatrixKernelClass = MatrixKernelClass;
    /// Chebyshev interpolation order
    constexpr static const int order = ORDER; 
    /// Class of the octree to work on
    using _OctreeClass = OctreeClass;


    /// Handler to deal with all symmetries: Stores permutation indices and
    /// vectors to reduce 343 different interactions to 16 only.
    struct SymmetryHandler;

    /// Needed for handling all symmetries
    const FSmartPointer<MatrixKernelClass,FSmartPointerMemory> MatrixKernel;
    const FSmartPointer<SymmetryHandler,  FSmartPointerMemory> SymHandler;

    /// The tree the kernel is working on. Needed to attain cells in the P2P
    /// operator (we only get the particles containers otherwise)
    const OctreeClass* const _tree;

    /// The tree height
    const unsigned int _treeHeight;

    /// count permuted local and multipole expansions
    unsigned int* countExp;

    /// Flops count for each operator of the FMM. 
    unsigned long long flopsP2M = 0,
        flopsM2M = 0,
        flopsM2L = 0,
        flopsL2L = 0,
        flopsL2P = 0,
        flopsP2P = 0;

    /// Operators count.
    unsigned long long countP2M = 0,
        countM2M = 0,
        countM2L = 0,
        countL2L = 0,
        countL2P = 0,
        countP2P = 0;



    /// start flop counters 
    FSize countFlopsM2MorL2L() const
        { return 3 * nnodes * (2*ORDER-1); }

    FSize countFlopsM2L(const unsigned int nexp, const unsigned int rank) const
        { return nexp * (4*nnodes*rank - rank - nnodes); }

    FSize countFlopsP2P() const
        { return 34; }

    FSize countFlopsP2Pmutual() const
        { return 39; }

    FSize countFlopsP2M(const FSize N) const    {
        const FSize first  = N * (18 + (ORDER-2) * 6 + (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1) * 2)));
        const FSize W2 = 3 * ORDER*(2*(ORDER-1)-1);
        const FSize W4 = 3 * (ORDER*(ORDER-1)*(2*(ORDER-1)-1) + ORDER*ORDER*(2*(ORDER-1)-1));
        const FSize W8 = 3 * (2*(ORDER-1)-1) * (ORDER*(ORDER-1)*(ORDER-1) + ORDER*ORDER*(ORDER-1) + nnodes);
        return first + W2 + W4 + W8 + nnodes*11;
    }

    FSize countFlopsL2PTotal(const FSize N) const    {
        const unsigned W0 = nnodes;
        const unsigned W2 = 3 * (ORDER-1)*ORDER*ORDER * 2*ORDER;
        const unsigned W4 = 3 * ORDER*(ORDER-1)*(ORDER-1) * 2*ORDER;
        const unsigned W8 = (ORDER-1)*(ORDER-1)*(ORDER-1) * (2*ORDER-1);
        const FSize second = N * (38 + (ORDER-2)*15 + (ORDER-1)*((ORDER-1) * (27 + (ORDER-1) * 16))) + 6;
        return W0 + W2 + W4 + W8 + second;
    }
    // end flop counters

public:
    /**
     * The constructor initializes all constant attributes and it reads the
     * precomputed and compressed M2L operators from a binary file (an
     * runtime_error is thrown if the required file is not valid).
     */
    FChebSymCostKernel(OctreeClass* tree,
                          const FReal Epsilon)
        : MatrixKernel(new MatrixKernelClass()),
          SymHandler(new SymmetryHandler(MatrixKernel.getPtr(), Epsilon)),
          _tree(tree),
          _treeHeight(_tree->getHeight())
        {
            countExp = new unsigned int [343];
        }
    


    /** Copy constructor */
    FChebSymCostKernel(const FChebSymCostKernel& other) :
        SymHandler(other.SymHandler),
        _tree(other._tree),
        _treeHeight(other._treeHeight)
        {
        countExp = new unsigned int [343];
    }

    /** Destructor */
    ~FChebSymCostKernel() {
        delete [] countExp;

    }
    
    void printResults(std::ostream& os) const {
        os << "\n==================================================" 
           << "\n- P2M Flops:" << flopsP2M 
           << "\n- M2M Flops:" << flopsM2M
           << "\n- M2L Flops:" << flopsM2L
           << "\n- L2L Flops:" << flopsL2L
           << "\n- L2P Flops:" << flopsL2P
           << "\n- P2P Flops:" << flopsP2P
           << "\n- Overall Flops = " << flopsP2M + flopsM2M + flopsM2L
            + flopsL2L + flopsL2P + flopsP2P
           << "\n==================================================\n"
           << std::endl;

        os << "P2P count: " << countP2P << std::endl;
        os << "P2M count: " << countP2M << std::endl;
        os << "M2M count: " << countM2M << std::endl;
        os << "M2L count: " << countM2L << std::endl;
        os << "L2L count: " << countL2L << std::endl;
        os << "L2P count: " << countL2P << std::endl;
        
    }
    
    
    void P2M(CellClass* const cell, const ContainerClass* const SourceParticles) override {
        FSize tmpCost = countFlopsP2M(SourceParticles->getNbParticles());
        flopsP2M += tmpCost;
        cell->addCost(tmpCost);
        countP2M++;
    }



    void M2M(CellClass* const FRestrict cell,
             const CellClass*const FRestrict *const FRestrict ChildCells,
             const int /*TreeLevel*/) override {
        FSize flops = 0;
        for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
            if (ChildCells[ChildIndex])    flops += countFlopsM2MorL2L();
        flopsM2M += flops;
        cell->addCost(flops);
        countM2M++;
    }




    void M2L(CellClass* const FRestrict cell,
             const CellClass* /*SourceCells*/[],
             const int positions[],
             const int size,
             const int /* TreeLevel */) override
        {
            FSize flops = 0;
            // count how ofter each of the 16 interactions is used
            memset(countExp, 0, sizeof(int) * 343);
            for (int idx=0; idx<size; ++idx)
                  countExp[SymHandler->pindices[positions[idx]]]++;
            // multiply (mat-mat-mul)
            for (int pidx=0; pidx<343; ++pidx)
                if (countExp[pidx])
                    flops += countFlopsM2L(countExp[pidx], SymHandler->LowRank[pidx])
                        + countExp[pidx]*nnodes;
            flopsM2L += flops;
            cell->addCost(flops);
            countM2L++;
        }


    void L2L(const CellClass* const FRestrict /* not needed */,
             CellClass* FRestrict *const FRestrict ChildCells,
             const int /* TreeLevel*/) override {
        FSize flops = 0;
        FSize tmpCost = 0;
        for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex)
            if (ChildCells[ChildIndex]) {
                tmpCost = countFlopsM2MorL2L() + nnodes;
                flops += tmpCost;
                ChildCells[ChildIndex]->addCost(tmpCost);
            }

        flopsL2L += flops;
        countL2L++;
    }



    void L2P(const CellClass* const cell,
             ContainerClass* const TargetParticles) override {
        //// 1.a) apply Sx
        //flopsL2P += countFlopsP2MorL2P(TargetParticlesParticles->getNbParticles()) + TargetParticles->getNbParticles();
        //// 1.b) apply Px (grad Sx)
        //flopsL2P += countFlopsL2PGradient(TargetParticlesParticles->getNbParticles()) + 3 * TargetParticles->getNbParticles();
        
        // or
        
        // 2) apply Sx and Px (grad Sx)
        FSize tmpCost = 0;
        tmpCost =  countFlopsL2PTotal(TargetParticles->getNbParticles()) + 4 * TargetParticles->getNbParticles();
        flopsL2P += tmpCost;
        cell->addCost(tmpCost);
        countL2P++;
    }

    

    void P2P(const FTreeCoordinate& LeafCellCoordinate, // needed for periodic boundary conditions
             ContainerClass* const FRestrict TargetParticles,
             const ContainerClass* const FRestrict SourceParticles,
             ContainerClass* const NeighborSourceParticles[],
             const int positions[],
             const int size) override {
        FSize tmpCost = 0;
        FSize srcPartCount = SourceParticles->getNbParticles();
        FSize tgtPartCount = TargetParticles->getNbParticles();

        if ( TargetParticles != SourceParticles ) {
            tmpCost +=  countFlopsP2P() * tgtPartCount * srcPartCount;

            for ( int idx = 0; idx < size; ++idx ) {
                    tmpCost += 
                        countFlopsP2P() 
                        * tgtPartCount
                        * NeighborSourceParticles[idx]->getNbParticles();
            }
        } else {
            tmpCost +=
                countFlopsP2Pmutual()
                * ((tgtPartCount * tgtPartCount
                    + tgtPartCount)
                   / 2);
            for (int idx=0; idx < size && positions[idx]<=13; ++idx)
            {
                    tmpCost += 
                        countFlopsP2Pmutual()
                        * tgtPartCount
                        * NeighborSourceParticles[idx]->getNbParticles();                    
                }
        }

        flopsP2P += tmpCost;

        CellClass* cell = _tree->getCell(
            LeafCellCoordinate.getMortonIndex(),
            _treeHeight - 1);

        cell->addNearCost(tmpCost);
        countP2P++;
    }

    void P2POuter(const FTreeCoordinate& LeafCellCoordinate, // needed for periodic boundary conditions
             ContainerClass* const FRestrict TargetParticles,
             ContainerClass* const NeighborSourceParticles[],
             const int positions[],
             const int size) override {
        FSize tmpCost = 0;
        FSize tgtPartCount = TargetParticles->getNbParticles();

        {
            for (int idx=0; idx < size && positions[idx]<=13; ++idx)
            {
                    tmpCost +=
                        countFlopsP2Pmutual()
                        * tgtPartCount
                        * NeighborSourceParticles[idx]->getNbParticles();
                }
        }

        flopsP2P += tmpCost;

        CellClass* cell = _tree->getCell(
            LeafCellCoordinate.getMortonIndex(),
            _treeHeight - 1);

        cell->addNearCost(tmpCost);
        countP2P++;
    }
};



/**
 * Handler to deal with all symmetries: Stores permutation indices and vectors
 * to reduce 343 different interactions to 16 only.
 */
template < typename FReal, class CellClass, class ContainerClass,
           class MatrixKernelClass, int ORDER, class OctreeClass>
struct FChebSymCostKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, OctreeClass>
::SymmetryHandler
{
    // M2L operators
    FReal*    K[343];
    int LowRank[343];
        
    // permutation vectors and permutated indices
    unsigned int pvectors[343][nnodes];
    unsigned int pindices[343];


    // compute rank
    unsigned int getRank(const FReal singular_values[], const double eps)
        {
            FReal nrm2(0.);
            for (unsigned int k=0; k<nnodes; ++k)
                nrm2 += singular_values[k] * singular_values[k];
        
            FReal nrm2k(0.);
            for (unsigned int k=nnodes; k>0; --k) {
                nrm2k += singular_values[k-1] * singular_values[k-1];
                if (nrm2k > eps*eps * nrm2)    return k;
            }
            throw std::runtime_error("rank cannot be larger than nnodes");
            return 0;
        }
    

    /** Constructor */
    SymmetryHandler(const MatrixKernelClass *const MatrixKernel, const double Epsilon)
        {
            // init all 343 item to zero, because effectively only 16 exist
            for (unsigned int t=0; t<343; ++t) {
                K[t] = nullptr;
                LowRank[t] = 0;
            }
            
            // set permutation vector and indices
            const FInterpSymmetries<ORDER> Symmetries;
            for (int i=-3; i<=3; ++i)
                for (int j=-3; j<=3; ++j)
                    for (int k=-3; k<=3; ++k)
                        if (abs(i)>1 || abs(j)>1 || abs(k)>1) {
                            const unsigned int idx = ((i+3) * 7 + (j+3)) * 7 + (k+3);
                            pindices[idx] = Symmetries.getPermutationArrayAndIndex(i,j,k, pvectors[idx]);
                        }

            // precompute 16 M2L operators
            this->precomputeSVD(MatrixKernel, Epsilon);
        }



    /** Destructor */
    ~SymmetryHandler()
        {
            for (unsigned int t=0; t<343; ++t)
                if (K[  t]!=nullptr) delete [] K[  t];
        }



private:
    void precomputeSVD(const MatrixKernelClass *const MatrixKernel, const double Epsilon)
        {
            // interpolation points of source (Y) and target (X) cell
            FPoint<FReal> X[nnodes], Y[nnodes];
            // set roots of target cell (X)
            FChebTensor<FReal, ORDER>::setRoots(FPoint<FReal>(0.,0.,0.), FReal(2.), X);
            // temporary matrix
            FReal* U = new FReal [nnodes*nnodes];

            // needed for the SVD
            const unsigned int LWORK = 2 * (3*nnodes + nnodes);
            FReal *const WORK = new FReal [LWORK];
            FReal *const VT = new FReal [nnodes*nnodes];
            FReal *const S = new FReal [nnodes];
        
            unsigned int counter = 0;
            for (int i=2; i<=3; ++i) {
                for (int j=0; j<=i; ++j) {
                    for (int k=0; k<=j; ++k) {

                        // assemble matrix
                        const FPoint<FReal> cy(FReal(2.*i), FReal(2.*j), FReal(2.*k));
                        FChebTensor<FReal, ORDER>::setRoots(cy, FReal(2.), Y);
                        for (unsigned int n=0; n<nnodes; ++n)
                            for (unsigned int m=0; m<nnodes; ++m)
                                U[n*nnodes + m] = MatrixKernel->evaluate(X[m], Y[n]);

                        // applying weights ////////////////////////////////////////
                        FReal weights[nnodes];
                        FChebTensor<FReal,ORDER>::setRootOfWeights(weights);
                        for (unsigned int n=0; n<nnodes; ++n) {
                            FBlas::scal(nnodes, weights[n], U + n,  nnodes); // scale rows
                            FBlas::scal(nnodes, weights[n], U + n * nnodes); // scale cols
                        }
                        //////////////////////////////////////////////////////////        

                        // truncated singular value decomposition of matrix
                        const unsigned int info    = FBlas::gesvd(nnodes, nnodes, U, S, VT, nnodes, LWORK, WORK);
                        if (info!=0) throw std::runtime_error("SVD did not converge with " + std::to_string(info));
                        const unsigned int rank = this->getRank(S, Epsilon);

                        // store 
                        const unsigned int idx = (i+3)*7*7 + (j+3)*7 + (k+3);
                        assert(K[idx]==nullptr);
                        K[idx] = new FReal [2*rank*nnodes];
                        LowRank[idx] = rank;
                        for (unsigned int r=0; r<rank; ++r)
                            FBlas::scal(nnodes, S[r], U + r*nnodes);
                        FBlas::copy(rank*nnodes, U,  K[idx]);
                        for (unsigned int r=0; r<rank; ++r)
                            FBlas::copy(nnodes, VT + r, nnodes, K[idx] + rank*nnodes + r*nnodes, 1);

                        // un-weighting ////////////////////////////////////////////
                        for (unsigned int n=0; n<nnodes; ++n) {
                            // scale rows
                            FBlas::scal(rank, FReal(1.) / weights[n], K[idx] + n,               nnodes);
                            // scale rows
                            FBlas::scal(rank, FReal(1.) / weights[n], K[idx] + rank*nnodes + n, nnodes);
                        }
                        //////////////////////////////////////////////////////////        

                        //std::cout << "(" << i << "," << j << "," << k << ") " << idx <<
                        //    ", low rank = " << rank << std::endl;

                        counter++;
                    }
                }
            }
            std::cout << "num interactions = " << counter << std::endl;
            delete [] U;
            delete [] WORK;
            delete [] VT;
            delete [] S;
        }

};


#endif 

// [--END--]
