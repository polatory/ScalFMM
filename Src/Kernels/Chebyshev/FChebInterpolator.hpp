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
#ifndef FCHEBINTERPOLATOR_HPP
#define FCHEBINTERPOLATOR_HPP


#include "../Interpolation/FInterpMapping.hpp"
#include "../Interpolation/FInterpMatrixKernel.hpp" //PB
#include "FChebTensor.hpp"
#include "FChebRoots.hpp"

#include "Utils/FBlas.hpp"



/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */

/**
 * @class FChebInterpolator
 *
 * The class @p FChebInterpolator defines the anterpolation (M2M) and
 * interpolation (L2L) concerning operations.
 * 
 * PB: MatrixKernelClass is passed as template in order to inform interpolator 
 * of the size of the vectorial interpolators. Default is the scalar 
 * matrix kernel class of type ONE_OVER_R (NRHS=NLHS=1).
 */
template <int ORDER, class MatrixKernelClass = struct FInterpMatrixKernelR>
class FChebInterpolator : FNoCopyable
{
  // compile time constants and types
  enum {nnodes = TensorTraits<ORDER>::nnodes,
        nRhs = MatrixKernelClass::NRHS,
        nLhs = MatrixKernelClass::NLHS,
        nPV = MatrixKernelClass::NPV};
  typedef FChebRoots< ORDER>  BasisType;
  typedef FChebTensor<ORDER> TensorType;

protected: // PB for OptiDis

  FReal T_of_roots[ORDER][ORDER];
  FReal T[ORDER * (ORDER-1)];
  unsigned int node_ids[nnodes][3];

  // 8 Non-leaf (i.e. M2M/L2L) interpolators 
  // x1 per level if box is extended
  // only 1 is required for all levels if extension is 0
  FReal*** ChildParentInterpolator;

  // Tree height (needed by M2M/L2L if cell width is extended)
  const int TreeHeight;
  // Root cell width (only used by M2M/L2L)
  const FReal RootCellWidth;
  // Cell width extension (only used by M2M/L2L, kernel handles extension for P2M/L2P)
  const FReal CellWidthExtension;


    // permutations (only needed in the tensor product interpolation case)
    unsigned int perm[3][nnodes];

    ////////////////////////////////////////////////////////////////////
    // needed for P2M
    struct IMN2MNI {
        enum {size = ORDER * (ORDER-1) * (ORDER-1)};
        unsigned int imn[size], mni[size];
        IMN2MNI() {
            unsigned int counter = 0;
            for (unsigned int i=0; i<ORDER; ++i) {
                for (unsigned int m=0; m<ORDER-1; ++m) {
                    for (unsigned int n=0; n<ORDER-1; ++n) {
                        imn[counter] = n*(ORDER-1)*ORDER + m*ORDER + i;
                        mni[counter] = i*(ORDER-1)*(ORDER-1) + n*(ORDER-1) + m;
                        counter++;
                    }
                }
            }
        }
    } perm0;

    struct JNI2NIJ {
        enum {size = ORDER * ORDER * (ORDER-1)};
        unsigned int jni[size], nij[size];
        JNI2NIJ() {
            unsigned int counter = 0;
            for (unsigned int i=0; i<ORDER; ++i) {
                for (unsigned int j=0; j<ORDER; ++j) {
                    for (unsigned int n=0; n<ORDER-1; ++n) {
                        jni[counter] = i*(ORDER-1)*ORDER + n*ORDER + j;
                        nij[counter] = j*ORDER*(ORDER-1) + i*(ORDER-1) + n;
                        counter++;
                    }
                }
            }
        }
    } perm1;

    struct KIJ2IJK {
        enum {size = ORDER * ORDER * ORDER};
        unsigned int kij[size], ijk[size];
        KIJ2IJK() {
            unsigned int counter = 0;
            for (unsigned int i=0; i<ORDER; ++i) {
                for (unsigned int j=0; j<ORDER; ++j) {
                    for (unsigned int k=0; k<ORDER; ++k) {
                        kij[counter] = j*ORDER*ORDER + i*ORDER + k;
                        ijk[counter] = k*ORDER*ORDER + j*ORDER + i;
                        counter++;
                    }
                }
            }
        }
    } perm2;
    ////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////
    // needed for L2P
    struct IJK2JKI {
        enum {size = ORDER * ORDER * ORDER};
        unsigned int ijk[size], jki[size];
        IJK2JKI() {
            unsigned int counter = 0;
            for (unsigned int i=0; i<ORDER; ++i) {
                for (unsigned int j=0; j<ORDER; ++j) {
                    for (unsigned int k=0; k<ORDER; ++k) {
                        ijk[counter] = k*ORDER*ORDER + j*ORDER + i;
                        jki[counter] = i*ORDER*ORDER + k*ORDER + j;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[jki[i]] = IN[ijk[i]]; }
    } perm3;

    struct IJK2KIJ {
        enum {size = ORDER * ORDER * ORDER};
        unsigned int ijk[size], kij[size];
        IJK2KIJ() {
            unsigned int counter = 0;
            for (unsigned int i=0; i<ORDER; ++i) {
                for (unsigned int j=0; j<ORDER; ++j) {
                    for (unsigned int k=0; k<ORDER; ++k) {
                        ijk[counter] = k*ORDER*ORDER + j*ORDER + i;
                        kij[counter] = j*ORDER*ORDER + i*ORDER + k;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[kij[i]] = IN[ijk[i]]; }
    } perm4;

    struct LJK2JKL {
        enum {size = (ORDER-1) * ORDER * ORDER};
        unsigned int ljk[size], jkl[size];
        LJK2JKL() {
            unsigned int counter = 0;
            for (unsigned int l=0; l<ORDER-1; ++l) {
                for (unsigned int j=0; j<ORDER; ++j) {
                    for (unsigned int k=0; k<ORDER; ++k) {
                        ljk[counter] = k*ORDER*(ORDER-1) + j*(ORDER-1) + l;
                        jkl[counter] = l*ORDER*ORDER + k*ORDER + j;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[jkl[i]] = IN[ljk[i]]; }
    } perm5;

    struct LJK2KLJ {
        enum {size = (ORDER-1) * ORDER * ORDER};
        unsigned int ljk[size], klj[size];
        LJK2KLJ() {
            unsigned int counter = 0;
            for (unsigned int l=0; l<ORDER-1; ++l) {
                for (unsigned int j=0; j<ORDER; ++j) {
                    for (unsigned int k=0; k<ORDER; ++k) {
                        ljk[counter] = k*ORDER*(ORDER-1) + j*(ORDER-1) + l;
                        klj[counter] = j*(ORDER-1)*ORDER + l*ORDER + k;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[klj[i]] = IN[ljk[i]]; }
    } perm6;

    struct MKI2KIM {
        enum {size = (ORDER-1) * ORDER * ORDER};
        unsigned int mki[size], kim[size];
        MKI2KIM() {
            unsigned int counter = 0;
            for (unsigned int m=0; m<ORDER-1; ++m) {
                for (unsigned int k=0; k<ORDER; ++k) {
                    for (unsigned int i=0; i<ORDER; ++i) {
                        mki[counter] = i*ORDER*(ORDER-1) + k*(ORDER-1) + m;
                        kim[counter] = m*ORDER*ORDER + i*ORDER + k;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[kim[i]] = IN[mki[i]]; }
    } perm7;

    struct MKL2KLM {
        enum {size = (ORDER-1) * ORDER * (ORDER-1)};
        unsigned int mkl[size], klm[size];
        MKL2KLM() {
            unsigned int counter = 0;
            for (unsigned int m=0; m<ORDER-1; ++m) {
                for (unsigned int k=0; k<ORDER; ++k) {
                    for (unsigned int l=0; l<ORDER-1; ++l) {
                        mkl[counter] = l*ORDER*(ORDER-1) + k*(ORDER-1) + m;
                        klm[counter] = m*(ORDER-1)*ORDER + l*ORDER + k;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[klm[i]] = IN[mkl[i]]; }
    } perm8;

    struct NLM2LMN {
        enum {size = (ORDER-1) * (ORDER-1) * (ORDER-1)};
        unsigned int nlm[size], lmn[size];
        NLM2LMN() {
            unsigned int counter = 0;
            for (unsigned int n=0; n<ORDER-1; ++n) {
                for (unsigned int l=0; l<ORDER-1; ++l) {
                    for (unsigned int m=0; m<ORDER-1; ++m) {
                        nlm[counter] = m*(ORDER-1)*(ORDER-1) + l*(ORDER-1) + n;
                        lmn[counter] = n*(ORDER-1)*(ORDER-1) + m*(ORDER-1) + l;
                        counter++;
                    }
                }
            }
        }
        void permute(const FReal *const IN, FReal *const OUT) const
        { for (unsigned int i=0; i<size; ++i) OUT[lmn[i]] = IN[nlm[i]]; }
    } perm9;

    ////////////////////////////////////////////////////////////////////



    /**
     * Initialize the child - parent - interpolator, it is basically the matrix
     * S which is precomputed and reused for all M2M and L2L operations, ie for
     * all non leaf inter/anterpolations.
     * This is a sub-optimal version : complexity p^6.
     * This function handles extended cells.
     */
  void initM2MandL2L(const int TreeLevel, const FReal ParentWidth)
  {
    FPoint ChildRoots[nnodes];

    // Ratio of extended cell widths (definition: child ext / parent ext)
    const FReal ExtendedCellRatio = 
      FReal(FReal(ParentWidth)/FReal(2.) + CellWidthExtension) / FReal(ParentWidth + CellWidthExtension);

    // Child cell center and width
    FPoint ChildCenter;
    const FReal ChildWidth(2.*ExtendedCellRatio);

    // loop: child cells
    for (unsigned int child=0; child<8; ++child) {

      // allocate memory
      ChildParentInterpolator[TreeLevel][child] = new FReal [nnodes * nnodes];

      // set child info
      FChebTensor<ORDER>::setRelativeChildCenter(child, ChildCenter, ExtendedCellRatio);
      FChebTensor<ORDER>::setRoots(ChildCenter, ChildWidth, ChildRoots);

      // assemble child - parent - interpolator
      assembleInterpolator(nnodes, ChildRoots, ChildParentInterpolator[TreeLevel][child]);
    }
  }

  
  /**
   * Initialize the child - parent - interpolator, it is basically the matrix
   * S which is precomputed and reused for all M2M and L2L operations, ie for
   * all non leaf inter/anterpolations.
   * This is a more optimal version : complexity p^4.
   * This function handles extended cells.
   *
   */
  void initTensorM2MandL2L(const int TreeLevel, const FReal ParentWidth)
  {
    FReal ChildCoords[3][ORDER]; 
    FPoint ChildCenter;
       
    // Ratio of extended cell widths (definition: child ext / parent ext)
    const FReal ExtendedCellRatio = 
      FReal(FReal(ParentWidth)/FReal(2.) + CellWidthExtension) / FReal(ParentWidth + CellWidthExtension);

    // Child cell width
    const FReal ChildWidth(2.*ExtendedCellRatio);

    // loop: child cells
    for (unsigned int child=0; child<8; ++child) {

      // set child info 
      FChebTensor<ORDER>::setRelativeChildCenter(child, ChildCenter, ExtendedCellRatio);
      FChebTensor<ORDER>::setPolynomialsRoots(ChildCenter, ChildWidth, ChildCoords);
      // allocate memory
      ChildParentInterpolator[TreeLevel][child] = new FReal [3 * ORDER*ORDER];
      assembleInterpolator(ORDER, ChildCoords[0], ChildParentInterpolator[TreeLevel][child]);
      assembleInterpolator(ORDER, ChildCoords[1], ChildParentInterpolator[TreeLevel][child] + 1 * ORDER*ORDER);
      assembleInterpolator(ORDER, ChildCoords[2], ChildParentInterpolator[TreeLevel][child] + 2 * ORDER*ORDER);
    }


    // init permutations
    for (unsigned int i=0; i<ORDER; ++i) {
      for (unsigned int j=0; j<ORDER; ++j) {
        for (unsigned int k=0; k<ORDER; ++k) {
          const unsigned int index = k*ORDER*ORDER + j*ORDER + i;
          perm[0][index] = k*ORDER*ORDER + j*ORDER + i;
          perm[1][index] = i*ORDER*ORDER + k*ORDER + j;
          perm[2][index] = j*ORDER*ORDER + i*ORDER + k;
        }
      }
    }

  }



public:
    /**
     * Constructor: Initialize the Chebyshev polynomials at the Chebyshev
     * roots/interpolation point
     *
     * PB: Input parameters ONLY affect the computation of the M2M/L2L ops.
     * These parameters are ONLY required in the context of extended bbox.
     * If no M2M/L2L is required then the interpolator can be built with 
     * the default ctor.
     */
  explicit FChebInterpolator(const int inTreeHeight=3,
                             const FReal inRootCellWidth=FReal(1.), 
                             const FReal inCellWidthExtension=FReal(0.)) 
  : TreeHeight(inTreeHeight), 
    RootCellWidth(inRootCellWidth),
    CellWidthExtension(inCellWidthExtension)
    {
        // initialize chebyshev polynomials of root nodes: T_o(x_j)
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        T_of_roots[o][j] = FReal(BasisType::T(o, FReal(BasisType::roots[j])));

        // initialize chebyshev polynomials of root nodes: T_o(x_j)
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        T[(o-1)*ORDER + j] = FReal(BasisType::T(o, FReal(BasisType::roots[j])));


        // initialize root node ids
        TensorType::setNodeIds(node_ids);

        // initialize interpolation operator for M2M/L2L (non leaf operations)

        // allocate 8 arrays per level
        ChildParentInterpolator = new FReal**[TreeHeight];
        for (unsigned int l=0; l<static_cast<unsigned int>(TreeHeight); ++l){
          ChildParentInterpolator[l] = new FReal*[8];
          for (unsigned int c=0; c<8; ++c)
            ChildParentInterpolator[l][c]=nullptr;        
        }

        // Set number of non-leaf ios that actually need to be computed
        unsigned int reducedTreeHeight; // = 2 + nb of computed nl ios
        if(CellWidthExtension==0.) // if no cell extension, then ...
          reducedTreeHeight = 3; // cmp only 1 non-leaf io
        else 
          reducedTreeHeight = TreeHeight; // cmp 1 non-leaf io per level

        // Init non-leaf interpolators
        FReal CellWidth = RootCellWidth / FReal(2.); // at level 1
        CellWidth /= FReal(2.);                      // at level 2
        
        for (unsigned int l=2; l<reducedTreeHeight; ++l) {
          
          //this -> initM2MandL2L(l,CellWidth);     // non tensor-product interpolation
          this -> initTensorM2MandL2L(l,CellWidth); // tensor-product interpolation

          // update cell width
          CellWidth /= FReal(2.);                    // at level l+1 
        }
    }


    /**
     * Destructor: Delete dynamically allocated memory for M2M and L2L operator
     */
    ~FChebInterpolator()
    {
      for (unsigned int l=0; l<static_cast<unsigned int>(TreeHeight); ++l)
        for (unsigned int child=0; child<8; ++child)
          if(ChildParentInterpolator[l][child] != nullptr)
            delete [] ChildParentInterpolator[l][child];
    }


    /**
     * Assembles the interpolator \f$S_\ell\f$ of size \f$N\times
     * \ell^3\f$. Here local points is meant as points whose global coordinates
     * have already been mapped to the reference interval [-1,1].
     *
     * @param[in] NumberOfLocalPoints
     * @param[in] LocalPoints
     * @param[out] Interpolator
     */
    void assembleInterpolator(const unsigned int NumberOfLocalPoints,
                  const FPoint *const LocalPoints,
                  FReal *const Interpolator) const
    {
        // values of chebyshev polynomials of source particle: T_o(x_i)
        FReal T_of_x[ORDER][3];
        // loop: local points (mapped in [-1,1])
        for (unsigned int m=0; m<NumberOfLocalPoints; ++m) {
            // evaluate chebyshev polynomials at local points
            for (unsigned int o=1; o<ORDER; ++o) {
                T_of_x[o][0] = BasisType::T(o, LocalPoints[m].getX());
                T_of_x[o][1] = BasisType::T(o, LocalPoints[m].getY());
                T_of_x[o][2] = BasisType::T(o, LocalPoints[m].getZ());
            }

            // assemble interpolator
            for (unsigned int n=0; n<nnodes; ++n) {
                //Interpolator[n*nnodes + m] = FReal(1.);
                Interpolator[n*NumberOfLocalPoints + m] = FReal(1.);
                for (unsigned int d=0; d<3; ++d) {
                    const unsigned int j = node_ids[n][d];
                    FReal S_d = FReal(1.) / ORDER;
                    for (unsigned int o=1; o<ORDER; ++o)
                        S_d += FReal(2.) / ORDER * T_of_x[o][d] * T_of_roots[o][j];
                    //Interpolator[n*nnodes + m] *= S_d;
                    Interpolator[n*NumberOfLocalPoints + m] *= S_d;
                }

            }

        }

    }


    void assembleInterpolator(const unsigned int M, const FReal *const x, FReal *const S) const
    {
        // values of chebyshev polynomials of source particle: T_o(x_i)
        FReal T_of_x[ORDER];

        // loop: local points (mapped in [-1,1])
        for (unsigned int m=0; m<M; ++m) {
            // evaluate chebyshev polynomials at local points
            for (unsigned int o=1; o<ORDER; ++o)
                T_of_x[o] = BasisType::T(o, x[m]);

            for (unsigned int n=0; n<ORDER; ++n) {
                S[n*M + m] = FReal(1.) / ORDER;
                for (unsigned int o=1; o<ORDER; ++o)
                    S[n*M + m] += FReal(2.) / ORDER * T_of_x[o] * T_of_roots[o][n];
            }

        }

    }



    const unsigned int * getPermutationsM2ML2L(unsigned int i) const
    { return perm[i]; }






    /**
     * Particle to moment: application of \f$S_\ell(y,\bar y_n)\f$
     * (anterpolation, it is the transposed interpolation)
     */
    template <class ContainerClass>
    void applyP2M(const FPoint& center,
                                const FReal width,
                                FReal *const multipoleExpansion,
                                const ContainerClass *const sourceParticles) const;



    /**
     * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ (interpolation)
     */
    template <class ContainerClass>
    void applyL2P(const FPoint& center,
                                const FReal width,
                                const FReal *const localExpansion,
                                ContainerClass *const localParticles) const;


    /**
     * Local to particle operation: application of \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
     */
    template <class ContainerClass>
    void applyL2PGradient(const FPoint& center,
                                                const FReal width,
                                                const FReal *const localExpansion,
                                                ContainerClass *const localParticles) const;

    /**
     * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ and
     * \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
     */
    template <class ContainerClass>
    void applyL2PTotal(const FPoint& center,
                                         const FReal width,
                                         const FReal *const localExpansion,
                                         ContainerClass *const localParticles) const;


    /*
    void applyM2M(const unsigned int ChildIndex,
                                const FReal *const ChildExpansion,
                                FReal *const ParentExpansion) const
    {
        FBlas::gemtva(nnodes, nnodes, FReal(1.),
                                    ChildParentInterpolator[ChildIndex],
                                    const_cast<FReal*>(ChildExpansion), ParentExpansion);
    }

    void applyL2L(const unsigned int ChildIndex,
                                const FReal *const ParentExpansion,
                                FReal *const ChildExpansion) const
    {
        FBlas::gemva(nnodes, nnodes, FReal(1.),
                                 ChildParentInterpolator[ChildIndex],
                                 const_cast<FReal*>(ParentExpansion), ChildExpansion);
    }
    */


    void applyM2M(const unsigned int ChildIndex,
                  const FReal *const ChildExpansion,
                  FReal *const ParentExpansion,
                  const unsigned int TreeLevel = 2) const
    {
        FReal Exp[nnodes], PermExp[nnodes];
        // ORDER*ORDER*ORDER * (2*ORDER-1)
        FBlas::gemtm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                                 ChildParentInterpolator[TreeLevel][ChildIndex], ORDER,
                                 const_cast<FReal*>(ChildExpansion), ORDER, PermExp, ORDER);

        for (unsigned int n=0; n<nnodes; ++n)	Exp[n] = PermExp[perm[1][n]];
        // ORDER*ORDER*ORDER * (2*ORDER-1)
        FBlas::gemtm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                                 ChildParentInterpolator[TreeLevel][ChildIndex] + 2 * ORDER*ORDER, ORDER,
                                 Exp, ORDER, PermExp, ORDER);

        for (unsigned int n=0; n<nnodes; ++n)	Exp[perm[1][n]] = PermExp[perm[2][n]];
        // ORDER*ORDER*ORDER * (2*ORDER-1)
        FBlas::gemtm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                                 ChildParentInterpolator[TreeLevel][ChildIndex] + 1 * ORDER*ORDER, ORDER,
                                 Exp, ORDER, PermExp, ORDER);

        for (unsigned int n=0; n<nnodes; ++n)	ParentExpansion[perm[2][n]] += PermExp[n];
    }


    void applyL2L(const unsigned int ChildIndex,
                  const FReal *const ParentExpansion,
                  FReal *const ChildExpansion,
                  const unsigned int TreeLevel = 2) const
    {
        FReal Exp[nnodes], PermExp[nnodes];
        // ORDER*ORDER*ORDER * (2*ORDER-1)
        FBlas::gemm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                                ChildParentInterpolator[TreeLevel][ChildIndex], ORDER,
                                const_cast<FReal*>(ParentExpansion), ORDER, PermExp, ORDER);

        for (unsigned int n=0; n<nnodes; ++n)	Exp[n] = PermExp[perm[1][n]];
        // ORDER*ORDER*ORDER * (2*ORDER-1)
        FBlas::gemm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                                ChildParentInterpolator[TreeLevel][ChildIndex] + 2 * ORDER*ORDER, ORDER,
                                Exp, ORDER, PermExp, ORDER);

        for (unsigned int n=0; n<nnodes; ++n)	Exp[perm[1][n]] = PermExp[perm[2][n]];
        // ORDER*ORDER*ORDER * (2*ORDER-1)
        FBlas::gemm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                                ChildParentInterpolator[TreeLevel][ChildIndex] + 1 * ORDER*ORDER, ORDER,
                                Exp, ORDER, PermExp, ORDER);

        for (unsigned int n=0; n<nnodes; ++n)	ChildExpansion[perm[2][n]] += PermExp[n];
    }
    // total flops count: 3 * ORDER*ORDER*ORDER * (2*ORDER-1)
};







/**
 * Particle to moment: application of \f$S_\ell(y,\bar y_n)\f$
 * (anterpolation, it is the transposed interpolation)
 */
template <int ORDER, class MatrixKernelClass>
template <class ContainerClass>
inline void FChebInterpolator<ORDER,MatrixKernelClass>::applyP2M(const FPoint& center,
                                                                 const FReal width,
                                                                 FReal *const multipoleExpansion,
                                                                 const ContainerClass *const inParticles) const
{

    // allocate stuff
    const map_glob_loc map(center, width);
    FPoint localPosition;

    // Init W
    FReal W1[nRhs];
    FReal W2[nRhs][3][ ORDER-1];
    FReal W4[nRhs][3][(ORDER-1)*(ORDER-1)];
    FReal W8[nRhs][   (ORDER-1)*(ORDER-1)*(ORDER-1)];
    for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){
      W1[idxRhs] = FReal(0.);
      for(unsigned int i=0; i<(ORDER-1); ++i) W2[idxRhs][0][i] = W2[idxRhs][1][i] = W2[idxRhs][2][i] = FReal(0.);
      for(unsigned int i=0; i<(ORDER-1)*(ORDER-1); ++i)	W4[idxRhs][0][i] = W4[idxRhs][1][i] = W4[idxRhs][2][i] = FReal(0.);
      for(unsigned int i=0; i<(ORDER-1)*(ORDER-1)*(ORDER-1); ++i)	W8[idxRhs][i] = FReal(0.);
    }


    // loop over source particles
//    const FReal*const physicalValues = inParticles->getPhysicalValues();
    const FReal*const positionsX = inParticles->getPositions()[0];
    const FReal*const positionsY = inParticles->getPositions()[1];
    const FReal*const positionsZ = inParticles->getPositions()[2];

    for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++idxPart){
        // map global position to [-1,1]
        map(FPoint(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]), localPosition); // 15 flops

        FReal T_of_x[3][ORDER];
        T_of_x[0][0] = FReal(1.); T_of_x[0][1] = localPosition.getX();
        T_of_x[1][0] = FReal(1.); T_of_x[1][1] = localPosition.getY();
        T_of_x[2][0] = FReal(1.); T_of_x[2][1] = localPosition.getZ();
        const FReal x2 = FReal(2.) * T_of_x[0][1]; // 1 flop
        const FReal y2 = FReal(2.) * T_of_x[1][1]; // 1 flop
        const FReal z2 = FReal(2.) * T_of_x[2][1]; // 1 flop
        for (unsigned int j=2; j<ORDER; ++j) {
            T_of_x[0][j] = x2 * T_of_x[0][j-1] - T_of_x[0][j-2]; // 2 flops
            T_of_x[1][j] = y2 * T_of_x[1][j-1] - T_of_x[1][j-2]; // 2 flops
            T_of_x[2][j] = z2 * T_of_x[2][j-1] - T_of_x[2][j-2]; // 2 flops
        }
        for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){
          const FReal*const physicalValues = inParticles->getPhysicalValues(idxRhs);

          const FReal weight = physicalValues[idxPart];
          W1[idxRhs] += weight; // 1 flop
          for (unsigned int i=1; i<ORDER; ++i) {
            const FReal wx = weight * T_of_x[0][i]; // 1 flop
            const FReal wy = weight * T_of_x[1][i]; // 1 flop
            const FReal wz = weight * T_of_x[2][i]; // 1 flop
            W2[idxRhs][0][i-1] += wx; // 1 flop
            W2[idxRhs][1][i-1] += wy; // 1 flop
            W2[idxRhs][2][i-1] += wz; // 1 flop
            for (unsigned int j=1; j<ORDER; ++j) {
              const FReal wxy = wx * T_of_x[1][j]; // 1 flop
              const FReal wxz = wx * T_of_x[2][j]; // 1 flop
              const FReal wyz = wy * T_of_x[2][j]; // 1 flop
              W4[idxRhs][0][(j-1)*(ORDER-1) + (i-1)] += wxy; // 1 flop
              W4[idxRhs][1][(j-1)*(ORDER-1) + (i-1)] += wxz; // 1 flop
              W4[idxRhs][2][(j-1)*(ORDER-1) + (i-1)] += wyz; // 1 flop
              for (unsigned int k=1; k<ORDER; ++k) {
                const FReal wxyz = wxy * T_of_x[2][k]; // 1 flop
                W8[idxRhs][(k-1)*(ORDER-1)*(ORDER-1) + (j-1)*(ORDER-1) + (i-1)] += wxyz; // 1 flop
              } // flops: (ORDER-1) * 2
            } // flops: (ORDER-1) * (6 + (ORDER-1) * 2)
          } // flops: (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1) * 2))
        } // flops: ... * NRHS 
    } // flops: N * (18 + (ORDER-2) * 6 + (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1) * 2)))

    ////////////////////////////////////////////////////////////////////


    for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){

    // loop over interpolation points
    FReal F2[3][ORDER];
    FReal F4[3][ORDER*ORDER];
    FReal F8[   ORDER*ORDER*ORDER];
    {
        // compute W2: 3 * ORDER*(2*(ORDER-1)-1) flops
        FBlas::gemv(ORDER, ORDER-1, FReal(1.), const_cast<FReal*>(T), W2[idxRhs][0], F2[0]);
        FBlas::gemv(ORDER, ORDER-1, FReal(1.), const_cast<FReal*>(T), W2[idxRhs][1], F2[1]);
        FBlas::gemv(ORDER, ORDER-1, FReal(1.), const_cast<FReal*>(T), W2[idxRhs][2], F2[2]);

        // compute W4: 3 * [ORDER*(ORDER-1)*(2*(ORDER-1)-1) + ORDER*ORDER*(2*(ORDER-1)-1)]
        FReal C[ORDER * (ORDER-1)];
        FBlas::gemmt(ORDER, ORDER-1, ORDER-1, FReal(1.), const_cast<FReal*>(T), ORDER, W4[idxRhs][0], ORDER-1, C,     ORDER);
        FBlas::gemmt(ORDER, ORDER-1, ORDER,   FReal(1.), const_cast<FReal*>(T), ORDER, C,     ORDER,   F4[0], ORDER);
        FBlas::gemmt(ORDER, ORDER-1, ORDER-1, FReal(1.), const_cast<FReal*>(T), ORDER, W4[idxRhs][1], ORDER-1, C,     ORDER);
        FBlas::gemmt(ORDER, ORDER-1, ORDER,   FReal(1.), const_cast<FReal*>(T), ORDER, C,     ORDER,   F4[1], ORDER);
        FBlas::gemmt(ORDER, ORDER-1, ORDER-1, FReal(1.), const_cast<FReal*>(T), ORDER, W4[idxRhs][2], ORDER-1, C,     ORDER);
        FBlas::gemmt(ORDER, ORDER-1, ORDER,   FReal(1.), const_cast<FReal*>(T), ORDER, C,     ORDER,   F4[2], ORDER);

        // compute W8: 3 * (2*(ORDER-1)-1) * [ORDER*(ORDER-1)*(ORDER-1) + ORDER*ORDER*(ORDER-1) + ORDER*ORDER*ORDER]
        FReal D[ORDER * (ORDER-1) * (ORDER-1)];
        FBlas::gemm(ORDER, ORDER-1, (ORDER-1)*(ORDER-1), FReal(1.),	const_cast<FReal*>(T), ORDER, W8[idxRhs], ORDER-1, D, ORDER);
        FReal E[(ORDER-1) * (ORDER-1) * ORDER];
        for (unsigned int s=0; s<perm0.size; ++s)	E[perm0.mni[s]] = D[perm0.imn[s]];
        FReal F[ORDER * (ORDER-1) * ORDER];
        FBlas::gemm(ORDER, ORDER-1, ORDER*(ORDER-1), FReal(1.), const_cast<FReal*>(T), ORDER, E, ORDER-1, F, ORDER);
        FReal G[(ORDER-1) * ORDER * ORDER];
        for (unsigned int s=0; s<perm1.size; ++s)	G[perm1.nij[s]] = F[perm1.jni[s]];
        FReal H[ORDER * ORDER * ORDER];
        FBlas::gemm(ORDER, ORDER-1, ORDER*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER, G, ORDER-1, H, ORDER);
        for (unsigned int s=0; s<perm2.size; ++s)	F8[perm2.ijk[s]] = H[perm2.kij[s]];
    }

    // assemble multipole expansions
    for (unsigned int i=0; i<ORDER; ++i) {
        for (unsigned int j=0; j<ORDER; ++j) {
            for (unsigned int k=0; k<ORDER; ++k) {
                const unsigned int idx = k*ORDER*ORDER + j*ORDER + i;
                multipoleExpansion[2*idxRhs*nnodes + idx] += (W1[idxRhs] +
                                                                     FReal(2.) * (F2[0][i] + F2[1][j] + F2[2][k]) +
                                                                     FReal(4.) * (F4[0][j*ORDER+i] + F4[1][k*ORDER+i] + F4[2][k*ORDER+j]) +
                                                                     FReal(8.) *  F8[idx]) / nnodes; // 11 * ORDER*ORDER*ORDER flops
            }
        }
    }

    } // NRHS

}


///**
// * Particle to moment: application of \f$S_\ell(y,\bar y_n)\f$
// * (anterpolation, it is the transposed interpolation)
// */
//template <int ORDER>
//template <class ContainerClass>
//inline void FChebInterpolator<ORDER>::applyP2M(const FPoint& center,
//																							 const FReal width,
//																							 FReal *const multipoleExpansion,
//																							 const ContainerClass *const sourceParticles) const
//{
//	// set all multipole expansions to zero
//	FBlas::setzero(nnodes, multipoleExpansion);
//
//	// allocate stuff
//	const map_glob_loc map(center, width);
//	FPoint localPosition;
//	FReal T_of_x[ORDER][3];
//	FReal S[3], c1;
//	//
//	FReal xpx,ypy,zpz ;
//	c1 = FReal(8.) / nnodes ; // 1 flop
//	// loop over source particles
//	typename ContainerClass::ConstBasicIterator iter(*sourceParticles);
//	while(iter.hasNotFinished()){
//
//		// map global position to [-1,1]
//		map(iter.data().getPosition(), localPosition); // 15 flops
//
//		// evaluate chebyshev polynomials of source particle: T_o(x_i)
//		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
//		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
//		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
//		xpx = FReal(2.) * localPosition.getX() ; // 1 flop
//		ypy = FReal(2.) * localPosition.getY() ; // 1 flop
//		zpz = FReal(2.) * localPosition.getZ() ; // 1 flop
//
//		for (unsigned int o=2; o<ORDER; ++o) {
//			T_of_x[o][0] = xpx * T_of_x[o-1][0] - T_of_x[o-2][0]; // 2 flops
//			T_of_x[o][1] = ypy * T_of_x[o-1][1] - T_of_x[o-2][1];	// 2 flops
//			T_of_x[o][2] = zpz * T_of_x[o-1][2] - T_of_x[o-2][2]; // 2 flops
//		} // flops: (ORDER-1) * 6
//		
//		// anterpolate
//		const FReal sourceValue = iter.data().getPhysicalValue();
//		for (unsigned int n=0; n<nnodes; ++n) {
//			const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};
//			S[0] = FReal(0.5) + T_of_x[1][0] * T_of_roots[1][j[0]]; // 2 flops 
//			S[1] = FReal(0.5) + T_of_x[1][1] * T_of_roots[1][j[1]]; // 2 flops
//			S[2] = FReal(0.5) + T_of_x[1][2] * T_of_roots[1][j[2]]; // 2 flops
//			for (unsigned int o=2; o<ORDER; ++o) {
//				S[0] += T_of_x[o][0] * T_of_roots[o][j[0]]; // 2 flops
//				S[1] += T_of_x[o][1] * T_of_roots[o][j[1]]; // 2 flops
//				S[2] += T_of_x[o][2] * T_of_roots[o][j[2]]; // 2 flops
//			} // flops: (ORDER-2) * 6
//
//			// gather contributions
//			multipoleExpansion[n]	+= c1 *	S[0] * S[1] * S[2] *	sourceValue; // 4 flops
//		} // flops: ORDER*ORDER*ORDER * (10 + (ORDER-2) * 6)
//
//		// increment source iterator
//		iter.gotoNext();
//	} // flops: M * (18 + (ORDER-1) * 6 + ORDER*ORDER*ORDER * (10 + (ORDER-2) * 6))
//}



/**
 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ (interpolation)
 */
template <int ORDER, class MatrixKernelClass>
template <class ContainerClass>
inline void FChebInterpolator<ORDER,MatrixKernelClass>::applyL2P(const FPoint& center,
                                             const FReal width,
                                             const FReal *const localExpansion,
                                             ContainerClass *const inParticles) const
{
    FReal f1[nLhs];
    FReal W2[nLhs][3][ ORDER-1];
    FReal W4[nLhs][3][(ORDER-1)*(ORDER-1)];
    FReal W8[nLhs][   (ORDER-1)*(ORDER-1)*(ORDER-1)];
    { 
      for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){

      // sum over interpolation points
        f1[idxLhs] = FReal(0.);
        for(unsigned int i=0; i<ORDER-1; ++i)	                   W2[idxLhs][0][i] = W2[idxLhs][1][i] = W2[idxLhs][2][i] = FReal(0.);
        for(unsigned int i=0; i<(ORDER-1)*(ORDER-1); ++i)        W4[idxLhs][0][i] = W4[idxLhs][1][i] = W4[idxLhs][2][i] = FReal(0.);
        for(unsigned int i=0; i<(ORDER-1)*(ORDER-1)*(ORDER-1); ++i)	W8[idxLhs][i] = FReal(0.);

        for (unsigned int idx=0; idx<nnodes; ++idx) {
            const unsigned int i = node_ids[idx][0];
            const unsigned int j = node_ids[idx][1];
            const unsigned int k = node_ids[idx][2];

            f1[idxLhs] += localExpansion[2*idxLhs*nnodes + idx]; // 1 flop

            for (unsigned int l=0; l<ORDER-1; ++l) {
                const FReal wx = T[l*ORDER+i] * localExpansion[2*idxLhs*nnodes + idx]; // 1 flops
                const FReal wy = T[l*ORDER+j] * localExpansion[2*idxLhs*nnodes + idx]; // 1 flops
                const FReal wz = T[l*ORDER+k] * localExpansion[2*idxLhs*nnodes + idx]; // 1 flops
                W2[idxLhs][0][l] += wx; // 1 flops
                W2[idxLhs][1][l] += wy; // 1 flops
                W2[idxLhs][2][l] += wz; // 1 flops
                for (unsigned int m=0; m<ORDER-1; ++m) {
                    const FReal wxy = wx * T[m*ORDER + j]; // 1 flops
                    const FReal wxz = wx * T[m*ORDER + k]; // 1 flops
                    const FReal wyz = wy * T[m*ORDER + k]; // 1 flops
                    W4[idxLhs][0][m*(ORDER-1)+l] += wxy; // 1 flops
                    W4[idxLhs][1][m*(ORDER-1)+l] += wxz; // 1 flops
                    W4[idxLhs][2][m*(ORDER-1)+l] += wyz; // 1 flops
                    for (unsigned int n=0; n<ORDER-1; ++n) {
                        const FReal wxyz = wxy * T[n*ORDER + k]; // 1 flops
                        W8[idxLhs][n*(ORDER-1)*(ORDER-1) + m*(ORDER-1) + l]	+= wxyz; // 1 flops
                    } // (ORDER-1) * 2 flops
                } // (ORDER-1) * (6 + (ORDER-1)*2) flops
            } // (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1)*2)) flops
        } // ORDER*ORDER*ORDER * (1 + (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1)*2))) flops
      } // NLHS
    }


    // loop over particles
    const map_glob_loc map(center, width);
    FPoint localPosition;

    //const FReal*const physicalValues = inParticles->getPhysicalValues();
    const FReal*const positionsX = inParticles->getPositions()[0];
    const FReal*const positionsY = inParticles->getPositions()[1];
    const FReal*const positionsZ = inParticles->getPositions()[2];
    //FReal*const forcesX = inParticles->getForcesX();
    //FReal*const forcesY = inParticles->getForcesY();
    //FReal*const forcesZ = inParticles->getForcesZ();
//    FReal*const potentials = inParticles->getPotentials();

    for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++ idxPart){

      // map global position to [-1,1]
      map(FPoint(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]), localPosition); // 15 flops

      FReal T_of_x[3][ORDER];
      {
        T_of_x[0][0] = FReal(1.); T_of_x[0][1] = localPosition.getX();
        T_of_x[1][0] = FReal(1.); T_of_x[1][1] = localPosition.getY();
        T_of_x[2][0] = FReal(1.); T_of_x[2][1] = localPosition.getZ();
        const FReal x2 = FReal(2.) * T_of_x[0][1]; // 1 flop
        const FReal y2 = FReal(2.) * T_of_x[1][1]; // 1 flop
        const FReal z2 = FReal(2.) * T_of_x[2][1]; // 1 flop
        for (unsigned int j=2; j<ORDER; ++j) {
          T_of_x[0][j] = x2 * T_of_x[0][j-1] - T_of_x[0][j-2]; // 2 flops
          T_of_x[1][j] = y2 * T_of_x[1][j-1] - T_of_x[1][j-2]; // 2 flops
          T_of_x[2][j] = z2 * T_of_x[2][j-1] - T_of_x[2][j-2]; // 2 flops
        }
      }

      for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){      
        // distribution over potential components:
        // We sum the multidim contribution of PhysValue
        // This was originally done at M2L step but moved here 
        // because their storage is required by the force computation.
        // In fact : f_{ik}(x)=w_j(x) \nabla_{x_i} K_{ij}(x,y)w_j(y))
        const unsigned int idxPot = idxLhs / nPV; 

        FReal*const potentials = inParticles->getPotentials(idxPot);

        // interpolate and increment target value
        FReal targetValue = potentials[idxPart];
        {
          FReal f2, f4, f8;
          {
            f2 = f4 = f8 = FReal(0.);
            for (unsigned int l=1; l<ORDER; ++l) {
              f2 +=
                T_of_x[0][l] * W2[idxLhs][0][l-1] +
                T_of_x[1][l] * W2[idxLhs][1][l-1] +
                T_of_x[2][l] * W2[idxLhs][2][l-1]; // 6 flops
              for (unsigned int m=1; m<ORDER; ++m) {
                f4 +=
                  T_of_x[0][l] * T_of_x[1][m] * W4[idxLhs][0][(m-1)*(ORDER-1)+(l-1)] +
                  T_of_x[0][l] * T_of_x[2][m] * W4[idxLhs][1][(m-1)*(ORDER-1)+(l-1)] +
                  T_of_x[1][l] * T_of_x[2][m] * W4[idxLhs][2][(m-1)*(ORDER-1)+(l-1)]; // 9 flops
                for (unsigned int n=1; n<ORDER; ++n) {
                  f8 +=
                    T_of_x[0][l] * T_of_x[1][m] * T_of_x[2][n] *
                    W8[idxLhs][(n-1)*(ORDER-1)*(ORDER-1) + (m-1)*(ORDER-1) + (l-1)];
                } // ORDER * 4 flops
              } // ORDER * (9 + ORDER * 4) flops
            } // ORDER * (ORDER * (9 + ORDER * 4)) flops
          }
          targetValue = (f1[idxLhs] + FReal(2.)*f2 + FReal(4.)*f4 + FReal(8.)*f8) / nnodes; // 7 flops
        } // 7 + ORDER * (ORDER * (9 + ORDER * 4)) flops

          // set potential
        potentials[idxPart] += (targetValue);
      } // NLHS
    } // N * (7 + ORDER * (ORDER * (9 + ORDER * 4))) flops
}


//	FReal F2[3][ORDER-1];
//	FBlas::gemtv(ORDER, ORDER-1, FReal(1.), const_cast<FReal*>(T), const_cast<FReal*>(localExpansion), F2[0]);
//	for (unsigned int i=1; i<ORDER*ORDER; ++i)
//		FBlas::gemtva(ORDER, ORDER-1, FReal(1.), const_cast<FReal*>(T),
//									const_cast<FReal*>(localExpansion) + ORDER*i, F2[0]);
//	for (unsigned int i=0; i<ORDER-1; ++i)
//		std::cout << W2[0][i] << "\t" << F2[0][i] << std::endl;

//	FReal F2[(ORDER-1) * ORDER*ORDER];
//	FBlas::gemtm(ORDER, ORDER-1, ORDER*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER,
//							 const_cast<FReal*>(localExpansion), ORDER, F2, ORDER-1);
//	FReal F[ORDER-1]; FBlas::setzero(ORDER-1, F);
//	for (unsigned int i=0; i<ORDER-1; ++i)
//		for (unsigned int j=0; j<ORDER*ORDER; ++j) F[i] += F2[j*(ORDER-1) + i];
//	for (unsigned int i=0; i<ORDER-1; ++i)
//		std::cout << W2[0][i] << "\t" << F[i] << std::endl;


///**
// * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ (interpolation)
// */
//template <int ORDER>
//template <class ContainerClass>
//inline void FChebInterpolator<ORDER>::applyL2P(const FPoint& center,
//																							 const FReal width,
//																							 const FReal *const localExpansion,
//																							 ContainerClass *const localParticles) const
//{
//	// allocate stuff
//	const map_glob_loc map(center, width);
//	FPoint localPosition;
//	FReal T_of_x[ORDER][3];
//	FReal xpx,ypy,zpz ;
//	FReal S[3],c1;
//	//
//	c1 = FReal(8.) / nnodes ;
//	typename ContainerClass::BasicIterator iter(*localParticles);
//	while(iter.hasNotFinished()){
//			
//		// map global position to [-1,1]
//		map(iter.data().getPosition(), localPosition); // 15 flops
//
//		// evaluate chebyshev polynomials of source particle: T_o(x_i)
//		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
//		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
//		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
//		xpx = FReal(2.) * localPosition.getX() ; // 1 flop
//		ypy = FReal(2.) * localPosition.getY() ; // 1 flop
//		zpz = FReal(2.) * localPosition.getZ() ; // 1 flop
//		for (unsigned int o=2; o<ORDER; ++o) {
//			T_of_x[o][0] = xpx * T_of_x[o-1][0] - T_of_x[o-2][0]; // 2 flop
//			T_of_x[o][1] = ypy * T_of_x[o-1][1] - T_of_x[o-2][1]; // 2 flop
//			T_of_x[o][2] = zpz * T_of_x[o-1][2] - T_of_x[o-2][2]; // 2 flop
//		} // (ORDER-2) * 6 flops
//
//		// interpolate and increment target value
//		FReal targetValue = iter.data().getPotential();
//		for (unsigned int n=0; n<nnodes; ++n) {
//			const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};
//			S[0] = T_of_x[1][0] * T_of_roots[1][j[0]]; // 1 flops
//			S[1] = T_of_x[1][1] * T_of_roots[1][j[1]]; // 1 flops
//			S[2] = T_of_x[1][2] * T_of_roots[1][j[2]]; // 1 flops
//			for (unsigned int o=2; o<ORDER; ++o) {
//				S[0] += T_of_x[o][0] * T_of_roots[o][j[0]]; // 2 flops
//				S[1] += T_of_x[o][1] * T_of_roots[o][j[1]]; // 2 flops
//				S[2] += T_of_x[o][2] * T_of_roots[o][j[2]]; // 2 flops
//			} // (ORDER-2) * 6 flops 
//			// gather contributions
//			S[0] += FReal(0.5); // 1 flops
//			S[1] += FReal(0.5); // 1 flops
//			S[2] += FReal(0.5); // 1 flops
//			targetValue	+= S[0] * S[1] * S[2] * localExpansion[n]; // 4 flops
//		} // ORDER*ORDER*ORDER * (10 + (ORDER-2) * 6) flops
//		// scale
//		targetValue *= c1; // 1 flops
//
//		// set potential
//		iter.data().setPotential(targetValue);
//
//		// increment target iterator
//		iter.gotoNext();
//	} // N * ORDER*ORDER*ORDER * (10 + (ORDER-2) * 6) flops
//}






/**
 * Local to particle operation: application of \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
 */
template <int ORDER, class MatrixKernelClass>
template <class ContainerClass>
inline void FChebInterpolator<ORDER,MatrixKernelClass>::applyL2PGradient(const FPoint& center,
                                                                         const FReal width,
                                                                         const FReal *const localExpansion,
                                                                         ContainerClass *const inParticles) const
{
    ////////////////////////////////////////////////////////////////////
    // TENSOR-PRODUCT INTERPOLUTION NOT IMPLEMENTED YET HERE!!! ////////
    ////////////////////////////////////////////////////////////////////

    // setup local to global mapping
    const map_glob_loc map(center, width);
    FPoint Jacobian;
    map.computeJacobian(Jacobian);
    const FReal jacobian[3] = {Jacobian.getX(), Jacobian.getY(), Jacobian.getZ()};
    FPoint localPosition;
    FReal T_of_x[ORDER][3];
    FReal U_of_x[ORDER][3];
    FReal P[3];

//    const FReal*const physicalValues = inParticles->getPhysicalValues();
    const FReal*const positionsX = inParticles->getPositions()[0];
    const FReal*const positionsY = inParticles->getPositions()[1];
    const FReal*const positionsZ = inParticles->getPositions()[2];
//    FReal*const forcesX = inParticles->getForcesX();
//    FReal*const forcesY = inParticles->getForcesY();
//    FReal*const forcesZ = inParticles->getForcesZ();
    //FReal*const potentials = inParticles->getPotentials();

    for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++ idxPart){

        // map global position to [-1,1]
        map(FPoint(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]), localPosition);

        // evaluate chebyshev polynomials of source particle
        // T_0(x_i) and T_1(x_i)
        T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
        T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
        T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
        // U_0(x_i) and U_1(x_i)
        U_of_x[0][0] = FReal(1.);	U_of_x[1][0] = localPosition.getX() * FReal(2.);
        U_of_x[0][1] = FReal(1.);	U_of_x[1][1] = localPosition.getY() * FReal(2.);
        U_of_x[0][2] = FReal(1.);	U_of_x[1][2] = localPosition.getZ() * FReal(2.);
        for (unsigned int o=2; o<ORDER; ++o) {
            // T_o(x_i)
            T_of_x[o][0] = FReal(2.)*localPosition.getX()*T_of_x[o-1][0] - T_of_x[o-2][0];
            T_of_x[o][1] = FReal(2.)*localPosition.getY()*T_of_x[o-1][1] - T_of_x[o-2][1];
            T_of_x[o][2] = FReal(2.)*localPosition.getZ()*T_of_x[o-1][2] - T_of_x[o-2][2];
            // U_o(x_i)
            U_of_x[o][0] = FReal(2.)*localPosition.getX()*U_of_x[o-1][0] - U_of_x[o-2][0];
            U_of_x[o][1] = FReal(2.)*localPosition.getY()*U_of_x[o-1][1] - U_of_x[o-2][1];
            U_of_x[o][2] = FReal(2.)*localPosition.getZ()*U_of_x[o-1][2] - U_of_x[o-2][2];
        }

        // scale, because dT_o/dx = oU_{o-1}
        for (unsigned int o=2; o<ORDER; ++o) {
            U_of_x[o-1][0] *= FReal(o);
            U_of_x[o-1][1] *= FReal(o);
            U_of_x[o-1][2] *= FReal(o);
        }

        // apply P and increment forces
        FReal forces[nLhs][3]; 
        for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs)
          for (unsigned int i=0; i<3; ++i)
            forces[idxLhs][i] = FReal(0.);

        for (unsigned int n=0; n<nnodes; ++n) {

            // tensor indices of chebyshev nodes
            const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};

            // f0 component //////////////////////////////////////
            P[0] = U_of_x[0][0] * T_of_roots[1][j[0]];
            P[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
            P[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
            for (unsigned int o=2; o<ORDER; ++o) {
                P[0] += U_of_x[o-1][0] * T_of_roots[o][j[0]];
                P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]];
                P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]];
            }
            P[0] *= FReal(2.);
            P[1] *= FReal(2.); P[1] += FReal(1.);
            P[2] *= FReal(2.); P[2] += FReal(1.);
            for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs)
              forces[idxLhs][0]	+= P[0] * P[1] * P[2] * localExpansion[2*idxLhs*nnodes + n];

            // f1 component //////////////////////////////////////
            P[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
            P[1] = U_of_x[0][1] * T_of_roots[1][j[1]];
            P[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
            for (unsigned int o=2; o<ORDER; ++o) {
                P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]];
                P[1] += U_of_x[o-1][1] * T_of_roots[o][j[1]];
                P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]];
            }
            P[0] *= FReal(2.); P[0] += FReal(1.);
            P[1] *= FReal(2.);
            P[2] *= FReal(2.); P[2] += FReal(1.);
            for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs)
            forces[idxLhs][1]	+= P[0] * P[1] * P[2] * localExpansion[2*idxLhs*nnodes + n];

            // f2 component //////////////////////////////////////
            P[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
            P[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
            P[2] = U_of_x[0][2] * T_of_roots[1][j[2]];
            for (unsigned int o=2; o<ORDER; ++o) {
                P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]];
                P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]];
                P[2] += U_of_x[o-1][2] * T_of_roots[o][j[2]];
            }
            P[0] *= FReal(2.); P[0] += FReal(1.);
            P[1] *= FReal(2.); P[1] += FReal(1.);
            P[2] *= FReal(2.);
            for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs)
              forces[idxLhs][2]	+= P[0] * P[1] * P[2] * localExpansion[2*idxLhs*nnodes + n];
        }

        for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
          const unsigned int idxPot = idxLhs / nPV; 
          const unsigned int idxPV  = idxLhs % nPV; 

          // scale forces
          forces[idxLhs][0] *= jacobian[0] / nnodes;
          forces[idxLhs][1] *= jacobian[1] / nnodes;
          forces[idxLhs][2] *= jacobian[2] / nnodes;

          // get pointers to PhysValues and force components
          const FReal*const physicalValues = inParticles->getPhysicalValues(idxPV);
          FReal*const forcesX = inParticles->getForcesX(idxPot);
          FReal*const forcesY = inParticles->getForcesY(idxPot);
          FReal*const forcesZ = inParticles->getForcesZ(idxPot);

          // set computed forces
          forcesX[idxPart] += forces[idxLhs][0] * physicalValues[idxPart];
          forcesY[idxPart] += forces[idxLhs][1] * physicalValues[idxPart];
          forcesZ[idxPart] += forces[idxLhs][2] * physicalValues[idxPart];
        }
    }
}


/**
 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ and
 * \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
 */
template <int ORDER, class MatrixKernelClass>
template <class ContainerClass>
inline void FChebInterpolator<ORDER,MatrixKernelClass>::applyL2PTotal(const FPoint& center,
                                                                      const FReal width,
                                                                      const FReal *const localExpansion,
                                                                      ContainerClass *const inParticles) const
{
    FReal f1[nLhs];
    FReal W2[nLhs][3][ ORDER-1];
    FReal W4[nLhs][3][(ORDER-1)*(ORDER-1)];
    FReal W8[nLhs][   (ORDER-1)*(ORDER-1)*(ORDER-1)];

    //{ // sum over interpolation points
    //	f1 = FReal(0.);
    //	for(unsigned int i=0; i<ORDER-1; ++i)	                   W2[0][i] = W2[1][i] = W2[2][i] = FReal(0.);
    //	for(unsigned int i=0; i<(ORDER-1)*(ORDER-1); ++i)        W4[0][i] = W4[1][i] = W4[2][i] = FReal(0.);
    //	for(unsigned int i=0; i<(ORDER-1)*(ORDER-1)*(ORDER-1); ++i)	W8[i] = FReal(0.);
    //
    //	for (unsigned int idx=0; idx<nnodes; ++idx) {
    //		const unsigned int i = node_ids[idx][0];
    //		const unsigned int j = node_ids[idx][1];
    //		const unsigned int k = node_ids[idx][2];
    //
    //		f1 += localExpansion[idx]; // 1 flop
    //
    //		for (unsigned int l=0; l<ORDER-1; ++l) {
    //			const FReal wx = T[l*ORDER+i] * localExpansion[idx]; // 1 flops
    //			const FReal wy = T[l*ORDER+j] * localExpansion[idx]; // 1 flops
    //			const FReal wz = T[l*ORDER+k] * localExpansion[idx]; // 1 flops
    //			W2[0][l] += wx; // 1 flops
    //			W2[1][l] += wy; // 1 flops
    //			W2[2][l] += wz; // 1 flops
    //			for (unsigned int m=0; m<ORDER-1; ++m) {
    //				const FReal wxy = wx * T[m*ORDER + j]; // 1 flops
    //				const FReal wxz = wx * T[m*ORDER + k]; // 1 flops
    //				const FReal wyz = wy * T[m*ORDER + k]; // 1 flops
    //				W4[0][m*(ORDER-1)+l] += wxy; // 1 flops
    //				W4[1][m*(ORDER-1)+l] += wxz; // 1 flops
    //				W4[2][m*(ORDER-1)+l] += wyz; // 1 flops
    //				for (unsigned int n=0; n<ORDER-1; ++n) {
    //					const FReal wxyz = wxy * T[n*ORDER + k]; // 1 flops
    //					W8[n*(ORDER-1)*(ORDER-1) + m*(ORDER-1) + l]	+= wxyz; // 1 flops
    //				} // (ORDER-1) * 2 flops
    //			} // (ORDER-1) * (6 + (ORDER-1)*2) flops
    //		} // (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1)*2)) flops
    //
    //	} // ORDER*ORDER*ORDER * (1 + (ORDER-1) * (6 + (ORDER-1) * (6 + (ORDER-1)*2))) flops
    //
    //}

    {
        // for W2
        FReal lE[nnodes];
        FReal F2[(ORDER-1) * ORDER*ORDER];
        // for W4
        FReal F4[ORDER * ORDER*(ORDER-1)];
        FReal G4[(ORDER-1) * ORDER*(ORDER-1)];
        // for W8
        FReal G8[ORDER * (ORDER-1)*(ORDER-1)];

        // sum local expansions
        for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
          f1[idxLhs] = FReal(0.);
          for (unsigned int idx=0; idx<nnodes; ++idx)	f1[idxLhs] += localExpansion[2*idxLhs*nnodes + idx]; // 1 flop
        
        //////////////////////////////////////////////////////////////////
        // IMPORTANT: NOT CHANGE ORDER OF COMPUTATIONS!!! ////////////////
        //////////////////////////////////////////////////////////////////

        // W2[0] ///////////////// (ORDER-1)*ORDER*ORDER * 2*ORDER
        FBlas::gemtm(ORDER, ORDER-1, ORDER*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER,
                                 const_cast<FReal*>(localExpansion) + 2*idxLhs*nnodes, ORDER, F2, ORDER-1);
        for (unsigned int l=0; l<ORDER-1; ++l) { W2[idxLhs][0][l] = F2[l];
            for (unsigned int j=1; j<ORDER*ORDER; ++j) W2[idxLhs][0][l] += F2[j*(ORDER-1) + l];	}
        // W4[0] ///////////////// ORDER*(ORDER-1)*(ORDER-1) + 2*ORDER
        perm5.permute(F2, F4);
        FBlas::gemtm(ORDER, ORDER-1, ORDER*(ORDER-1), FReal(1.), const_cast<FReal*>(T), ORDER, F4, ORDER, G4, ORDER-1);
        for (unsigned int l=0; l<ORDER-1; ++l)
            for (unsigned int m=0; m<ORDER-1; ++m) { W4[idxLhs][0][m*(ORDER-1)+l] = G4[l*ORDER*(ORDER-1) + m];
                for (unsigned int k=1; k<ORDER; ++k) W4[idxLhs][0][m*(ORDER-1)+l] += G4[l*ORDER*(ORDER-1) + k*(ORDER-1) + m];	}
        // W8 //////////////////// (ORDER-1)*(ORDER-1)*(ORDER-1) * (2*ORDER-1)
        perm8.permute(G4, G8);
        FReal F8[(ORDER-1)*(ORDER-1)*(ORDER-1)];
        FBlas::gemtm(ORDER, ORDER-1, (ORDER-1)*(ORDER-1), FReal(1.), const_cast<FReal*>(T), ORDER, G8, ORDER, F8, ORDER-1);
        perm9.permute(F8, W8[idxLhs]);
        // W4[1] ///////////////// ORDER*(ORDER-1)*(ORDER-1) + 2*ORDER
        perm6.permute(F2, F4);
        FBlas::gemtm(ORDER, ORDER-1, (ORDER-1)*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER, F4, ORDER, G4, ORDER-1);
        for (unsigned int l=0; l<ORDER-1; ++l)
            for (unsigned int n=0; n<ORDER-1; ++n) { W4[idxLhs][1][n*(ORDER-1)+l] = G4[l*(ORDER-1) + n];
                for (unsigned int j=1; j<ORDER; ++j) W4[idxLhs][1][n*(ORDER-1)+l] += G4[j*(ORDER-1)*(ORDER-1) + l*(ORDER-1) + n];	}
        // W2[1] ///////////////// (ORDER-1)*ORDER*ORDER * 2*ORDER
        perm3.permute(localExpansion + 2*idxLhs*nnodes, lE);
        FBlas::gemtm(ORDER, ORDER-1, ORDER*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER, lE, ORDER, F2, ORDER-1);
        for (unsigned int i=0; i<ORDER-1; ++i) { W2[idxLhs][1][i] = F2[i];
            for (unsigned int j=1; j<ORDER*ORDER; ++j) W2[idxLhs][1][i] += F2[j*(ORDER-1) + i]; }
        // W4[2] ///////////////// ORDER*(ORDER-1)*(ORDER-1) + 2*ORDER
        perm7.permute(F2, F4);
        FBlas::gemtm(ORDER, ORDER-1, (ORDER-1)*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER, F4, ORDER, G4, ORDER-1);
        for (unsigned int m=0; m<ORDER-1; ++m)
            for (unsigned int n=0; n<ORDER-1; ++n) { W4[idxLhs][2][n*(ORDER-1)+m] = G4[m*ORDER*(ORDER-1) + n];
                for (unsigned int i=1; i<ORDER; ++i) W4[idxLhs][2][n*(ORDER-1)+m] += G4[m*ORDER*(ORDER-1) + i*(ORDER-1) + n];	}
        // W2[2] ///////////////// (ORDER-1)*ORDER*ORDER * 2*ORDER
        perm4.permute(localExpansion + 2*idxLhs*nnodes, lE);
        FBlas::gemtm(ORDER, ORDER-1, ORDER*ORDER, FReal(1.), const_cast<FReal*>(T), ORDER, lE, ORDER, F2, ORDER-1);
        for (unsigned int i=0; i<ORDER-1; ++i) { W2[idxLhs][2][i] = F2[i];
            for (unsigned int j=1; j<ORDER*ORDER; ++j) W2[idxLhs][2][i] += F2[j*(ORDER-1) + i]; }
    }

    }// NLHS

    // loop over particles
    const map_glob_loc map(center, width);
    FPoint Jacobian;
    map.computeJacobian(Jacobian); // 6 flops
    const FReal jacobian[3] = {Jacobian.getX(), Jacobian.getY(), Jacobian.getZ()};
    FPoint localPosition;

//    const FReal*const physicalValues = inParticles->getPhysicalValues();
    const FReal*const positionsX = inParticles->getPositions()[0];
    const FReal*const positionsY = inParticles->getPositions()[1];
    const FReal*const positionsZ = inParticles->getPositions()[2];
//    FReal*const forcesX = inParticles->getForcesX();
//    FReal*const forcesY = inParticles->getForcesY();
//    FReal*const forcesZ = inParticles->getForcesZ();
//    FReal*const potentials = inParticles->getPotentials();

    for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++ idxPart){

        // map global position to [-1,1]
        map(FPoint(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]), localPosition); // 15 flops

        FReal U_of_x[3][ORDER];
        FReal T_of_x[3][ORDER];
        {
            T_of_x[0][0] = FReal(1.); T_of_x[0][1] = localPosition.getX();
            T_of_x[1][0] = FReal(1.); T_of_x[1][1] = localPosition.getY();
            T_of_x[2][0] = FReal(1.); T_of_x[2][1] = localPosition.getZ();
            const FReal x2 = FReal(2.) * T_of_x[0][1]; // 1 flop
            const FReal y2 = FReal(2.) * T_of_x[1][1]; // 1 flop
            const FReal z2 = FReal(2.) * T_of_x[2][1]; // 1 flop
            U_of_x[0][0] = FReal(1.);	U_of_x[0][1] = x2;
            U_of_x[1][0] = FReal(1.);	U_of_x[1][1] = y2;
            U_of_x[2][0] = FReal(1.);	U_of_x[2][1] = z2;
            for (unsigned int j=2; j<ORDER; ++j) {
                T_of_x[0][j] = x2 * T_of_x[0][j-1] - T_of_x[0][j-2]; // 2 flops
                T_of_x[1][j] = y2 * T_of_x[1][j-1] - T_of_x[1][j-2]; // 2 flops
                T_of_x[2][j] = z2 * T_of_x[2][j-1] - T_of_x[2][j-2]; // 2 flops
                U_of_x[0][j] = x2 * U_of_x[0][j-1] - U_of_x[0][j-2]; // 2 flops
                U_of_x[1][j] = y2 * U_of_x[1][j-1] - U_of_x[1][j-2]; // 2 flops
                U_of_x[2][j] = z2 * U_of_x[2][j-1] - U_of_x[2][j-2]; // 2 flops
            }
            // scale, because dT_j/dx = jU_{j-1}
            for (unsigned int j=2; j<ORDER; ++j) {
                U_of_x[0][j-1] *= FReal(j); // 1 flops
                U_of_x[1][j-1] *= FReal(j); // 1 flops
                U_of_x[2][j-1] *= FReal(j); // 1 flops
            }
        } // 3 + (ORDER-2)*15

        // apply P and increment forces
        FReal potential[nLhs];
        FReal forces[nLhs][3]; 
        for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
          potential[idxLhs]= FReal(0.);
          for (unsigned int i=0; i<3; ++i)
            forces[idxLhs][i] = FReal(0.);
        }

        for( int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){

        {
            FReal f2[4], f4[4], f8[4];
            for (unsigned int i=0; i<4; ++i) f2[i] = f4[i] = f8[i] = FReal(0.);
            {
                for (unsigned int l=1; l<ORDER; ++l) {
                    const FReal w2[3] = {W2[idxLhs][0][l-1], W2[idxLhs][1][l-1], W2[idxLhs][2][l-1]};
                    f2[0] += T_of_x[0][l  ] * w2[0] + T_of_x[1][l] * w2[1] + T_of_x[2][l] * w2[2]; // 6 flops
                    f2[1] += U_of_x[0][l-1] * w2[0]; // 2 flops
                    f2[2] += U_of_x[1][l-1] * w2[1]; // 2 flops
                    f2[3] += U_of_x[2][l-1] * w2[2]; // 2 flops
                    for (unsigned int m=1; m<ORDER; ++m) {
                        const unsigned int w4idx = (m-1)*(ORDER-1)+(l-1);
                        const FReal w4[3] = {W4[idxLhs][0][w4idx], W4[idxLhs][1][w4idx], W4[idxLhs][2][w4idx]};
                        f4[0] +=
                            T_of_x[0][l] * T_of_x[1][m] * w4[0] +
                            T_of_x[0][l] * T_of_x[2][m] * w4[1] +
                            T_of_x[1][l] * T_of_x[2][m] * w4[2]; // 9 flops
                        f4[1] += U_of_x[0][l-1] * T_of_x[1][m]   * w4[0] + U_of_x[0][l-1] * T_of_x[2][m]   * w4[1]; // 6 flops
                        f4[2] += T_of_x[0][l]   * U_of_x[1][m-1] * w4[0] + U_of_x[1][l-1] * T_of_x[2][m]   * w4[2]; // 6 flops
                        f4[3] += T_of_x[0][l]   * U_of_x[2][m-1] * w4[1] + T_of_x[1][l]   * U_of_x[2][m-1] * w4[2]; // 6 flops
                        for (unsigned int n=1; n<ORDER; ++n) {
                            const FReal w8 = W8[idxLhs][(n-1)*(ORDER-1)*(ORDER-1) + (m-1)*(ORDER-1) + (l-1)];
                            f8[0] += T_of_x[0][l]   * T_of_x[1][m]   * T_of_x[2][n]   * w8; // 4 flops
                            f8[1] += U_of_x[0][l-1] * T_of_x[1][m]   * T_of_x[2][n]   * w8; // 4 flops
                            f8[2] += T_of_x[0][l]   * U_of_x[1][m-1] * T_of_x[2][n]   * w8; // 4 flops
                            f8[3] += T_of_x[0][l]   * T_of_x[1][m]   * U_of_x[2][n-1] * w8; // 4 flops
                        } // (ORDER-1) * 16 flops
                    } // (ORDER-1) * (27 + (ORDER-1) * 16) flops
                } // (ORDER-1) * ((ORDER-1) * (27 + (ORDER-1) * 16)) flops
            }
            potential[idxLhs] = (f1[idxLhs] + FReal(2.)*f2[0] + FReal(4.)*f4[0] + FReal(8.)*f8[0]) / nnodes; // 7 flops
            forces[idxLhs][0] = (     FReal(2.)*f2[1] + FReal(4.)*f4[1] + FReal(8.)*f8[1]) * jacobian[0] / nnodes; // 7 flops
            forces[idxLhs][1] = (     FReal(2.)*f2[2] + FReal(4.)*f4[2] + FReal(8.)*f8[2]) * jacobian[1] / nnodes; // 7 flops
            forces[idxLhs][2] = (     FReal(2.)*f2[3] + FReal(4.)*f4[3] + FReal(8.)*f8[3]) * jacobian[2] / nnodes; // 7 flops
        } // 28 + (ORDER-1) * ((ORDER-1) * (27 + (ORDER-1) * 16)) flops

        const  int idxPot = idxLhs / nPV;
        const  int idxPV  = idxLhs % nPV;

        // get potentials, physValues and forces components 
        const FReal*const physicalValues = inParticles->getPhysicalValues(idxPV);
        FReal*const forcesX = inParticles->getForcesX(idxPot);
        FReal*const forcesY = inParticles->getForcesY(idxPot);
        FReal*const forcesZ = inParticles->getForcesZ(idxPot);
        FReal*const potentials = inParticles->getPotentials(idxPot);

        // set computed potential
        potentials[idxPart] += (potential[idxLhs]); // 1 flop

        // set computed forces
        forcesX[idxPart] += forces[idxLhs][0] * physicalValues[idxPart];
        forcesY[idxPart] += forces[idxLhs][1] * physicalValues[idxPart];
        forcesZ[idxPart] += forces[idxLhs][2] * physicalValues[idxPart]; // 6 flops

        }// NLHS

    } // N * (38 + (ORDER-2)*15 + (ORDER-1)*((ORDER-1) * (27 + (ORDER-1) * 16))) + 6 flops
}


///**
// * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ and
// * \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
// */
//template <int ORDER>
//template <class ContainerClass>
//inline void FChebInterpolator<ORDER>::applyL2PTotal(const FPoint& center,
//																										const FReal width,
//																										const FReal *const localExpansion,
//																										ContainerClass *const localParticles) const
//{
//	// setup local to global mapping
//	const map_glob_loc map(center, width);
//	FPoint Jacobian;
//	map.computeJacobian(Jacobian); // 6 flops
//	const FReal jacobian[3] = {Jacobian.getX(), Jacobian.getY(), Jacobian.getZ()}; 
//	FPoint localPosition;
//	FReal T_of_x[ORDER][3];
//	FReal U_of_x[ORDER][3];
//	FReal P[6];
//	//
//	FReal xpx,ypy,zpz ;
//	FReal c1 = FReal(8.0) / nnodes; // 1 flop
//	//
//	typename ContainerClass::BasicIterator iter(*localParticles);
//	while(iter.hasNotFinished()){
//			
//		// map global position to [-1,1]
//		map(iter.data().getPosition(), localPosition); // 15 flops
//			
//		// evaluate chebyshev polynomials of source particle
//		// T_0(x_i) and T_1(x_i)
//		xpx = FReal(2.) * localPosition.getX(); // 1 flop
//		ypy = FReal(2.) * localPosition.getY(); // 1 flop
//		zpz = FReal(2.) * localPosition.getZ(); // 1 flop
//		//
//		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
//		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
//		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
//		U_of_x[0][0] = FReal(1.);	U_of_x[1][0] = xpx;
//		U_of_x[0][1] = FReal(1.);	U_of_x[1][1] = ypy;
//		U_of_x[0][2] = FReal(1.);	U_of_x[1][2] = zpz;
//		for (unsigned int o=2; o<ORDER; ++o) {
//			T_of_x[o][0] = xpx * T_of_x[o-1][0] - T_of_x[o-2][0]; // 2 flops 
//			T_of_x[o][1] = ypy * T_of_x[o-1][1] - T_of_x[o-2][1]; // 2 flops
//			T_of_x[o][2] = zpz * T_of_x[o-1][2] - T_of_x[o-2][2]; // 2 flops
//			U_of_x[o][0] = xpx * U_of_x[o-1][0] - U_of_x[o-2][0]; // 2 flops
//			U_of_x[o][1] = ypy * U_of_x[o-1][1] - U_of_x[o-2][1]; // 2 flops
//			U_of_x[o][2] = zpz * U_of_x[o-1][2] - U_of_x[o-2][2]; // 2 flops
//		}
//
//		// scale, because dT_o/dx = oU_{o-1}
//		for (unsigned int o=2; o<ORDER; ++o) {
//			U_of_x[o-1][0] *= FReal(o); // 1 flops
//			U_of_x[o-1][1] *= FReal(o); // 1 flops
//			U_of_x[o-1][2] *= FReal(o); // 1 flops
//		}
//
//		// apply P and increment forces
//		FReal potential = FReal(0.);
//		FReal forces[3] = {FReal(0.), FReal(0.), FReal(0.)};
//		//
//		// Optimization:
//		//   Here we compute 1/2 S and 1/2 P  rather S and F like in the paper
//		for (unsigned int n=0; n<nnodes; ++n) {
//		  
//		  // tensor indices of chebyshev nodes
//		  const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};
//		  //
//		  P[0] = FReal(0.5) + T_of_x[1][0] * T_of_roots[1][j[0]]; // 2 flops 
//		  P[1] = FReal(0.5) + T_of_x[1][1] * T_of_roots[1][j[1]]; // 2 flops
//		  P[2] = FReal(0.5) + T_of_x[1][2] * T_of_roots[1][j[2]]; // 2 flops
//		  P[3] = U_of_x[0][0] * T_of_roots[1][j[0]]; // 1 flop
//		  P[4] = U_of_x[0][1] * T_of_roots[1][j[1]]; // 1 flop
//		  P[5] = U_of_x[0][2] * T_of_roots[1][j[2]]; // 1 flop
//		  for (unsigned int o=2; o<ORDER; ++o) {
//		    P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]]; // 2 flop
//		    P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]]; // 2 flop
//		    P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]]; // 2 flop
//		    P[3] += U_of_x[o-1][0] * T_of_roots[o][j[0]]; // 2 flop
//		    P[4] += U_of_x[o-1][1] * T_of_roots[o][j[1]]; // 2 flop
//		    P[5] += U_of_x[o-1][2] * T_of_roots[o][j[2]]; // 2 flop
//		  }
//		  //
//		  potential	+= P[0] * P[1] * P[2] * localExpansion[n]; // 4 flops
//		  forces[0]	+= P[3] * P[1] * P[2] * localExpansion[n]; // 4 flops
//		  forces[1]	+= P[0] * P[4] * P[2] * localExpansion[n]; // 4 flops
//		  forces[2]	+= P[0] * P[1] * P[5] * localExpansion[n]; // 4 flops
//		}
//		//
//		potential *= c1 ; // 1 flop
//		forces[0] *= jacobian[0] *c1; // 2 flops 
//		forces[1] *= jacobian[1] *c1; // 2 flops
//		forces[2] *= jacobian[2] *c1; // 2 flops
//		// set computed potential
//		iter.data().incPotential(potential); // 1 flop
//		
//		// set computed forces
//		iter.data().incForces(forces[0] * iter.data().getPhysicalValue(),
//													forces[1] * iter.data().getPhysicalValue(),
//													forces[2] * iter.data().getPhysicalValue()); // 6 flops
//
//		// increment iterator
//		iter.gotoNext();
//	}
//}


#endif








        ////struct IMN2MNI {
        ////	enum {size = ORDER * (ORDER-1) * (ORDER-1)};
        ////	unsigned int imn[size], mni[size];
        ////	IMN2MNI() {
        ////		unsigned int counter = 0;
        ////		for (unsigned int i=0; i<ORDER; ++i) {
        ////			for (unsigned int m=0; m<ORDER-1; ++m) {
        ////				for (unsigned int n=0; n<ORDER-1; ++n) {
        ////					imn[counter] = n*(ORDER-1)*ORDER + m*ORDER + i;
        ////					mni[counter] = i*(ORDER-1)*(ORDER-1) + n*(ORDER-1) + m;
        ////					counter++;
        ////				}
        ////			}
        ////		}
        ////	}
        ////} perm0;
        //
        ////for (unsigned int i=0; i<ORDER; ++i) {
        ////	for (unsigned int m=0; m<ORDER-1; ++m) {
        ////		for (unsigned int n=0; n<ORDER-1; ++n) {
        ////			const unsigned int a = n*(ORDER-1)*ORDER + m*ORDER + i;
        ////			const unsigned int b = i*(ORDER-1)*(ORDER-1) + n*(ORDER-1) + m;
        ////			E[b] = D[a];
        ////		}
        ////	}
        ////}

        ////struct JNI2NIJ {
        ////	enum {size = ORDER * ORDER * (ORDER-1)};
        ////	unsigned int jni[size], nij[size];
        ////	JNI2NIJ() {
        ////		unsigned int counter = 0;
        ////		for (unsigned int i=0; i<ORDER; ++i) {
        ////			for (unsigned int j=0; j<ORDER; ++j) {
        ////				for (unsigned int n=0; n<ORDER-1; ++n) {
        ////					jni[counter] = i*(ORDER-1)*ORDER + n*ORDER + j;
        ////					nij[counter] = j*ORDER*(ORDER-1) + i*(ORDER-1) + n;
        ////					counter++;
        ////				}
        ////			}
        ////		}
        ////	}
        ////} perm1;
        //
        ////for (unsigned int i=0; i<ORDER; ++i) {
        ////	for (unsigned int j=0; j<ORDER; ++j) {
        ////		for (unsigned int n=0; n<ORDER-1; ++n) {
        ////			const unsigned int a = i*(ORDER-1)*ORDER + n*ORDER + j;
        ////			const unsigned int b = j*ORDER*(ORDER-1) + i*(ORDER-1) + n;
        ////			G[b] = F[a];
        ////		}
        ////	}
        ////}

        ////struct KIJ2IJK {
        ////	enum {size = ORDER * ORDER * ORDER};
        ////	unsigned int kij[size], ijk[size];
        ////	KIJ2IJK() {
        ////		unsigned int counter = 0;
        ////		for (unsigned int i=0; i<ORDER; ++i) {
        ////			for (unsigned int j=0; j<ORDER; ++j) {
        ////				for (unsigned int k=0; k<ORDER; ++k) {
        ////					kij[counter] = j*ORDER*ORDER + i*ORDER + k;
        ////					ijk[counter] = k*ORDER*ORDER + j*ORDER + i;
        ////					counter++;
        ////				}
        ////			}
        ////		}
        ////	}
        ////} perm2;
        //
        ////for (unsigned int i=0; i<ORDER; ++i) {
        ////	for (unsigned int j=0; j<ORDER; ++j) {
        ////		for (unsigned int k=0; k<ORDER; ++k) {
        ////			const unsigned int a = j*ORDER*ORDER + i*ORDER + k;
        ////			const unsigned int b = k*ORDER*ORDER + j*ORDER + i;
        ////			F8[b] = H[a];
        ////		}
        ////	}
        ////}

        //FReal T_of_y[ORDER * (ORDER-1)];
        //for (unsigned int o=1; o<ORDER; ++o)
        //	for (unsigned int j=0; j<ORDER; ++j)
        //		T_of_y[(o-1)*ORDER + j] = FReal(FChebRoots<ORDER>::T(o, FReal(FChebRoots<ORDER>::roots[j])));

    //struct SumP2M {
    //	unsigned int f2[3][nnodes], f4[3][nnodes];
    //	SumP2M() {
    //		for (unsigned int i=0; i<ORDER; ++i) {
    //			for (unsigned int j=0; j<ORDER; ++j) {
    //				for (unsigned int k=0; k<ORDER; ++k) {
    //					const unsigned int idx = k*ORDER*ORDER + j*ORDER + i;
    //					f2[0][idx] = i;
    //					f2[1][idx] = j;
    //					f2[2][idx] = k;
    //					f4[0][idx] = j*ORDER+i;
    //					f4[1][idx] = k*ORDER+i;
    //					f4[2][idx] = k*ORDER+j;
    //				}
    //			}
    //		}
    //	}
    //} idx0;
    //
    //for (unsigned int i=0; i<nnodes; ++i)
    //	multipoleExpansion[i] = (W1 +
    //													 FReal(2.) * (F2[0][idx0.f2[0][i]] + F2[1][idx0.f2[1][i]] + F2[2][idx0.f2[2][i]]) +
    //													 FReal(4.) * (F4[0][idx0.f4[0][i]] + F4[1][idx0.f4[1][i]] + F4[2][idx0.f4[2][i]]) +
    //													 FReal(8.) *  F8[i]) / nnodes;
