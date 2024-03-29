#ifndef FADAPTUNIFKERNEL_HPP
#define FADAPTUNIFKERNEL_HPP
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
// Keep in private GIT
// @SCALFMM_PRIVATE
#include "Utils/FGlobal.hpp"

#include "Utils/FPoint.hpp"

#include "Adaptive/FAdaptiveCell.hpp"
#include "Adaptive/FAdaptiveKernelWrapper.hpp"
#include "Adaptive/FAbstractAdaptiveKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"
#include "Kernels/Uniform/FUnifM2LHandler.hpp"

class FTreeCoordinate;

// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// ================

// for verbosity only!!!
//#define COUNT_BLOCKED_INTERACTIONS

// if timings should be logged
//#define LOG_TIMINGS

/**
 * @author O. Coulaud
 * @class FAdaptUnifKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Lagrange interpolation based FMM operators.
 * It implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are
 * required by the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */

template<class FReal, class CellClass, class ContainerClass, class MatrixKernelClass, int ORDER, int NVALS = 1>
class FAdaptiveUnifKernel : public FUnifKernel<FReal,CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
        , public FAbstractAdaptiveKernel<CellClass, ContainerClass> {
    //
    typedef FUnifKernel<FReal,CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	KernelBaseClass;

    enum {order = ORDER,
          nnodes = TensorTraits<ORDER>::nnodes};

    /// Needed for M2L operator
    //  // If we choose to  pre-assemble adaptive M2L operators
    //  // then we need to provide an adaptive M2L handler
    //  // and the transfer in Fourier space is straightforward.
    //  typedef FUnifM2LHandler<ORDER,MatrixKernelClass::Type> M2LHandlerClass;
    //  const M2LHandlerClass M2LHandler;

    const MatrixKernelClass *const MatrixKernel;
    int sminM, sminL;
public:

    using KernelBaseClass::P2M;
    using KernelBaseClass::M2M;
    using KernelBaseClass::M2L;
    using KernelBaseClass::finishedLevelM2L;
    using KernelBaseClass::L2L;
    using KernelBaseClass::L2P;
    using KernelBaseClass::P2P;
    using KernelBaseClass::P2PRemote;
    //	/**
    //	 * The constructor initializes all constant attributes and it reads the
    //	 * precomputed and compressed M2L operators from a binary file (an
    //	 * runtime_error is thrown if the required file is not valid).
    //	 */
    FAdaptiveUnifKernel(const int inTreeHeight, const FReal inBoxWidth,
                        const FPoint<FReal>& inBoxCenter, const MatrixKernelClass *const inMatrixKernel, const int &minM, const int &minL) : KernelBaseClass(inTreeHeight, inBoxWidth, inBoxCenter, inMatrixKernel)
      /*, M2LHandler(inMatrixKernel, inTreeHeight, inBoxWidth)*/, MatrixKernel(inMatrixKernel),sminM(minM),sminL(minM)
    {}
    //	/** Copy constructor */
    FAdaptiveUnifKernel(const FAdaptiveUnifKernel& other)
        : KernelBaseClass(other)/*, M2LHandler(other.M2LHandler)*/, MatrixKernel(other.MatrixKernel),sminM(other.sminM),sminL(other.sminL)
    {	}

    //
    //	/** Destructor */
    ~FAdaptiveUnifKernel()
    {
        //this->~KernelBaseClass() ;
    }
    void P2M(CellClass* const pole, const int cellLevel, const ContainerClass* const particles) override {

        const FPoint<FReal> CellCenter(KernelBaseClass::getCellCenter(pole->getCoordinate(),cellLevel));
        const FReal BoxWidthAtLevel = KernelBaseClass::BoxWidth / FMath::pow(2.0,cellLevel);

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1) apply Sy
            KernelBaseClass::Interpolator->applyP2M(CellCenter, BoxWidthAtLevel,
                                                    pole->getMultipole(idxRhs), particles);
            //      // 2) apply Discrete Fourier Transform
            //      M2LHandler.applyZeroPaddingAndDFT(pole->getMultipole(idxRhs),
            //                                        pole->getTransformedMultipole(idxRhs));

        }

    }

    void M2M(CellClass* const pole, const int poleLevel, const CellClass* const subCell, const int subCellLevel) override {

        const FPoint<FReal> subCellCenter(KernelBaseClass::getCellCenter(subCell->getCoordinate(),subCellLevel));
        const FReal subCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0,subCellLevel)));

        const FPoint<FReal> poleCellCenter(KernelBaseClass::getCellCenter(pole->getCoordinate(),poleLevel));
        const FReal poleCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0, poleLevel)));

        //    ////////////////////////////////////////////////////////////////////////////
        //    /// p^6 version
        //    // allocate memory
        //    FReal* subChildParentInterpolator = new FReal [nnodes * nnodes];
        //
        //    // set child info
        //    FPoint<FReal> ChildRoots[nnodes], localChildRoots[nnodes];
        //    FUnifTensor<FReal,ORDER>::setRoots(subCellCenter, subCellWidth, ChildRoots);
        //
        //    // map global position of roots to local position in parent cell
        //    const map_glob_loc<FReal> map(poleCellCenter, poleCellWidth);
        //    for (unsigned int n=0; n<nnodes; ++n)
        //      map(ChildRoots[n], localChildRoots[n]);
        //
        //    // assemble child - parent - interpolator
        //    KernelBaseClass::Interpolator->assembleInterpolator(nnodes, localChildRoots, subChildParentInterpolator);


        ////////////////////////////////////////////////////////////////////////////
        /// p^4 version

        // Set sub-child coords
        FReal globalChildCoords[3][ORDER];
        FUnifTensor<FReal,order>::setPolynomialsRoots(subCellCenter, subCellWidth, globalChildCoords);

        // Map global position of sub-child nodes to [-1,1]
        FReal localChildCoords[3][ORDER];
        const map_glob_loc<FReal> map(poleCellCenter, poleCellWidth);
        FPoint<FReal> localChildPoints;
        for (unsigned int n=0; n<ORDER; ++n) {
            map(FPoint<FReal>(globalChildCoords[0][n],globalChildCoords[1][n],globalChildCoords[2][n]), localChildPoints);
            localChildCoords[0][n] = localChildPoints.getX();
            localChildCoords[1][n] = localChildPoints.getY();
            localChildCoords[2][n] = localChildPoints.getZ();
        }

        // assemble interpolator
        FReal* subChildParentInterpolator = new FReal [3 * ORDER*ORDER];
        KernelBaseClass::Interpolator->assembleInterpolator(ORDER, localChildCoords[0], subChildParentInterpolator);
        KernelBaseClass::Interpolator->assembleInterpolator(ORDER, localChildCoords[1], subChildParentInterpolator + 1 * ORDER*ORDER);
        KernelBaseClass::Interpolator->assembleInterpolator(ORDER, localChildCoords[2], subChildParentInterpolator + 2 * ORDER*ORDER);

        // get permutation operators
        unsigned int perm[3][nnodes];
        for (unsigned int i=0;i<3; ++i)
            for (unsigned int n=0; n<nnodes; ++n)
                perm[i][n] = KernelBaseClass::Interpolator->getPermutationsM2ML2L(i)[n];

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

            // 1) apply Sy (using tensor product M2M with an interpolator computed on the fly)

            // Do NOT reset multipole expansion !
            //FBlas::scal(nnodes, FReal(0.), pole->getMultipole(idxRhs));

            //      /// p^6 version
            //      FBlas::gemtva(nnodes, nnodes, FReal(1.),
            //                    subChildParentInterpolator,
            //                    const_cast<FReal*>(subCell->getMultipole(idxRhs)), pole->getMultipole(idxRhs));

            /// p^4 version
            FReal Exp[nnodes], PermExp[nnodes];
            // ORDER*ORDER*ORDER * (2*ORDER-1)
            FBlas::gemtm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                         subChildParentInterpolator, ORDER,
                         const_cast<FReal*>(subCell->getMultipole(idxRhs)), ORDER, PermExp, ORDER);

            for (unsigned int n=0; n<nnodes; ++n)	Exp[n] = PermExp[perm[1][n]];
            // ORDER*ORDER*ORDER * (2*ORDER-1)
            FBlas::gemtm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                         subChildParentInterpolator + 2 * ORDER*ORDER, ORDER,
                         Exp, ORDER, PermExp, ORDER);

            for (unsigned int n=0; n<nnodes; ++n)	Exp[perm[1][n]] = PermExp[perm[2][n]];
            // ORDER*ORDER*ORDER * (2*ORDER-1)
            FBlas::gemtm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                         subChildParentInterpolator + 1 * ORDER*ORDER, ORDER,
                         Exp, ORDER, PermExp, ORDER);

            for (unsigned int n=0; n<nnodes; ++n)	pole->getMultipole(idxRhs)[perm[2][n]] += PermExp[n];



            //      // 2) Apply Discete Fourier Transform
            //      M2LHandler.applyZeroPaddingAndDFT(ParentCell->getMultipole(idxRhs),
            //                                         ParentCell->getTransformedMultipole(idxRhs));
        }
    }

    void P2L(CellClass* const local, const int localLevel, const ContainerClass* const particles) override {

        // Target cell: local
        const FReal localCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0, localLevel)));
        const FPoint<FReal> localCellCenter(KernelBaseClass::getCellCenter(local->getCoordinate(),localLevel));
        //   std::cout << "   call P2L  localLevel "<< localLevel << "  localCellCenter "<< localCellCenter <<std::endl;
        // interpolation points of target (X) cell
        FPoint<FReal> X[nnodes];
        FUnifTensor<FReal,order>::setRoots(localCellCenter, localCellWidth, X);

        // read positions
        const FReal*const positionsX = particles->getPositions()[0];
        const FReal*const positionsY = particles->getPositions()[1];
        const FReal*const positionsZ = particles->getPositions()[2];

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

            // read physicalValue
            const FReal*const physicalValues = particles->getPhysicalValues();

            // apply P2L
            for (FSize idxPart=0; idxPart<particles->getNbParticles(); ++idxPart){

                const FPoint<FReal> y = FPoint<FReal>(positionsX[idxPart],
                                                      positionsY[idxPart],
                                                      positionsZ[idxPart]);

                for (unsigned int m=0; m<nnodes; ++m)
                    local->getLocal(idxRhs)[m]+=MatrixKernel->evaluate(X[m], y) * physicalValues[idxPart];

            }

        }// NVALS

    }

    void M2L(CellClass* const local, const int localLevel, const CellClass* const pole, const int poleLevel) override {

        // Source cell: pole
        const FReal poleCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0, poleLevel)));
        const FPoint<FReal> poleCellCenter(KernelBaseClass::getCellCenter(pole->getCoordinate(),poleLevel));

        // Target cell: local
        const FReal localCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0, localLevel)));
        const FPoint<FReal> localCellCenter(KernelBaseClass::getCellCenter(local->getCoordinate(),localLevel));

        // interpolation points of source (Y) and target (X) cell
        FPoint<FReal> X[nnodes], Y[nnodes];
        FUnifTensor<FReal,order>::setRoots(poleCellCenter, poleCellWidth, Y);
        FUnifTensor<FReal,order>::setRoots(localCellCenter, localCellWidth, X);


        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

            // Dense M2L
            const FReal *const MultipoleExpansion = pole->getMultipole(idxRhs);

            for (unsigned int m=0; m<nnodes; ++m)
                for (unsigned int n=0; n<nnodes; ++n){
                    local->getLocal(idxRhs)[m]+=MatrixKernel->evaluate(X[m], Y[n]) * MultipoleExpansion[n];

                }
        }
    }

    void M2P(const CellClass* const pole, const int poleLevel, ContainerClass* const particles) override {

        // Source cell: pole
        const FReal poleCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0, poleLevel)));
        const FPoint<FReal> poleCellCenter(KernelBaseClass::getCellCenter(pole->getCoordinate(),poleLevel));

        // interpolation points of source (Y) cell
        FPoint<FReal> Y[nnodes];
        FUnifTensor<FReal,order>::setRoots(poleCellCenter, poleCellWidth, Y);

        // read positions
        const FReal*const positionsX = particles->getPositions()[0];
        const FReal*const positionsY = particles->getPositions()[1];
        const FReal*const positionsZ = particles->getPositions()[2];

        // get potential
        FReal*const physVal = particles->getPhysicalValues(/*idxPot*/);
        FReal*const potentials = particles->getPotentials(/*idxPot*/);
        FReal*const fx = particles->getForcesX(/*idxPot*/);
        FReal*const fy = particles->getForcesY(/*idxPot*/);
        FReal*const fz = particles->getForcesZ(/*idxPot*/);

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

            const FReal *const MultipoleExpansion = pole->getMultipole(idxRhs);

            // apply M2P
            for ( FSize idxPart=0; idxPart<particles->getNbParticles(); ++idxPart){

                const FPoint<FReal> x = FPoint<FReal>(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]);

                for (unsigned int n=0; n<nnodes; ++n){

                    FReal Kxy[1];
                    FReal dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(x,Y[n],Kxy,dKxy);

                    potentials[idxPart] += Kxy[0] * MultipoleExpansion[/*idxLhs*nnodes+*/n];
                    fx[idxPart] += -dKxy[0] * physVal[idxPart] * MultipoleExpansion[/*idxLhs*nnodes+*/n];
                    fy[idxPart] += -dKxy[1] * physVal[idxPart] * MultipoleExpansion[/*idxLhs*nnodes+*/n];
                    fz[idxPart] += -dKxy[2] * physVal[idxPart] * MultipoleExpansion[/*idxLhs*nnodes+*/n];

                }

            }// Particles

        }// NVALS

    }

    void L2L(const CellClass* const local, const int localLevel, CellClass* const subCell, const int subCellLevel) override {

        const FPoint<FReal> subCellCenter(KernelBaseClass::getCellCenter(subCell->getCoordinate(),subCellLevel));
        const FReal subCellWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0,subCellLevel)));

        const FPoint<FReal> localCenter(KernelBaseClass::getCellCenter(local->getCoordinate(),localLevel));
        const FReal localWidth(KernelBaseClass::BoxWidth / FReal(FMath::pow(2.0,localLevel)));

        //    ////////////////////////////////////////////////////////////////////////////
        //    /// p^6 version
        //    // allocate memory
        //    FReal* subChildParentInterpolator = new FReal [nnodes * nnodes];
        //
        //    // set child info
        //    FPoint<FReal> ChildRoots[nnodes], localChildRoots[nnodes];
        //    FUnifTensor<FReal,ORDER>::setRoots(subCellCenter, subCellWidth, ChildRoots);
        //
        //    // map global position of roots to local position in parent cell
        //    const map_glob_loc<FReal> map(localCenter, localWidth);
        //    for (unsigned int n=0; n<nnodes; ++n)
        //      map(ChildRoots[n], localChildRoots[n]);
        //
        //    // assemble child - parent - interpolator
        //    KernelBaseClass::Interpolator->assembleInterpolator(nnodes, localChildRoots, subChildParentInterpolator);

        ////////////////////////////////////////////////////////////////////////////
        /// p^4 version
        // Set sub-child coords
        FReal globalChildCoords[3][ORDER];
        FUnifTensor<FReal,order>::setPolynomialsRoots(subCellCenter, subCellWidth, globalChildCoords);

        // Map global position of sub-child nodes to [-1,1]
        FReal localChildCoords[3][ORDER];
        const map_glob_loc<FReal> map(localCenter, localWidth);
        FPoint<FReal> localChildPoints;
        for (unsigned int n=0; n<ORDER; ++n) {
            map(FPoint<FReal>(globalChildCoords[0][n],globalChildCoords[1][n],globalChildCoords[2][n]), localChildPoints);
            localChildCoords[0][n] = localChildPoints.getX();
            localChildCoords[1][n] = localChildPoints.getY();
            localChildCoords[2][n] = localChildPoints.getZ();
        }

        // assemble interpolator
        FReal* subChildParentInterpolator = new FReal [3 * ORDER*ORDER];
        KernelBaseClass::Interpolator->assembleInterpolator(ORDER, localChildCoords[0], subChildParentInterpolator);
        KernelBaseClass::Interpolator->assembleInterpolator(ORDER, localChildCoords[1], subChildParentInterpolator + 1 * ORDER*ORDER);
        KernelBaseClass::Interpolator->assembleInterpolator(ORDER, localChildCoords[2], subChildParentInterpolator + 2 * ORDER*ORDER);

        // get permutation operators
        unsigned int perm[3][nnodes];
        for (unsigned int i=0;i<3; ++i)
            for (unsigned int n=0; n<nnodes; ++n)
                perm[i][n] = KernelBaseClass::Interpolator->getPermutationsM2ML2L(i)[n];

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            //      // 1) Apply Inverse Discete Fourier Transform
            //      M2LHandler.unapplyZeroPaddingAndDFT(local->getTransformedLocal(idxRhs),
            //                                           const_cast<CellClass*>(local)->getLocal(idxRhs));

            // 2) apply Sx

            //      /// p^6 version
            //      FBlas::gemva(nnodes, nnodes, FReal(1.),
            //                   subChildParentInterpolator,
            //                   const_cast<FReal*>(local->getLocal(idxRhs)), subCell->getLocal(idxRhs));

            /// p^4 version
            FReal Exp[nnodes], PermExp[nnodes];
            // ORDER*ORDER*ORDER * (2*ORDER-1)
            FBlas::gemm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                        subChildParentInterpolator, ORDER,
                        const_cast<FReal*>(local->getLocal(idxRhs)), ORDER, PermExp, ORDER);

            for (unsigned int n=0; n<nnodes; ++n)	Exp[n] = PermExp[perm[1][n]];
            // ORDER*ORDER*ORDER * (2*ORDER-1)
            FBlas::gemm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                        subChildParentInterpolator + 2 * ORDER*ORDER, ORDER,
                        Exp, ORDER, PermExp, ORDER);

            for (unsigned int n=0; n<nnodes; ++n)	Exp[perm[1][n]] = PermExp[perm[2][n]];
            // ORDER*ORDER*ORDER * (2*ORDER-1)
            FBlas::gemm(ORDER, ORDER, ORDER*ORDER, FReal(1.),
                        subChildParentInterpolator + 1 * ORDER*ORDER, ORDER,
                        Exp, ORDER, PermExp, ORDER);

            for (unsigned int n=0; n<nnodes; ++n)	subCell->getLocal(idxRhs)[perm[2][n]] += PermExp[n];
            // total flops count: 3 * ORDER*ORDER*ORDER * (2*ORDER-1)

        }

    }

    void L2P(const CellClass* const local, const int cellLevel, ContainerClass* const particles)  override {

        const FPoint<FReal> CellCenter(KernelBaseClass::getCellCenter(local->getCoordinate(),cellLevel));
        const FReal LocalBoxWidth = KernelBaseClass::BoxWidth / FMath::pow(2.0,cellLevel);

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

            //      // 1)  Apply Inverse Discete Fourier Transform
            //      M2LHandler.unapplyZeroPaddingAndDFT(local->getTransformedLocal(idxRhs),
            //                                          const_cast<CellClass*>(local)->getLocal(idxRhs));

            // 2.a) apply Sx
            KernelBaseClass::Interpolator->applyL2P(CellCenter, LocalBoxWidth, local->getLocal(idxRhs), particles);

            // 2.b) apply Px (grad Sx)
            KernelBaseClass::Interpolator->applyL2PGradient(CellCenter, LocalBoxWidth, local->getLocal(idxRhs), particles);

        }
    }

    void P2P(ContainerClass* target, const ContainerClass* sources)  override {
        ContainerClass* sourcesArray[27] = { const_cast<ContainerClass*> (sources) };
        DirectInteractionComputer<FReal,MatrixKernelClass::NCMP, NVALS>::template P2PRemote(target,sourcesArray,1,MatrixKernel);
    }

    bool preferP2M(const ContainerClass* const particles) override {
        return particles->getNbParticles() >this->sminM;
    }
    bool preferP2M(const int /*atLevel*/, const ContainerClass*const particles[], const int nbContainers) override {
        FSize counterParticles = 0;
        for(FSize idxContainer = 0 ; idxContainer < nbContainers ; ++idxContainer){
            counterParticles += particles[idxContainer]->getNbParticles();
        }
        //		std::cout << " Part("<<counterParticles<< ") ";
        return counterParticles >this->sminM;
    }
};

//
//template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
//class FAdaptUnifKernel
//		: public FUnifKernel<FReal,CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
//{
//	typedef FUnifKernel<FReal,CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	KernelBaseClass;
//
//#ifdef LOG_TIMINGS
//	FTic time;
//	FReal t_m2l_1, t_m2l_2, t_m2l_3;
//#endif
//
//public:
//	/**
//	 * The constructor initializes all constant attributes and it reads the
//	 * precomputed and compressed M2L operators from a binary file (an
//	 * runtime_error is thrown if the required file is not valid).
//	 */
//	FAdaptUnifKernel(const int inTreeHeight,
//			const FReal inBoxWidth,
//			const FPoint<FReal>& inBoxCenter)
//: KernelBaseClass(inTreeHeight, inBoxWidth, inBoxCenter)
//{
//
//#ifdef LOG_TIMINGS
//		t_m2l_1 = FReal(0.);
//		t_m2l_2 = FReal(0.);
//		t_m2l_3 = FReal(0.);
//#endif
//}
//
//
//	/** Copy constructor */
//	FAdaptUnifKernel(const FAdaptUnifKernel& other)
//	: KernelBaseClass(other)
//	{	}
//
//
//
//	/** Destructor */
//	~FAdaptUnifKernel()
//	{
//		this->~KernelBaseClass() ;
//#ifdef LOG_TIMINGS
//		std::cout << "- Permutation took " << t_m2l_1 << "s"
//				<< "\n- GEMMT and GEMM took " << t_m2l_2 << "s"
//				<< "\n- Unpermutation took " << t_m2l_3 << "s"
//				<< std::endl;
//#endif
//	}
//
//
//	void P2MAdapt(CellClass* const ParentCell,  const int &level)
//	{
//		const FPoint<FReal> LeafCellCenter(KernelBaseClass::getLeafCellCenter(ParentCell->getCoordinate()));
//		const FReal BoxWidth = KernelBaseClass::BoxWidthLeaf*FMath::pow(2.0,KernelBaseClass::TreeHeight-level);
//		//
//		for(int i = 0 ; i <ParentCell->getLeavesSize(); ++i ){
//			//
//			for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
//				KernelBaseClass::Interpolator->applyP2M(LeafCellCenter, BoxWidth,
//						ParentCell->getMultipole(idxRhs), ParentCell->getLeaf(i)->getSrc());
//			}
//		}
//	}
//	void M2MAdapt(CellClass* const FRestrict ParentCell, const int &TreeLevel, const int &numberOfM2M,
//			const int * FRestrict ChildLevel , const CellClass*const FRestrict *const FRestrict ChildCells)
//	{
//		for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
//			//            // apply Sy
//			for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
//				if (ChildCells[ChildIndex]){
//					//			KernelBaseClass::Interpolator->applyM2M(ChildIndex, ChildCells[ChildIndex]->getMultipole(idxRhs), ParentCell->getMultipole(idxRhs));
//				}
//			}
//		}
//	}
//
//
//
//	void M2L(CellClass* const FRestrict TargetCell,
//			const CellClass* SourceCells[343],
//			const int /*NumSourceCells*/,
//			const int TreeLevel)
//	{
//
//	}
//
//
//	void L2L(const CellClass* const FRestrict ParentCell,
//			CellClass* FRestrict *const FRestrict ChildCells,
//			const int /*TreeLevel*/)
//	{
//		//        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
//		//            // apply Sx
//		//            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
//		//                if (ChildCells[ChildIndex]){
//		//                    AbstractBaseClass::Interpolator->applyL2L(ChildIndex, ParentCell->getLocal(idxRhs), ChildCells[ChildIndex]->getLocal(idxRhs));
//		//                }
//		//            }
//		//        }
//	}
//
//	void L2P(const CellClass* const LeafCell,
//			ContainerClass* const TargetParticles)
//	{
//		KernelBaseClass::L2P(LeafCell,TargetParticles) ;
//	}
//
//	//    void P2P(const FTreeCoordinate& /* LeafCellCoordinate */, // needed for periodic boundary conditions
//	//                     ContainerClass* const FRestrict TargetParticles,
//	//                     const ContainerClass* const FRestrict /*SourceParticles*/,
//	//                     ContainerClass* const NeighborSourceParticles[27],
//	//                     const int /* size */)
//	//    {
//	//        DirectInteractionComputer<MatrixKernelClass::Identifier, NVALS>::P2P(TargetParticles,NeighborSourceParticles);
//	//    }
//	//
//	//
//	//    void P2PRemote(const FTreeCoordinate& /*inPosition*/,
//	//                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
//	//                   ContainerClass* const inNeighbors[27], const int /*inSize*/){
//	//       DirectInteractionComputer<MatrixKernelClass::Identifier, NVALS>::P2PRemote(inTargets,inNeighbors,27);
//	//    }
//
//};
//
//






#endif //FADAPTUNIFKERNELS_HPP

// [--END--]
