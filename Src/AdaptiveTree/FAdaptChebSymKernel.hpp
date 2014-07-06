#ifndef FADAPTCHEBSYMKERNEL_HPP
#define FADAPTCHEBSYMKERNEL_HPP
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
#include "Utils/FTrace.hpp"
#include "Utils/FPoint.hpp"

#include "Adaptative/FAdaptiveCell.hpp"
#include "Adaptative/FAdaptiveKernelWrapper.hpp"
#include "Adaptative/FAbstractAdaptiveKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

class FTreeCoordinate;

// ==== CMAKE =====
// @FUSE_BLAS
// ================

// for verbosity only!!!
//#define COUNT_BLOCKED_INTERACTIONS

// if timings should be logged
//#define LOG_TIMINGS

/**
 * @author O. Coulaud
 * @class FAdaptChebSymKernel
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

template< class CellClass, class ContainerClass, class MatrixKernelClass, int ORDER, int NVALS = 1>
class FAdaptiveChebSymKernel : FChebSymKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
, public FAbstractAdaptiveKernel<CellClass, ContainerClass> {
	//
	typedef FChebSymKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	KernelBaseClass;
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
	FAdaptiveChebSymKernel(const int inTreeHeight, const FReal inBoxWidth,
			const FPoint& inBoxCenter) : KernelBaseClass(inTreeHeight, inBoxWidth, inBoxCenter)
	{}
	//	/** Copy constructor */
	FAdaptiveChebSymKernel(const FAdaptiveChebSymKernel& other)
		: KernelBaseClass(other)
		{	}

	//
	//	/** Destructor */
		~FAdaptiveChebSymKernel()
		{
			//this->~KernelBaseClass() ;
		}
	void P2M(CellClass* const pole, const int cellLevel, const ContainerClass* const particles) override {
		//pole->setDataUp(pole->getDataUp() + particles->getNbParticles());
		//
		const FPoint LeafCellCenter(KernelBaseClass::getLeafCellCenter(pole->getCoordinate()));
		const FReal BoxWidth = KernelBaseClass::BoxWidthLeaf*FMath::pow(2.0,KernelBaseClass::TreeHeight-cellLevel);
		//

		for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
			KernelBaseClass::Interpolator->applyP2M(LeafCellCenter, BoxWidth,
					pole->getMultipole(idxRhs), particles);
		}

	}

	void M2M(CellClass* const pole, const int /*poleLevel*/, const CellClass* const subCell, const int /*subCellLevel*/) override {
	//	pole->setDataUp(pole->getDataUp() + subCell->getDataUp());
	}

	void P2L(CellClass* const local, const int /*localLevel*/, const ContainerClass* const particles) override {
	//	local->setDataDown(local->getDataDown() + particles->getNbParticles());
	}

	void M2L(CellClass* const local, const int /*localLevel*/, const CellClass* const aNeighbor, const int /*neighborLevel*/) override {
	//	local->setDataDown(local->getDataDown() + aNeighbor->getDataUp());
	}

	void M2P(const CellClass* const pole, const int /*poleLevel*/, ContainerClass* const particles) override {
//		long long int*const particlesAttributes = particles->getDataDown();
//		for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
//			particlesAttributes[idxPart] += pole->getDataUp();
//		}
	}

	void L2L(const CellClass* const local, const int /*localLevel*/, CellClass* const subCell, const int /*subCellLevel*/) override {
	//	subCell->setDataDown(local->getDataDown() + subCell->getDataDown());
	}

	void L2P(const CellClass* const local, const int /*cellLevel*/, ContainerClass* const particles)  override {
//		long long int*const particlesAttributes = particles->getDataDown();
//		for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
//			particlesAttributes[idxPart] += local->getDataDown();
//		}
	}

	void P2P(ContainerClass* target, const ContainerClass* sources)  override {
//		long long int*const particlesAttributes = target->getDataDown();
//		for(int idxPart = 0 ; idxPart < target->getNbParticles() ; ++idxPart){
//			particlesAttributes[idxPart] += sources->getNbParticles();
//		}
	}

	bool preferP2M(const ContainerClass* const particles) override {
		return particles->getNbParticles() < 10;
	}
	bool preferP2M(const int /*atLevel*/, const ContainerClass*const particles[], const int nbContainers) override {
		int counterParticles = 0;
		for(int idxContainer = 0 ; idxContainer < nbContainers ; ++idxContainer){
			counterParticles += particles[idxContainer]->getNbParticles();
		}
		return counterParticles < 10;
	}
};

//
//template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
//class FAdaptChebSymKernel
//		: public FChebSymKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
//{
//	typedef FChebSymKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	KernelBaseClass;
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
//	FAdaptChebSymKernel(const int inTreeHeight,
//			const FReal inBoxWidth,
//			const FPoint& inBoxCenter)
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
//	FAdaptChebSymKernel(const FAdaptChebSymKernel& other)
//	: KernelBaseClass(other)
//	{	}
//
//
//
//	/** Destructor */
//	~FAdaptChebSymKernel()
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
//		const FPoint LeafCellCenter(KernelBaseClass::getLeafCellCenter(ParentCell->getCoordinate()));
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






#endif //FADAPTCHEBSYMKERNELS_HPP

// [--END--]
