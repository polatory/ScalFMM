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

#include "Kernels/Chebyshev/FChebSymKernel.hpp"

class FTreeCoordinate;


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
template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
class FAdaptChebSymKernel
		: public FChebSymKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
	typedef FChebSymKernel<CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>	KernelBaseClass;

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
	FAdaptChebSymKernel(const int inTreeHeight,
			const FReal inBoxWidth,
			const FPoint& inBoxCenter)
: KernelBaseClass(inTreeHeight, inBoxWidth, inBoxCenter)
{

#ifdef LOG_TIMINGS
		t_m2l_1 = FReal(0.);
		t_m2l_2 = FReal(0.);
		t_m2l_3 = FReal(0.);
#endif
}


	/** Copy constructor */
	FAdaptChebSymKernel(const FAdaptChebSymKernel& other)
	: KernelBaseClass(other)
	{	}



	/** Destructor */
	~FAdaptChebSymKernel()
	{
		this->~KernelBaseClass() ;
#ifdef LOG_TIMINGS
		std::cout << "- Permutation took " << t_m2l_1 << "s"
				<< "\n- GEMMT and GEMM took " << t_m2l_2 << "s"
				<< "\n- Unpermutation took " << t_m2l_3 << "s"
				<< std::endl;
#endif
	}


	void P2MAdapt(CellClass* const ParentCell,  const int &level)
	{
		const FPoint LeafCellCenter(KernelBaseClass::getLeafCellCenter(ParentCell->getCoordinate()));
		const FReal BoxWidth = KernelBaseClass::BoxWidthLeaf*FMath::pow(2.0,KernelBaseClass::TreeHeight-level);
		//
		for(int i = 0 ; i <ParentCell->getLeavesSize(); ++i ){
			//
			for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
				KernelBaseClass::Interpolator->applyP2M(LeafCellCenter, BoxWidth,
						ParentCell->getMultipole(idxRhs), ParentCell->getLeaf(i)->getSrc());
			}
		}
	}
	void M2MAdapt(CellClass* const FRestrict ParentCell, const int &TreeLevel, const int &numberOfM2M,
			const int * FRestrict ChildLevel , const CellClass*const FRestrict *const FRestrict ChildCells)
	{
		for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
			//            // apply Sy
			//            FBlas::scal(nnodes*2, FReal(0.), ParentCell->getMultipole(idxRhs));
			//            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
			//                if (ChildCells[ChildIndex]){
			//                    AbstractBaseClass::Interpolator->applyM2M(ChildIndex, ChildCells[ChildIndex]->getMultipole(idxRhs), ParentCell->getMultipole(idxRhs));
			//                }
			//            }
		}
	}



	void M2L(CellClass* const FRestrict TargetCell,
			const CellClass* SourceCells[343],
			const int /*NumSourceCells*/,
			const int TreeLevel)
	{

	}


	void L2L(const CellClass* const FRestrict ParentCell,
			CellClass* FRestrict *const FRestrict ChildCells,
			const int /*TreeLevel*/)
	{
		//        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
		//            // apply Sx
		//            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
		//                if (ChildCells[ChildIndex]){
		//                    AbstractBaseClass::Interpolator->applyL2L(ChildIndex, ParentCell->getLocal(idxRhs), ChildCells[ChildIndex]->getLocal(idxRhs));
		//                }
		//            }
		//        }
	}

	void L2P(const CellClass* const LeafCell,
			ContainerClass* const TargetParticles)
	{
		KernelBaseClass::L2P(LeafCell,TargetParticles) ;
	}

	//    void P2P(const FTreeCoordinate& /* LeafCellCoordinate */, // needed for periodic boundary conditions
	//                     ContainerClass* const FRestrict TargetParticles,
	//                     const ContainerClass* const FRestrict /*SourceParticles*/,
	//                     ContainerClass* const NeighborSourceParticles[27],
	//                     const int /* size */)
	//    {
	//        DirectInteractionComputer<MatrixKernelClass::Identifier, NVALS>::P2P(TargetParticles,NeighborSourceParticles);
	//    }
	//
	//
	//    void P2PRemote(const FTreeCoordinate& /*inPosition*/,
	//                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
	//                   ContainerClass* const inNeighbors[27], const int /*inSize*/){
	//       DirectInteractionComputer<MatrixKernelClass::Identifier, NVALS>::P2PRemote(inTargets,inNeighbors,27);
	//    }

};








#endif //FADAPTCHEBSYMKERNELS_HPP

// [--END--]
