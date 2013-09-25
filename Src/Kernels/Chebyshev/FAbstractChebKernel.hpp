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
#ifndef FABSTRACTCHEBKERNEL_HPP
#define FABSTRACTCHEBKERNEL_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FSmartPointer.hpp"

#include "../../Components/FAbstractKernels.hpp"

#include "../P2P/FP2P.hpp"

#include "./FChebInterpolator.hpp"

class FTreeCoordinate;
template <KERNEL_FUNCCTION_IDENTIFIER Identifier, int NVALS> struct DirectInteactionComputer;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FAbstractChebKernel
 * @brief
 * This kernels implement the Chebyshev interpolation based FMM operators. It
 * implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template < class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER, int NVALS = 1>
class FAbstractChebKernel : public FAbstractKernels< CellClass, ContainerClass>
{
protected:
  enum {nnodes = TensorTraits<ORDER>::nnodes};
	typedef FChebInterpolator<ORDER> InterpolatorClass;

	/// Needed for P2M, M2M, L2L and L2P operators
	const FSmartPointer<InterpolatorClass,FSmartPointerMemory> Interpolator;
	/// Needed for P2P operator
	const FSmartPointer<MatrixKernelClass,FSmartPointerMemory> MatrixKernel;
	/// Height of the entire oct-tree
	const unsigned int TreeHeight;
	/// Corner of oct-tree box
	const FPoint BoxCorner;
	/// Width of oct-tree box   
	const FReal BoxWidth;    
	/// Width of a leaf cell box 
	const FReal BoxWidthLeaf;

	/**
	 * Compute center of leaf cell from its tree coordinate.
	 * @param[in] Coordinate tree coordinate
	 * @return center of leaf cell
	 */
	const FPoint getLeafCellCenter(const FTreeCoordinate& Coordinate) const
	{
		return FPoint(BoxCorner.getX() + (FReal(Coordinate.getX()) + FReal(.5)) * BoxWidthLeaf,
											 BoxCorner.getY() + (FReal(Coordinate.getY()) + FReal(.5)) * BoxWidthLeaf,
											 BoxCorner.getZ() + (FReal(Coordinate.getZ()) + FReal(.5)) * BoxWidthLeaf);
	}

public:
	/**
	 * The constructor initializes all constant attributes and it reads the
	 * precomputed and compressed M2L operators from a binary file (an
	 * runtime_error is thrown if the required file is not valid).
	 */
	FAbstractChebKernel(const int inTreeHeight,
											const FPoint& inBoxCenter,
											const FReal inBoxWidth)
		: Interpolator(new InterpolatorClass()),
			MatrixKernel(new MatrixKernelClass()),
			TreeHeight(inTreeHeight),
			BoxCorner(inBoxCenter - inBoxWidth / FReal(2.)),
			BoxWidth(inBoxWidth),
			BoxWidthLeaf(BoxWidth / FReal(FMath::pow(2, inTreeHeight - 1)))
	{
		/* empty */
	}

    virtual ~FAbstractChebKernel(){
        // should not be used
    }

	const InterpolatorClass *const getPtrToInterpolator() const
	{ return Interpolator.getPtr(); }


	virtual void P2M(CellClass* const LeafCell,
									 const ContainerClass* const SourceParticles) = 0;


	virtual void M2M(CellClass* const FRestrict ParentCell,
									 const CellClass*const FRestrict *const FRestrict ChildCells,
									 const int TreeLevel) = 0;


	virtual void M2L(CellClass* const FRestrict TargetCell,
									 const CellClass* SourceCells[343],
									 const int NumSourceCells,
									 const int TreeLevel) = 0;


	virtual void L2L(const CellClass* const FRestrict ParentCell,
									 CellClass* FRestrict *const FRestrict ChildCells,
									 const int TreeLevel) = 0;


	virtual void L2P(const CellClass* const LeafCell,
									 ContainerClass* const TargetParticles) = 0;
	
	

	void P2P(const FTreeCoordinate& /* LeafCellCoordinate */, // needed for periodic boundary conditions
					 ContainerClass* const FRestrict TargetParticles,
                     const ContainerClass* const FRestrict /*SourceParticles*/,
					 ContainerClass* const NeighborSourceParticles[27],
					 const int /* size */)
	{
        DirectInteactionComputer<MatrixKernelClass::Identifier, NVALS>::P2P(TargetParticles,NeighborSourceParticles);
	}


    void P2PRemote(const FTreeCoordinate& /*inPosition*/,
                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
                   ContainerClass* const inNeighbors[27], const int /*inSize*/){
        DirectInteactionComputer<MatrixKernelClass::Identifier, NVALS>::P2PRemote(inTargets,inNeighbors,27);
    }

};


///////////////////////////////////////////////////////
// P2P Wrappers
///////////////////////////////////////////////////////

/*! Specialization for Laplace potential */
template <>
struct DirectInteactionComputer<ONE_OVER_R, 1>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14);
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        FP2P::FullRemote(inTargets,inNeighbors,inSize);
    }
};


/*! Specialization for Leonard-Jones potential */
template <>
struct DirectInteactionComputer<LEONARD_JONES_POTENTIAL, 1>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        FP2P::FullMutualLJ(TargetParticles,NeighborSourceParticles,14);
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        FP2P::FullRemoteLJ(inTargets,inNeighbors,inSize);
    }
};

///////////////////////////////////////////////////////
// In case of multi right hand side
///////////////////////////////////////////////////////

template <int NVALS>
struct DirectInteactionComputer<ONE_OVER_R, NVALS>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullMutual(TargetParticles,NeighborSourceParticles,14);
        }
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullRemote(inTargets,inNeighbors,inSize);
        }
    }
};


/*! Specialization for Leonard-Jones potential */
template <int NVALS>
struct DirectInteactionComputer<LEONARD_JONES_POTENTIAL, NVALS>
{
    template <typename ContainerClass>
    static void P2P(		 ContainerClass* const FRestrict TargetParticles,
                     ContainerClass* const NeighborSourceParticles[27]){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullMutualLJ(TargetParticles,NeighborSourceParticles,14);
        }
    }

    template <typename ContainerClass>
    static void P2PRemote( ContainerClass* const FRestrict inTargets,
                           ContainerClass* const inNeighbors[27],
                           const int inSize){
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FP2P::FullRemoteLJ(inTargets,inNeighbors,inSize);
        }
    }
};



#endif //FCHEBKERNELS_HPP

// [--END--]
