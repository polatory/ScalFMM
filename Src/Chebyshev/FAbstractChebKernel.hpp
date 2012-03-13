#ifndef FABSTRACTCHEBKERNEL_HPP
#define FABSTRACTCHEBKERNEL_HPP
// [--License--]

#include "../Utils/FGlobal.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FSmartPointer.hpp"

#include "../Components/FAbstractKernels.hpp"

#include "./FChebInterpolator.hpp"

class FTreeCoordinate;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FAbstractChebKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Chebyshev interpolation based FMM operators. It
 * implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam ParticleClass Type of particle
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template <class ParticleClass, class CellClass,	class ContainerClass,	class MatrixKernelClass, int ORDER>
class FAbstractChebKernel : public FAbstractKernels<ParticleClass, CellClass, ContainerClass>
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
	const F3DPosition BoxCorner;
	/// Width of oct-tree box   
	const FReal BoxWidth;    
	/// Width of a leaf cell box 
	const FReal BoxWidthLeaf;

	/**
	 * Compute center of leaf cell from its tree coordinate.
	 * @param[in] Coordinate tree coordinate
	 * @return center of leaf cell
	 */
	const F3DPosition getLeafCellCenter(const FTreeCoordinate& Coordinate) const
	{
		return F3DPosition(BoxCorner.getX() + (FReal(Coordinate.getX()) + FReal(.5)) * BoxWidthLeaf,
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
											const F3DPosition& inBoxCenter,
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
	
	

	void P2P(const FTreeCoordinate& LeafCellCoordinate,
					 ContainerClass* const FRestrict TargetParticles,
					 const ContainerClass* const FRestrict SourceParticles,
					 ContainerClass* const NeighborSourceParticles[27],
					 const int )
	{
		// loop: target particles
		typename ContainerClass::BasicIterator iTargets(*TargetParticles);
		while (iTargets.hasNotFinished()) {
			ParticleClass& Target = iTargets.data();
			const FReal wt = Target.getPhysicalValue();

			{ // loop: source particles (target leaf cell == source leaf cell)
				typename ContainerClass::ConstBasicIterator iSources(*SourceParticles);
				while (iSources.hasNotFinished()) {
					const ParticleClass& Source = iSources.data();
					// only if target and source are not identical
					if (&Target != &Source) {
						const FReal one_over_r = MatrixKernel->evaluate(Target.getPosition(), Source.getPosition());
						const FReal ws = Source.getPhysicalValue();
						// potential
						Target.incPotential(one_over_r * ws);
						// force
						F3DPosition force(Target.getPosition() - Source.getPosition());
						force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
						Target.incForces(force.getX(), force.getY(), force.getZ());
					}
					// progress sources
					iSources.gotoNext();
				}
			}
			
			{ // loop: source particles (target leaf cell != source leaf cell)
				for (unsigned int idx=0; idx<27; ++idx) {
					if (NeighborSourceParticles[idx]) {
						typename ContainerClass::ConstBasicIterator	iSources(*NeighborSourceParticles[idx]);
						while (iSources.hasNotFinished()) {
							const ParticleClass& Source = iSources.data();
							// target and source cannot be identical
							const FReal one_over_r = MatrixKernel->evaluate(Target.getPosition(), Source.getPosition());
							const FReal ws = Source.getPhysicalValue();
							// potential
							Target.incPotential(one_over_r * ws);
							F3DPosition force(Target.getPosition() - Source.getPosition());
							force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
							// force
							Target.incForces(force.getX(), force.getY(), force.getZ());
							// progress sources
							iSources.gotoNext();
						}
					}
				}
			}
			
			// progress targets
			iTargets.gotoNext();
		}
	}

};


#endif //FCHEBKERNELS_HPP

// [--END--]
