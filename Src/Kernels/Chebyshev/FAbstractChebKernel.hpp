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

#include "./FChebInterpolator.hpp"

class FTreeCoordinate;
template <KERNEL_FUNCCTION_IDENTIFIER Identifier> struct DirectInteactionComputer;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FAbstractChebKernel
 * @brief
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
					 const ContainerClass* const FRestrict SourceParticles,
					 ContainerClass* const NeighborSourceParticles[27],
					 const int /* size */)
	{
		// loop: target particles
		typename ContainerClass::BasicIterator iTargets(*TargetParticles);

		if (TargetParticles != SourceParticles) {

			while (iTargets.hasNotFinished()) {
				ParticleClass& Target = iTargets.data();
				
				{ // loop: source particles (target leaf cell == source leaf cell)
					typename ContainerClass::ConstBasicIterator iSources(*SourceParticles);
					while (iSources.hasNotFinished()) {
						const ParticleClass& Source = iSources.data();
						// only if target and source are not identical
						DirectInteactionComputer<MatrixKernelClass::Identifier>::compute(Target, Source);
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
								DirectInteactionComputer<MatrixKernelClass::Identifier>::compute(Target, Source);
								// progress sources
								iSources.gotoNext();
							}
						}
					}
				}
		
				// progress targets
				iTargets.gotoNext();
			}

		} else {

			while (iTargets.hasNotFinished()) {
				ParticleClass& Target = iTargets.data();
				
				{ // loop: source particles  (target leaf cell == source leaf cell)
					typename ContainerClass::BasicIterator iSources = iTargets;
					iSources.gotoNext();
					while (iSources.hasNotFinished()) {
						ParticleClass& Source = iSources.data();
						// only if target and source are not identical
						DirectInteactionComputer<MatrixKernelClass::Identifier>::computeMutual(Target, Source);
						// progress sources
						iSources.gotoNext();
					}
				}
				
				{ // loop: source particles (target leaf cell != source leaf cell)
					for (unsigned int idx=0; idx<=13; ++idx) {
						if (NeighborSourceParticles[idx]) {
							typename ContainerClass::BasicIterator iSources(*NeighborSourceParticles[idx]);
							while (iSources.hasNotFinished()) {
								ParticleClass& Source = iSources.data();
								// target and source cannot be identical
								DirectInteactionComputer<MatrixKernelClass::Identifier>::computeMutual(Target, Source);
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

	}

};


/*! Specialization for Laplace potential */
template <>
struct DirectInteactionComputer<ONE_OVER_R>
{
	template <typename ParticleClass>
	static void compute(ParticleClass& Target, const ParticleClass& Source) // 34 overall flops
	{
		FPoint xy(Source.getPosition() - Target.getPosition()); // 3 flops
		const FReal one_over_r = FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
																										 xy.getY()*xy.getY() +
																										 xy.getZ()*xy.getZ()); // 1 + 15 + 5 = 21 flops
		const FReal wt = Target.getPhysicalValue();
		const FReal ws = Source.getPhysicalValue();

		// laplace potential
		Target.incPotential(one_over_r * ws); // 2 flops

		// force
		xy *= ((ws*wt) * (one_over_r*one_over_r*one_over_r)); // 5 flops
		Target.incForces(xy.getX(), xy.getY(), xy.getZ()); // 3 flops
	}

	template <typename ParticleClass>
	static void computeMutual(ParticleClass& Target, ParticleClass& Source) // 39 overall flops
	{
		FPoint xy(Source.getPosition() - Target.getPosition()); // 3 flops
		const FReal one_over_r = FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
																										 xy.getY()*xy.getY() +
																										 xy.getZ()*xy.getZ()); // 1 + 15 + 5 = 21 flops
		const FReal wt = Target.getPhysicalValue();
		const FReal ws = Source.getPhysicalValue();

		// laplace potential
		Target.incPotential(one_over_r * ws); // 2 flops
		Source.incPotential(one_over_r * wt); // 2 flops

		// force
		xy *= ((ws*wt) * (one_over_r*one_over_r*one_over_r)); // 5 flops
		Target.incForces(  xy.getX(),    xy.getY(),    xy.getZ());  // 3 flops 
		Source.incForces((-xy.getX()), (-xy.getY()), (-xy.getZ())); // 3 flops
	}
};


/*! Specialization for Leonard-Jones potential */
template <>
struct DirectInteactionComputer<LEONARD_JONES_POTENTIAL>
{
	template <typename ParticleClass>
	static void compute(ParticleClass& Target, const ParticleClass& Source) // 39 overall flops
	{
		FPoint xy(Source.getPosition() - Target.getPosition()); // 3 flops
		const FReal one_over_r = FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
																										 xy.getY()*xy.getY() +
																										 xy.getZ()*xy.getZ()); // 1 + 15 + 5 = 21 flops
		const FReal wt = Target.getPhysicalValue();
		const FReal ws = Source.getPhysicalValue();

		// lenard-jones potential
		const FReal one_over_r3 = one_over_r * one_over_r * one_over_r;
		const FReal one_over_r6 = one_over_r3 * one_over_r3;
		Target.incPotential((one_over_r6*one_over_r6 - one_over_r6) * ws); // 2 flops

		// force
		const FReal one_over_r4 = one_over_r3 * one_over_r; // 1 flop
		xy *= ((ws*wt) * (FReal(12.)*one_over_r6*one_over_r4*one_over_r4 - FReal(6.)*one_over_r4*one_over_r4)); // 9 flops
		Target.incForces(xy.getX(), xy.getY(), xy.getZ()); // 3 flops
	}

	template <typename ParticleClass>
	static void computeMutual(ParticleClass& Target, ParticleClass& Source) // 44 overall flops
	{
		FPoint xy(Source.getPosition() - Target.getPosition()); // 3 flops
		const FReal one_over_r = FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
																										 xy.getY()*xy.getY() +
																										 xy.getZ()*xy.getZ()); // 1 + 15 + 5 = 21 flops
		const FReal wt = Target.getPhysicalValue();
		const FReal ws = Source.getPhysicalValue();

		// lenard-jones potential
		const FReal one_over_r3 = one_over_r * one_over_r * one_over_r;
		const FReal one_over_r6 = one_over_r3 * one_over_r3;
		Target.incPotential((one_over_r6*one_over_r6 - one_over_r6) * ws); // 2 flops
		Source.incPotential((one_over_r6*one_over_r6 - one_over_r6) * wt); // 2 flops

		// force
		const FReal one_over_r4 = one_over_r3 * one_over_r; // 1 flop
		xy *= ((ws*wt) * (FReal(12.)*one_over_r6*one_over_r4*one_over_r4 - FReal(6.)*one_over_r4*one_over_r4)); // 9 flops
		Target.incForces(  xy.getX(),    xy.getY(),    xy.getZ());  // 3 flops 
		Source.incForces((-xy.getX()), (-xy.getY()), (-xy.getZ())); // 3 flops
	}
};



#endif //FCHEBKERNELS_HPP

// [--END--]
