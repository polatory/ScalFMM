#ifndef FCHEBMATRIXKERNEL_HPP
#define FCHEBMATRIXKERNEL_HPP
// [--License--]

#include "../Utils/F3DPosition.hpp"
#include "../Utils/FNoCopyable.hpp"
#include "../Utils/FMath.hpp"


enum {ONE_OVER_R         = 1,
			ONE_OVER_R_SQUARED = 2};


/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @class FChebMatrixKernels
 * Please read the license
 */
struct FChebAbstractMatrixKernel : FNoCopyable
{
	virtual FReal evaluate(const F3DPosition&, const F3DPosition&) const = 0;
	virtual FReal getScaleFactor(FReal) const = 0;
};

/// One over r
struct FChebMatrixKernelR : FChebAbstractMatrixKernel
{
	enum {Identifier = ONE_OVER_R};
	FChebMatrixKernelR() {}
	FReal evaluate(const F3DPosition& x, const F3DPosition& y) const
	{
		const F3DPosition xy(x-y);
		return FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
																	 xy.getY()*xy.getY() +
																	 xy.getZ()*xy.getZ());
	}

	FReal getScaleFactor(FReal CellWidth) const
	{ return FReal(2.) / CellWidth; }
};

/// One over r^2
struct FChebMatrixKernelRR : FChebAbstractMatrixKernel
{
	enum {Identifier = ONE_OVER_R_SQUARED};
	FChebMatrixKernelRR() {}
	FReal evaluate(const F3DPosition& x, const F3DPosition& y) const
	{
		const F3DPosition xy(x-y);
		return FReal(1.) / FReal(xy.getX()*xy.getX() +
														 xy.getY()*xy.getY() +
														 xy.getZ()*xy.getZ());
	}

	FReal getScaleFactor(FReal CellWidth) const
	{ return FReal(4.) / (CellWidth*CellWidth); }
};


#endif // FCHEBMATRIXKERNEL_HPP

// [--END--]
