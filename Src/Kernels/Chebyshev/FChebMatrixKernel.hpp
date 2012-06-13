#ifndef FCHEBMATRIXKERNEL_HPP
#define FCHEBMATRIXKERNEL_HPP
// [--License--]

#include "../../Utils/FPoint.hpp"
#include "../../Utils/FNoCopyable.hpp"
#include "../../Utils/FMath.hpp"

// extendable
enum KERNEL_FUNCCTION_IDENTIFIER {ONE_OVER_R,
																	ONE_OVER_R_SQUARED,
																	LEONARD_JONES_POTENTIAL};

// probably not extedable :)
enum KERNEL_FUNCTION_TYPE {HOMOGENEOUS, NON_HOMOGENEOUS};


/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @class FChebMatrixKernels
 * Please read the license
 */
struct FChebAbstractMatrixKernel : FNoCopyable
{
	virtual FReal evaluate(const FPoint&, const FPoint&) const = 0;
	// I need both functions because required arguments are not always given
	virtual FReal getScaleFactor(const FReal, const int) const = 0;
	virtual FReal getScaleFactor(const FReal) const = 0;
};



/// One over r
struct FChebMatrixKernelR : FChebAbstractMatrixKernel
{
	static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
	static const KERNEL_FUNCCTION_IDENTIFIER Identifier = ONE_OVER_R;

	FChebMatrixKernelR() {}

	FReal evaluate(const FPoint& x, const FPoint& y) const
	{
		const FPoint xy(x-y);
		return FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
																	 xy.getY()*xy.getY() +
																	 xy.getZ()*xy.getZ());
	}

	FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
	{
		const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
		return getScaleFactor(CellWidth);
	}

	FReal getScaleFactor(const FReal CellWidth) const
	{
		return FReal(2.) / CellWidth;
	}
};



/// One over r^2
struct FChebMatrixKernelRR : FChebAbstractMatrixKernel
{
	static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
	static const KERNEL_FUNCCTION_IDENTIFIER Identifier = ONE_OVER_R_SQUARED;

	FChebMatrixKernelRR() {}

	FReal evaluate(const FPoint& x, const FPoint& y) const
	{
		const FPoint xy(x-y);
		return FReal(1.) / FReal(xy.getX()*xy.getX() +
														 xy.getY()*xy.getY() +
														 xy.getZ()*xy.getZ());
	}

	FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
	{
		const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
		return getScaleFactor(CellWidth);
	}

	FReal getScaleFactor(const FReal CellWidth) const
	{
		return FReal(4.) / CellWidth;
	}
};



/// One over r^12 - One over r^6
struct FChebMatrixKernelLJ : FChebAbstractMatrixKernel
{
	static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
	static const KERNEL_FUNCCTION_IDENTIFIER Identifier = LEONARD_JONES_POTENTIAL;

	FChebMatrixKernelLJ() {}

	FReal evaluate(const FPoint& x, const FPoint& y) const
	{
		const FPoint xy(x-y);
		const FReal r = xy.norm();
		const FReal r3 = r*r*r;
		const FReal one_over_r6 = FReal(1.) / (r3*r3);
		//return one_over_r6 * one_over_r6;
		//return one_over_r6;
		return one_over_r6 * one_over_r6 - one_over_r6;
	}

	FReal getScaleFactor(const FReal, const int) const
	{
		// return 1 because non homogeneous kernel functions cannot be scaled!!!
		return FReal(1.);
	}

	FReal getScaleFactor(const FReal) const
	{
		// return 1 because non homogeneous kernel functions cannot be scaled!!!
		return FReal(1.);
	}
};





#endif // FCHEBMATRIXKERNEL_HPP

// [--END--]
