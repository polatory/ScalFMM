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
#ifndef FCHEBMATRIXKERNEL_HPP
#define FCHEBMATRIXKERNEL_HPP

#include "../../Utils/FPoint.hpp"
#include "../../Utils/FNoCopyable.hpp"
#include "../../Utils/FMath.hpp"
#include "../../Utils/FBlas.hpp"


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
    virtual ~FChebAbstractMatrixKernel(){} // to remove warning
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





/*!  Functor which provides the interface to assemble a matrix based on the
  number of rows and cols and on the coordinates x and y and the type of the
  generating matrix-kernel function.
 */
template <typename MatrixKernel>
class EntryComputer
{
	const MatrixKernel Kernel;

	const unsigned int nx, ny;
	const FPoint *const px, *const py;

	const FReal *const weights;

public:
	explicit EntryComputer(const unsigned int _nx, const FPoint *const _px,
												 const unsigned int _ny, const FPoint *const _py,
												 const FReal *const _weights = NULL)
		: Kernel(),	nx(_nx), ny(_ny), px(_px), py(_py), weights(_weights) {}

//	template <typename Point>
//	void operator()(const unsigned int nx, const Point *const px,
//									const unsigned int ny, const Point *const py,
//									FReal *const data) const
//	{
//		for (unsigned int j=0; j<ny; ++j)
//			for (unsigned int i=0; i<nx; ++i)
//				data[j*nx + i] = Kernel.evaluate(px[i], py[j]);
//	}

	void operator()(const unsigned int xbeg, const unsigned int xend,
									const unsigned int ybeg, const unsigned int yend,
									FReal *const data) const
	{
		unsigned int idx = 0;
		if (weights) {
			for (unsigned int j=ybeg; j<yend; ++j)
				for (unsigned int i=xbeg; i<xend; ++i)
					data[idx++] = weights[i] * weights[j] * Kernel.evaluate(px[i], py[j]);
		} else {
			for (unsigned int j=ybeg; j<yend; ++j)
				for (unsigned int i=xbeg; i<xend; ++i)
					data[idx++] = Kernel.evaluate(px[i], py[j]);
		}

		/*
		// apply weighting matrices
		if (weights) {
			if ((xend-xbeg) == (yend-ybeg) && (xend-xbeg) == nx)
				for (unsigned int n=0; n<nx; ++n) {
					FBlas::scal(nx, weights[n], data + n,  nx); // scale rows
					FBlas::scal(nx, weights[n], data + n * nx); // scale cols
				}
			else if ((xend-xbeg) == 1 && (yend-ybeg) == ny)
				for (unsigned int j=0; j<ny; ++j)	data[j] *= weights[j];
			else if ((yend-ybeg) == 1 && (xend-xbeg) == nx)
				for (unsigned int i=0; i<nx; ++i)	data[i] *= weights[i];
		}
		*/

	}
};





#endif // FCHEBMATRIXKERNEL_HPP

// [--END--]
