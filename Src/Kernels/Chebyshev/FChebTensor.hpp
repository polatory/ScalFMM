#ifndef FCHEBTENSOR_HPP
#define FCHEBTENSOR_HPP

#include "../../Utils/FMath.hpp"

#include "./FChebRoots.hpp"
#include "./FChebMapping.hpp"


/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */

/**
 * @class TensorTraits
 *
 * The class @p TensorTraits gives the number of interpolation nodes per
 * cluster in 3D, depending on the interpolation order.
 *
 * @tparam ORDER interpolation order
 */
template <int ORDER> struct TensorTraits
{
	enum {nnodes = ORDER*ORDER*ORDER};
};



/**
 * @class FChebTensor
 *
 * The class FChebTensor provides function considering the tensor product
 * interpolation.
 *
 * @tparam ORDER interpolation order \f$\ell\f$
 */
template <int ORDER>
class FChebTensor : FNoCopyable
{
	enum {nnodes = TensorTraits<ORDER>::nnodes};
  typedef FChebRoots<ORDER> BasisType;
	
public:
	/**
	 * Sets the ids of the coordinates of all \f$\ell^3\f$ interpolation
	 * nodes
	 *
	 * @param[out] NodeIds ids of coordinates of interpolation nodes
	 */
	static
	void setNodeIds(unsigned int NodeIds[nnodes][3])
	{
		for (unsigned int n=0; n<nnodes; ++n) {
			NodeIds[n][0] =  n         % ORDER;
			NodeIds[n][1] = (n/ ORDER) % ORDER;
			NodeIds[n][2] =  n/(ORDER  * ORDER);
		}
	}


	/**
	 * Sets the roots of the Chebyshev quadrature weights defined as \f$w_i =
	 * \frac{\pi}{\ell}\sqrt{1-\bar x_i^2}\f$ with the Chebyshev roots \f$\bar
	 * x\f$.
	 *
	 * @param weights[out] the root of the weights \f$\sqrt{w_i}\f$
	 */
	static
	void setRootOfWeights(FReal weights[nnodes])
	{
		// weights in 1d
		FReal weights_1d[ORDER];
		for (unsigned int o=0; o<ORDER; ++o)
			weights_1d[o] = FMath::FPi/ORDER * FMath::Sqrt(FReal(1.)-FReal(BasisType::roots[o])*FReal(BasisType::roots[o]));
		// weights in 3d (tensor structure)
		unsigned int node_ids[nnodes][3];
		setNodeIds(node_ids);
		for (unsigned int n=0; n<nnodes; ++n) {
			weights[n] = FMath::Sqrt(weights_1d[node_ids[n][0]]*weights_1d[node_ids[n][1]]*weights_1d[node_ids[n][2]]);
		}
	}


	/**
	 * Sets the interpolation points in the cluster with @p center and @p width
	 *
	 * @param[in] center of cluster
	 * @param[in] width of cluster
	 * @param[out] rootPositions coordinates of interpolation points
	 */
	static
	void setRoots(const FPoint& center, const FReal width,
								FPoint rootPositions[nnodes])
	{
		unsigned int node_ids[nnodes][3];
		setNodeIds(node_ids);
		const map_loc_glob map(center, width);
		FPoint localPosition;
    for (unsigned int n=0; n<nnodes; ++n) {
			localPosition.setX(FReal(BasisType::roots[node_ids[n][0]]));
			localPosition.setY(FReal(BasisType::roots[node_ids[n][1]]));
			localPosition.setZ(FReal(BasisType::roots[node_ids[n][2]]));
			map(localPosition, rootPositions[n]);
		}
	}

	/**
	 * Set the relative child (width = 1) center according to the Morton index.
	 *
	 * @param[in] ChildIndex index of child according to Morton index
	 * @param[out] center
	 */
	static
	void setRelativeChildCenter(const unsigned int ChildIndex,
															FPoint& ChildCenter)
	{
		const int RelativeChildPositions[][3] = { {-1, -1, -1},
																							{-1, -1,  1},
																							{-1,  1, -1},
																							{-1,  1,  1},
																							{ 1, -1, -1},
																							{ 1, -1,  1},
																							{ 1,  1, -1},
																							{ 1,  1,  1} };
		ChildCenter.setX(FReal(RelativeChildPositions[ChildIndex][0]) / FReal(2.));
		ChildCenter.setY(FReal(RelativeChildPositions[ChildIndex][1]) / FReal(2.));
		ChildCenter.setZ(FReal(RelativeChildPositions[ChildIndex][2]) / FReal(2.));
	}
};





#endif
