#ifndef FCHEBTENSOR_HPP
#define FCHEBTENSOR_HPP


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
	 * Sets the interpolation points in the cluster with @p center and @p width
	 *
	 * @param[in] center of cluster
	 * @param[in] width of cluster
	 * @param[out] rootPositions coordinates of interpolation points
	 */
	static
	void setRoots(const F3DPosition& center, const FReal width,
								F3DPosition rootPositions[nnodes])
	{
		unsigned int node_ids[nnodes][3];
		setNodeIds(node_ids);
		const map_loc_glob map(center, width);
		F3DPosition localPosition;
    for (unsigned int n=0; n<nnodes; ++n) {
			localPosition.setX(FReal(BasisType::roots[node_ids[n][0]]));
			localPosition.setY(FReal(BasisType::roots[node_ids[n][1]]));
			localPosition.setZ(FReal(BasisType::roots[node_ids[n][2]]));
			map(localPosition, rootPositions[n]);
		}
	}
};





#endif
