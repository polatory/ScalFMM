#ifndef FCHEBINTERPOLATOR_HPP
#define FCHEBINTERPOLATOR_HPP


#include "./FChebMapping.hpp"
#include "./FChebTensor.hpp"
#include "./FChebRoots.hpp"

#include "../Utils/FBlas.hpp"



/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * Please read the license
 */

/**
 * @class FChebInterpolator
 *
 * The class @p FChebInterpolator defines the anterpolation (M2M) and
 * interpolation (L2L) concerning operations.
 */
template <int ORDER>
class FChebInterpolator : FNoCopyable
{
  // compile time constants and types
  enum {nnodes = TensorTraits<ORDER>::nnodes};
  typedef FChebRoots< ORDER>  BasisType;
  typedef FChebTensor<ORDER> TensorType;

  FReal T_of_roots[ORDER][ORDER];
	unsigned int node_ids[nnodes][3];
	FReal* ChildParentInterpolator[8];

	/**
	 * Initialize the child - parent - interpolator, it is basically the matrix
	 * S which is precomputed and reused for all M2M and L2L operations, ie for
	 * all non leaf inter/anterpolations.
	 */
	void initM2MandL2L()
	{
		F3DPosition ParentRoots[nnodes], ChildRoots[nnodes];
		const FReal ParentWidth(2.);
		const F3DPosition ParentCenter(0., 0., 0.);
		FChebTensor<ORDER>::setRoots(ParentCenter, ParentWidth, ParentRoots);

		F3DPosition ChildCenter;
		const FReal ChildWidth(1.);
		
		// loop: child cells
		for (unsigned int child=0; child<8; ++child) {

			// allocate memory
			ChildParentInterpolator[child] = new FReal [nnodes * nnodes];

			// set child info
			FChebTensor<ORDER>::setRelativeChildCenter(child, ChildCenter);
			FChebTensor<ORDER>::setRoots(ChildCenter, ChildWidth, ChildRoots);

			// assemble child - parent - interpolator
			assembleInterpolator(nnodes, ChildRoots, ChildParentInterpolator[child]);
		}
	}



public:
	/**
	 * Constructor: Initialize the Chebyshev polynomials at the Chebyshev
	 * roots/interpolation point
	 */
	explicit FChebInterpolator()
	{
		// initialize chebyshev polynomials of root nodes: T_o(x_j)
    for (unsigned int o=1; o<ORDER; ++o)
      for (unsigned int j=0; j<ORDER; ++j)
        T_of_roots[o][j] = FReal(BasisType::T(o, BasisType::roots[j]));

		// initialize root node ids
		TensorType::setNodeIds(node_ids);

		// initialize interpolation operator for non M2M and L2L (non leaf
		// operations)
		this -> initM2MandL2L();
	}

	
	/**
	 * Destructor: Delete dynamically allocated memory for M2M and L2L operator
	 */
	~FChebInterpolator()
	{
		for (unsigned int child=0; child<8; ++child)
			delete [] ChildParentInterpolator[child];
	}


	/**
	 * Assembles the interpolator \f$S_\ell\f$ of size \f$N\times
	 * \ell^3\f$. Here local points is meant as points whose global coordinates
	 * have already been mapped to the reference interval [-1,1].
	 *
	 * @param[in] NumberOfLocalPoints
	 * @param[in] LocalPoints
	 * @param[out] Interpolator
	 */
	void assembleInterpolator(const unsigned int NumberOfLocalPoints,
														const F3DPosition *const LocalPoints,
														FReal *const Interpolator) const
	{
		// values of chebyshev polynomials of source particle: T_o(x_i)
		FReal T_of_x[ORDER][3];

		// loop: local points (mapped in [-1,1])
		for (unsigned int m=0; m<NumberOfLocalPoints; ++m) {
			// evaluate chebyshev polynomials at local points
			for (unsigned int o=1; o<ORDER; ++o) {
				T_of_x[o][0] = BasisType::T(o, LocalPoints[m].getX());
				T_of_x[o][1] = BasisType::T(o, LocalPoints[m].getY());
				T_of_x[o][2] = BasisType::T(o, LocalPoints[m].getZ());
			}
			
			// assemble interpolator
			for (unsigned int n=0; n<nnodes; ++n) {
				Interpolator[n*nnodes + m] = FReal(1.);
				for (unsigned int d=0; d<3; ++d) {
					const unsigned int j = node_ids[n][d];
					FReal S_d = FReal(1.) / ORDER;
					for (unsigned int o=1; o<ORDER; ++o)
						S_d += FReal(2.) / ORDER * T_of_x[o][d] * T_of_roots[o][j];
					Interpolator[n*nnodes + m] *= S_d;
				}
			}
			
		}
		
	}

	
	/**
	 * Particle to moment: application of \f$S_\ell(y,\bar y_n)\f$
	 * (anterpolation, it is the transposed interpolation)
	 */
	template <class ContainerClass>
	void applyP2M(const F3DPosition& center,
								const FReal width,
								FReal *const multipoleExpansion,
								const ContainerClass *const sourceParticles) const
	{
		// setup global to local mapping
		const map_glob_loc map(center, width);
		F3DPosition localPosition;
		FReal T_of_x[ORDER][3];

		// set all multipole expansions to zero
		for (unsigned int n=0; n<nnodes; ++n) multipoleExpansion[n] = FReal(0.);

		// loop: source particles
		typename ContainerClass::ConstBasicIterator iter(*sourceParticles);
		while(iter.hasNotFinished()){

			// map global position to [-1,1]
			map(iter.data().getPosition(), localPosition);
			
			// get source value
			const FReal sourceValue = iter.data().getPhysicalValue();

			// evaluate chebyshev polynomials of source particle: T_o(x_i)
      for (unsigned int o=1; o<ORDER; ++o) {
        T_of_x[o][0] = BasisType::T(o, localPosition.getX());
        T_of_x[o][1] = BasisType::T(o, localPosition.getY());
        T_of_x[o][2] = BasisType::T(o, localPosition.getZ());
			}

			// anterpolate
			for (unsigned int n=0; n<nnodes; ++n) {
				//multipoleExpansion[n] = FReal(0.);
				FReal S = FReal(1.);
				for (unsigned int d=0; d<3; ++d) {
					const unsigned int j = node_ids[n][d];
					FReal S_d = FReal(1.) / ORDER;
					for (unsigned int o=1; o<ORDER; ++o)
						S_d += FReal(2.) / ORDER * T_of_x[o][d] * T_of_roots[o][j];
					S *= S_d;
				}
				multipoleExpansion[n] += S * sourceValue;
			}
			
			iter.gotoNext();
		} // end loop: source particles
	}


	
	/**
	 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ (interpolation)
	 */
	template <class ContainerClass>
	void applyL2P(const F3DPosition& center,
								const FReal width,
								const FReal *const localExpansion,
								ContainerClass *const localParticles) const
	{
		// setup local to global mapping
		const map_glob_loc map(center, width);
		F3DPosition localPosition;
		FReal T_of_x[ORDER][3];
		
		typename ContainerClass::BasicIterator iter(*localParticles);
		while(iter.hasNotFinished()){
			
			// map global position to [-1,1]
			map(iter.data().getPosition(), localPosition);
			
			// get target value
			FReal targetValue = iter.data().getPotential();

			// evaluate chebyshev polynomials of source particle: T_o(x_i)
      for (unsigned int o=1; o<ORDER; ++o) {
        T_of_x[o][0] = BasisType::T(o, localPosition.getX());
        T_of_x[o][1] = BasisType::T(o, localPosition.getY());
        T_of_x[o][2] = BasisType::T(o, localPosition.getZ());
			}

			// interpolate and increment target value
			for (unsigned int n=0; n<nnodes; ++n) {
				FReal S = FReal(1.);
				for (unsigned int d=0; d<3; ++d) {
					const unsigned int j = node_ids[n][d];
					FReal S_d = FReal(1.) / ORDER;
					for (unsigned int o=1; o<ORDER; ++o)
						S_d += FReal(2.) / ORDER * T_of_x[o][d] * T_of_roots[o][j];
					S *= S_d;
				}
				targetValue += S * localExpansion[n];
			}

			iter.data().setPotential(targetValue);

			iter.gotoNext();
		}
	}


	/**
	 * Local to particle operation: application of \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
	 */
	template <class ContainerClass>
	void applyL2PGradient(const F3DPosition& center,
												const FReal width,
												const FReal *const localExpansion,
												ContainerClass *const localParticles) const
	{
		// setup local to global mapping
		const map_glob_loc map(center, width);
		F3DPosition Jacobian;
		map.computeJacobian(Jacobian);
		const FReal jacobian[3] = {Jacobian.getX(), Jacobian.getY(), Jacobian.getZ()}; 
		F3DPosition localPosition;
		FReal T_of_x[ORDER][3];
		FReal U_of_x[ORDER][3];
		
		typename ContainerClass::BasicIterator iter(*localParticles);
		while(iter.hasNotFinished()){
			
			// map global position to [-1,1]
			map(iter.data().getPosition(), localPosition);
			
			// get target value
			FReal forces[3] = {iter.data().getForces().getX(),
												 iter.data().getForces().getY(),
												 iter.data().getForces().getZ()};

			// evaluate chebyshev polynomials of source particle
      for (unsigned int o=1; o<ORDER; ++o) {
				// T_o(x_i)
        T_of_x[o][0] = BasisType::T(o, localPosition.getX());
        T_of_x[o][1] = BasisType::T(o, localPosition.getY());
        T_of_x[o][2] = BasisType::T(o, localPosition.getZ());
				// T_o(x_i)
        U_of_x[o][0] = BasisType::U(o, localPosition.getX());
        U_of_x[o][1] = BasisType::U(o, localPosition.getY());
        U_of_x[o][2] = BasisType::U(o, localPosition.getZ());
			}

			// apply P
			for (unsigned int n=0; n<nnodes; ++n) {
				for (unsigned int i=0; i<3; ++i) {
					FReal P = FReal(1.);
					for (unsigned int d=0; d<3; ++d) {
						const unsigned int j = node_ids[n][d];
						FReal P_d;
						if (d==i) {
							P_d = 0.;
							for (unsigned int o=1; o<ORDER; ++o)
								P_d += 2. / ORDER * U_of_x[o][d] * T_of_roots[o][j]	* jacobian[d];
						} else {
							P_d = 1. / ORDER;
							for (unsigned int o=1; o<ORDER; ++o)
								P_d += 2. / ORDER * T_of_x[o][d] * T_of_roots[o][j];
						}
						P *= P_d;
					}
					// the minus sign comes due to the \f$- \nabla_x K(x,y) = \nabla_y K(x,y) = \bar K(x,y)\f$ 
					forces[i] -= P * localExpansion[n];
				}
			}
			iter.data().setForces(forces[0] * iter.data().getPhysicalValue(),
														forces[1] * iter.data().getPhysicalValue(),
														forces[2] * iter.data().getPhysicalValue());

			// increment iterator
			iter.gotoNext();
		}
	}

	
	void applyM2M(const unsigned int ChildIndex,
								const FReal *const ChildExpansion,
								FReal *const ParentExpansion) const
	{
		FBlas::gemtva(nnodes, nnodes, FReal(1.),
									ChildParentInterpolator[ChildIndex],
									const_cast<FReal *const>(ChildExpansion), ParentExpansion);
	}

	void applyL2L(const unsigned int ChildIndex,
								const FReal *const ParentExpansion,
								FReal *const ChildExpansion) const
	{
		FBlas::gemva(nnodes, nnodes, FReal(1.),
								 ChildParentInterpolator[ChildIndex],
								 const_cast<FReal *const>(ParentExpansion), ChildExpansion);
	}
	
};



#endif
