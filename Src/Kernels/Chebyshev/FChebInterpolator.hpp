#ifndef FCHEBINTERPOLATOR_HPP
#define FCHEBINTERPOLATOR_HPP


#include "./FChebMapping.hpp"
#include "./FChebTensor.hpp"
#include "./FChebRoots.hpp"

#include "../../Utils/FBlas.hpp"



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
		FPoint ParentRoots[nnodes], ChildRoots[nnodes];
		const FReal ParentWidth(2.);
		const FPoint ParentCenter(0., 0., 0.);
		FChebTensor<ORDER>::setRoots(ParentCenter, ParentWidth, ParentRoots);

		FPoint ChildCenter;
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
        T_of_roots[o][j] = FReal(BasisType::T(o, FReal(BasisType::roots[j])));

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
														const FPoint *const LocalPoints,
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
	void applyP2M(const FPoint& center,
								const FReal width,
								FReal *const multipoleExpansion,
								const ContainerClass *const sourceParticles) const;


	
	/**
	 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ (interpolation)
	 */
	template <class ContainerClass>
	void applyL2P(const FPoint& center,
								const FReal width,
								const FReal *const localExpansion,
								ContainerClass *const localParticles) const;


	/**
	 * Local to particle operation: application of \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
	 */
	template <class ContainerClass>
	void applyL2PGradient(const FPoint& center,
												const FReal width,
												const FReal *const localExpansion,
												ContainerClass *const localParticles) const;

	/**
	 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ and
	 * \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
	 */
	template <class ContainerClass>
	void applyL2PTotal(const FPoint& center,
										 const FReal width,
										 const FReal *const localExpansion,
										 ContainerClass *const localParticles) const;
	
	
	void applyM2M(const unsigned int ChildIndex,
								const FReal *const ChildExpansion,
								FReal *const ParentExpansion) const
	{
		FBlas::gemtva(nnodes, nnodes, FReal(1.),
									ChildParentInterpolator[ChildIndex],
									const_cast<FReal*>(ChildExpansion), ParentExpansion);
	}

	void applyL2L(const unsigned int ChildIndex,
								const FReal *const ParentExpansion,
								FReal *const ChildExpansion) const
	{
		FBlas::gemva(nnodes, nnodes, FReal(1.),
								 ChildParentInterpolator[ChildIndex],
								 const_cast<FReal*>(ParentExpansion), ChildExpansion);
	}
	
};







/**
 * Particle to moment: application of \f$S_\ell(y,\bar y_n)\f$
 * (anterpolation, it is the transposed interpolation)
 */
template <int ORDER>
template <class ContainerClass>
inline void FChebInterpolator<ORDER>::applyP2M(const FPoint& center,
																							 const FReal width,
																							 FReal *const multipoleExpansion,
																							 const ContainerClass *const sourceParticles) const
{
	// set all multipole expansions to zero
	FBlas::setzero(nnodes, multipoleExpansion);

	// allocate stuff
	const map_glob_loc map(center, width);
	FPoint localPosition;
	FReal T_of_x[ORDER][3];
	FReal S[3];

	// loop over source particles
	typename ContainerClass::ConstBasicIterator iter(*sourceParticles);
	while(iter.hasNotFinished()){

		// map global position to [-1,1]
		map(iter.data().getPosition(), localPosition);

		// evaluate chebyshev polynomials of source particle: T_o(x_i)
		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
		for (unsigned int o=2; o<ORDER; ++o) {
			T_of_x[o][0] = FReal(2.) * localPosition.getX() * T_of_x[o-1][0] - T_of_x[o-2][0];
			T_of_x[o][1] = FReal(2.) * localPosition.getY() * T_of_x[o-1][1] - T_of_x[o-2][1];
			T_of_x[o][2] = FReal(2.) * localPosition.getZ() * T_of_x[o-1][2] - T_of_x[o-2][2];
		}
		
		// anterpolate
		const FReal sourceValue = iter.data().getPhysicalValue();
		for (unsigned int n=0; n<nnodes; ++n) {
			const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};
			S[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
			S[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
			S[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				S[0] += T_of_x[o][0] * T_of_roots[o][j[0]];
				S[1] += T_of_x[o][1] * T_of_roots[o][j[1]];
				S[2] += T_of_x[o][2] * T_of_roots[o][j[2]];
			}
			// gather contributions
			S[0] *= FReal(2.); S[0] += FReal(1.);
			S[1] *= FReal(2.); S[1] += FReal(1.);
			S[2] *= FReal(2.); S[2] += FReal(1.);
			multipoleExpansion[n]	+= FReal(1.) / nnodes *	S[0] * S[1] * S[2] *	sourceValue;
		}
		// increment source iterator
		iter.gotoNext();
	}
}


/**
 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ (interpolation)
 */
template <int ORDER>
template <class ContainerClass>
inline void FChebInterpolator<ORDER>::applyL2P(const FPoint& center,
																							 const FReal width,
																							 const FReal *const localExpansion,
																							 ContainerClass *const localParticles) const
{
	// allocate stuff
	const map_glob_loc map(center, width);
	FPoint localPosition;
	FReal T_of_x[ORDER][3];
	FReal S[3];
		
	typename ContainerClass::BasicIterator iter(*localParticles);
	while(iter.hasNotFinished()){
			
		// map global position to [-1,1]
		map(iter.data().getPosition(), localPosition);

		// evaluate chebyshev polynomials of source particle: T_o(x_i)
		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
		for (unsigned int o=2; o<ORDER; ++o) {
			T_of_x[o][0] = FReal(2.) * localPosition.getX() * T_of_x[o-1][0] - T_of_x[o-2][0];
			T_of_x[o][1] = FReal(2.) * localPosition.getY() * T_of_x[o-1][1] - T_of_x[o-2][1];
			T_of_x[o][2] = FReal(2.) * localPosition.getZ() * T_of_x[o-1][2] - T_of_x[o-2][2];
		}

		// interpolate and increment target value
		FReal targetValue = iter.data().getPotential();
		for (unsigned int n=0; n<nnodes; ++n) {
			const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};
			S[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
			S[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
			S[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				S[0] += T_of_x[o][0] * T_of_roots[o][j[0]];
				S[1] += T_of_x[o][1] * T_of_roots[o][j[1]];
				S[2] += T_of_x[o][2] * T_of_roots[o][j[2]];
			}
			// gather contributions
			S[0] *= FReal(2.); S[0] += FReal(1.);
			S[1] *= FReal(2.); S[1] += FReal(1.);
			S[2] *= FReal(2.); S[2] += FReal(1.);
			targetValue	+= S[0] * S[1] * S[2] * localExpansion[n];
		}
		// scale
		targetValue /= nnodes;

		// set potential
		iter.data().setPotential(targetValue);
		// increment target iterator
		iter.gotoNext();
	}
}






/**
 * Local to particle operation: application of \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
 */
template <int ORDER>
template <class ContainerClass>
inline void FChebInterpolator<ORDER>::applyL2PGradient(const FPoint& center,
																											 const FReal width,
																											 const FReal *const localExpansion,
																											 ContainerClass *const localParticles) const
{
	// setup local to global mapping
	const map_glob_loc map(center, width);
	FPoint Jacobian;
	map.computeJacobian(Jacobian);
	const FReal jacobian[3] = {Jacobian.getX(), Jacobian.getY(), Jacobian.getZ()}; 
	FPoint localPosition;
	FReal T_of_x[ORDER][3];
	FReal U_of_x[ORDER][3];
	FReal P[3];

	typename ContainerClass::BasicIterator iter(*localParticles);
	while(iter.hasNotFinished()){
			
		// map global position to [-1,1]
		map(iter.data().getPosition(), localPosition);
			
		// evaluate chebyshev polynomials of source particle
		// T_0(x_i) and T_1(x_i)
		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
		// U_0(x_i) and U_1(x_i)
		U_of_x[0][0] = FReal(1.);	U_of_x[1][0] = localPosition.getX() * FReal(2.);
		U_of_x[0][1] = FReal(1.);	U_of_x[1][1] = localPosition.getY() * FReal(2.);
		U_of_x[0][2] = FReal(1.);	U_of_x[1][2] = localPosition.getZ() * FReal(2.);
		for (unsigned int o=2; o<ORDER; ++o) {
			// T_o(x_i)
			T_of_x[o][0] = FReal(2.)*localPosition.getX()*T_of_x[o-1][0] - T_of_x[o-2][0];
			T_of_x[o][1] = FReal(2.)*localPosition.getY()*T_of_x[o-1][1] - T_of_x[o-2][1];
			T_of_x[o][2] = FReal(2.)*localPosition.getZ()*T_of_x[o-1][2] - T_of_x[o-2][2];
			// U_o(x_i)
			U_of_x[o][0] = FReal(2.)*localPosition.getX()*U_of_x[o-1][0] - U_of_x[o-2][0];
			U_of_x[o][1] = FReal(2.)*localPosition.getY()*U_of_x[o-1][1] - U_of_x[o-2][1];
			U_of_x[o][2] = FReal(2.)*localPosition.getZ()*U_of_x[o-1][2] - U_of_x[o-2][2];
		}

		// scale, because dT_o/dx = oU_{o-1}
		for (unsigned int o=2; o<ORDER; ++o) {
			U_of_x[o-1][0] *= FReal(o);
			U_of_x[o-1][1] *= FReal(o);
			U_of_x[o-1][2] *= FReal(o);
		}

		// apply P and increment forces
		FReal forces[3] = {FReal(0.), FReal(0.), FReal(0.)};
		for (unsigned int n=0; n<nnodes; ++n) {
			
			// tensor indices of chebyshev nodes
			const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};

			// f0 component //////////////////////////////////////
			P[0] = U_of_x[0][0] * T_of_roots[1][j[0]];
			P[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
			P[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				P[0] += U_of_x[o-1][0] * T_of_roots[o][j[0]];
				P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]];
				P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]];
			}
			P[0] *= FReal(2.);
			P[1] *= FReal(2.); P[1] += FReal(1.);
			P[2] *= FReal(2.); P[2] += FReal(1.);
			forces[0]	+= P[0] * P[1] * P[2] * localExpansion[n];

			// f1 component //////////////////////////////////////
			P[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
			P[1] = U_of_x[0][1] * T_of_roots[1][j[1]];
			P[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]];
				P[1] += U_of_x[o-1][1] * T_of_roots[o][j[1]];
				P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]];
			}
			P[0] *= FReal(2.); P[0] += FReal(1.);
			P[1] *= FReal(2.); 
			P[2] *= FReal(2.); P[2] += FReal(1.);
			forces[1]	+= P[0] * P[1] * P[2] * localExpansion[n];

			// f2 component //////////////////////////////////////
			P[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
			P[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
			P[2] = U_of_x[0][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]];
				P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]];
				P[2] += U_of_x[o-1][2] * T_of_roots[o][j[2]];
			}
			P[0] *= FReal(2.); P[0] += FReal(1.);
			P[1] *= FReal(2.); P[1] += FReal(1.);
			P[2] *= FReal(2.);
			forces[2]	+= P[0] * P[1] * P[2] * localExpansion[n];
		}

		// scale forces
		forces[0] *= jacobian[0] / nnodes;
		forces[1] *= jacobian[1] / nnodes;
		forces[2] *= jacobian[2] / nnodes;

		// set computed forces
		iter.data().incForces(forces[0] * iter.data().getPhysicalValue(),
													forces[1] * iter.data().getPhysicalValue(),
													forces[2] * iter.data().getPhysicalValue());

		// increment iterator
		iter.gotoNext();
	}
}



/**
 * Local to particle operation: application of \f$S_\ell(x,\bar x_m)\f$ and
 * \f$\nabla_x S_\ell(x,\bar x_m)\f$ (interpolation)
 */
template <int ORDER>
template <class ContainerClass>
inline void FChebInterpolator<ORDER>::applyL2PTotal(const FPoint& center,
																										const FReal width,
																										const FReal *const localExpansion,
																										ContainerClass *const localParticles) const
{
	// setup local to global mapping
	const map_glob_loc map(center, width);
	FPoint Jacobian;
	map.computeJacobian(Jacobian);
	const FReal jacobian[3] = {Jacobian.getX(), Jacobian.getY(), Jacobian.getZ()}; 
	FPoint localPosition;
	FReal T_of_x[ORDER][3];
	FReal U_of_x[ORDER][3];
	FReal P[3];

	typename ContainerClass::BasicIterator iter(*localParticles);
	while(iter.hasNotFinished()){
			
		// map global position to [-1,1]
		map(iter.data().getPosition(), localPosition);
			
		// evaluate chebyshev polynomials of source particle
		// T_0(x_i) and T_1(x_i)
		T_of_x[0][0] = FReal(1.);	T_of_x[1][0] = localPosition.getX();
		T_of_x[0][1] = FReal(1.);	T_of_x[1][1] = localPosition.getY();
		T_of_x[0][2] = FReal(1.);	T_of_x[1][2] = localPosition.getZ();
		// U_0(x_i) and U_1(x_i)
		U_of_x[0][0] = FReal(1.);	U_of_x[1][0] = localPosition.getX() * FReal(2.);
		U_of_x[0][1] = FReal(1.);	U_of_x[1][1] = localPosition.getY() * FReal(2.);
		U_of_x[0][2] = FReal(1.);	U_of_x[1][2] = localPosition.getZ() * FReal(2.);
		for (unsigned int o=2; o<ORDER; ++o) {
			// T_o(x_i)
			T_of_x[o][0] = FReal(2.)*localPosition.getX()*T_of_x[o-1][0] - T_of_x[o-2][0];
			T_of_x[o][1] = FReal(2.)*localPosition.getY()*T_of_x[o-1][1] - T_of_x[o-2][1];
			T_of_x[o][2] = FReal(2.)*localPosition.getZ()*T_of_x[o-1][2] - T_of_x[o-2][2];
			// U_o(x_i)
			U_of_x[o][0] = FReal(2.)*localPosition.getX()*U_of_x[o-1][0] - U_of_x[o-2][0];
			U_of_x[o][1] = FReal(2.)*localPosition.getY()*U_of_x[o-1][1] - U_of_x[o-2][1];
			U_of_x[o][2] = FReal(2.)*localPosition.getZ()*U_of_x[o-1][2] - U_of_x[o-2][2];
		}

		// scale, because dT_o/dx = oU_{o-1}
		for (unsigned int o=2; o<ORDER; ++o) {
			U_of_x[o-1][0] *= FReal(o);
			U_of_x[o-1][1] *= FReal(o);
			U_of_x[o-1][2] *= FReal(o);
		}

		// apply P and increment forces
		FReal potential = FReal(0.);
		FReal forces[3] = {FReal(0.), FReal(0.), FReal(0.)};
		for (unsigned int n=0; n<nnodes; ++n) {
			
			// tensor indices of chebyshev nodes
			const unsigned int j[3] = {node_ids[n][0], node_ids[n][1], node_ids[n][2]};

			// f0 component //////////////////////////////////////
			P[0] = U_of_x[0][0] * T_of_roots[1][j[0]];
			P[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
			P[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
			// for potential
			FReal S0 = T_of_x[1][0] * T_of_roots[1][j[0]];
			for (unsigned int o=2; o<ORDER; ++o) {
				P[0] += U_of_x[o-1][0] * T_of_roots[o][j[0]];
				P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]];
				P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]];
				// for potential
				S0   += T_of_x[o][0] * T_of_roots[o][j[0]];
			}
			P[0] *= FReal(2.);
			P[1] *= FReal(2.); P[1] += FReal(1.);
			P[2] *= FReal(2.); P[2] += FReal(1.);
			forces[0]	+= P[0] * P[1] * P[2] * localExpansion[n];
			// for potential
			S0 *= FReal(2.); S0 += FReal(1.);
			potential	+= FReal(1.) / nnodes *	S0 * P[1] * P[2] * localExpansion[n];

			// f1 component //////////////////////////////////////
			P[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
			P[1] = U_of_x[0][1] * T_of_roots[1][j[1]];
			P[2] = T_of_x[1][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]];
				P[1] += U_of_x[o-1][1] * T_of_roots[o][j[1]];
				P[2] += T_of_x[o  ][2] * T_of_roots[o][j[2]];
			}
			P[0] *= FReal(2.); P[0] += FReal(1.);
			P[1] *= FReal(2.); 
			P[2] *= FReal(2.); P[2] += FReal(1.);
			forces[1]	+= P[0] * P[1] * P[2] * localExpansion[n];

			// f2 component //////////////////////////////////////
			P[0] = T_of_x[1][0] * T_of_roots[1][j[0]];
			P[1] = T_of_x[1][1] * T_of_roots[1][j[1]];
			P[2] = U_of_x[0][2] * T_of_roots[1][j[2]];
			for (unsigned int o=2; o<ORDER; ++o) {
				P[0] += T_of_x[o  ][0] * T_of_roots[o][j[0]];
				P[1] += T_of_x[o  ][1] * T_of_roots[o][j[1]];
				P[2] += U_of_x[o-1][2] * T_of_roots[o][j[2]];
			}
			P[0] *= FReal(2.); P[0] += FReal(1.);
			P[1] *= FReal(2.); P[1] += FReal(1.);
			P[2] *= FReal(2.);
			forces[2]	+= P[0] * P[1] * P[2] * localExpansion[n];
		}

		// scale forces
		forces[0] *= jacobian[0] / nnodes;
		forces[1] *= jacobian[1] / nnodes;
		forces[2] *= jacobian[2] / nnodes;

		// set computed potential
		iter.data().incPotential(potential);

		// set computed forces
		iter.data().incForces(forces[0] * iter.data().getPhysicalValue(),
													forces[1] * iter.data().getPhysicalValue(),
													forces[2] * iter.data().getPhysicalValue());

		// increment iterator
		iter.gotoNext();
	}
}


#endif
