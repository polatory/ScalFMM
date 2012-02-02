#ifndef FCHEBCELL_HPP
#define FCHEBCELL_HPP
// [--License--]

#include "../Extensions/FExtendPosition.hpp"
#include "../Extensions/FExtendMortonIndex.hpp"
#include "../Extensions/FExtendCoordinate.hpp"

#include "./FChebTensor.hpp"

/**
* @author Matthias Messner (matthias.messner@inria.fr)
* @class FChebCell
* Please read the license
*
* This class defines a cell used in the Chebyshev based FMM.
*/
template <int ORDER>
class FChebCell : public FExtendPosition,
									public FExtendMortonIndex,
									public FExtendCoordinate
{
	FReal multipole_exp[TensorTraits<ORDER>::nnodes]; //< Multipole expansion
	FReal     local_exp[TensorTraits<ORDER>::nnodes]; //< Local expansion
	
public:
	~FChebCell() {}
	
	/** Get Multipole */
	const FReal *const getMultipole() const
	{	return this->multipole_exp;	}
	/** Get Local */
	const FReal *const getLocal() const
	{	return this->local_exp;	}
	
	/** Get Multipole */
	FReal *const getMultipole()
	{	return this->multipole_exp;	}
	/** Get Local */
	FReal *const getLocal()
	{	return this->local_exp;	}
	
};


#endif //FCHEBCELL_HPP

// [--END--]