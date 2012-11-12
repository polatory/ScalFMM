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
#ifndef FBASICKERNELS_HPP
#define FBASICKERNELS_HPP


#include "FAbstractKernels.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels is empty and does nothing.
*/
template< class ParticleClass, class CellClass, class ContainerClass>
class FBasicKernels : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {
public:
    /** Default destructor */
    virtual ~FBasicKernels(){
    }

    /** Do nothing */
    virtual void P2M(CellClass* const , const ContainerClass* const ) {

    }

    /** Do nothing */
    virtual void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int ) {

    }

    /** Do nothing */
    virtual void M2L(CellClass* const FRestrict , const CellClass* [], const int , const int ) {

    }

    /** Do nothing */
    virtual void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) {

    }

    /** Do nothing */
    virtual void L2P(const CellClass* const , ContainerClass* const ){

    }


    /** Do nothing */
    virtual void P2P(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [27], const int ){

    }

    /** Do nothing */
    virtual void P2PRemote(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [27], const int ){

    }

};


#endif //FBASICKERNELS_HPP


