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
#ifndef FFmaPARTICLE_HPP
#define FFmaPARTICLE_HPP


#include "../Extensions/FExtendPosition.hpp"
#include "../Extensions/FExtendPhysicalValue.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaParticle
* Please read the license
*
* This class defines a particle for FMA loader.
* As defined in FFmaLoader it needs {FBasicParticle,FExtendPhysicalValue}
*/
class FFmaParticle : public FExtendPosition, public FExtendPhysicalValue {
public:
    /** Save the current cell in a buffer */
    void save(FBufferWriter& buffer) const{
        FExtendPosition::save(buffer);
        FExtendPhysicalValue::save(buffer);
    }
    /** Restore the current cell from a buffer */
    void restore(FBufferReader& buffer){
        FExtendPosition::restore(buffer);
        FExtendPhysicalValue::restore(buffer);
    }
};


#endif //FFmaPARTICLE_HPP


