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
#ifndef FBASICPARTICLE_HPP
#define FBASICPARTICLE_HPP


#include "../Extensions/FExtendPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticle
* Please read the license
*
* This class defines a basic particle used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractParticle without
* inheriting from it.
*/
class FBasicParticle : public FExtendPosition{
public:
    /** Save the current cell in a buffer */
    void save(FBufferWriter& buffer) const{
        FExtendPosition::save(buffer);
    }
    /** Restore the current cell from a buffer */
    void restore(FBufferReader& buffer){
        FExtendPosition::restore(buffer);
    }
};


#endif //FBASICPARTICLE_HPP


