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
#ifndef FTESTPARTICLE_HPP
#define FTESTPARTICLE_HPP


#include "FBasicParticle.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTestParticle
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know whit it is.
*
* Particles just need the data down (after and run with FTestKernel the
* value shoud be NB PARTICLES (-1)).
*/
class FTestParticle : public FBasicParticle {
protected:
    // To store data during downard pass
    long long int dataDown;
public:
    FTestParticle(): dataDown(0){
    }

    /** Default destructor */
    virtual ~FTestParticle(){
    }

    /** Get the down data */
    long long int getDataDown() const {
        return this->dataDown;
    }

    /** Set down data */
    void setDataDown(const long long int inData){
        this->dataDown = inData;
    }

    //////////////////////////////////////////////////

    /** Save the current cell in a buffer */
    void save(FBufferWriter& buffer) const{
        FBasicParticle::save(buffer);
        buffer << dataDown;
    }

    /** Restore the current cell from a buffer */
    void restore(FBufferReader& buffer){
        FBasicParticle::restore(buffer);
        buffer >> dataDown;
    }
};


#endif //FTESTPARTICLE_HPP


