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
#ifndef FABSTRACTSENDABLE_HPP
#define FABSTRACTSENDABLE_HPP

class FBufferReader;
class FBufferWriter;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractSendable
* Please read the license
*
* To make your cells are usable in the mpi fmm,
* they must provide this interface
*/
class FAbstractSendable {
protected:
    /** Empty Destructor */
    virtual ~FAbstractSendable(){}

    ///////////////////////////////////////////////
    // For Upward pass
    ///////////////////////////////////////////////

    /** Save your data */
    virtual void serializeUp(FBufferWriter&) const  = 0;
    /** Retrieve your data */
    virtual void deserializeUp(FBufferReader&) = 0;

    ///////////////////////////////////////////////
    // For Downward pass
    ///////////////////////////////////////////////

    /** Save your data */
    virtual void serializeDown(FBufferWriter&) const = 0;
    /** Retrieve your data */
    virtual void deserializeDown(FBufferReader&) = 0;
};


#endif //FABSTRACTSENDABLE_HPP


