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
#ifndef FABSTRACTSERIALIZABLE_HPP
#define FABSTRACTSERIALIZABLE_HPP

class FBufferReader;
class FBufferWriter;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractSerializable
* Please read the license
*
* This propose an interface to save and restore a class
*
* To make your particles are usable in the mpi fmm,
* they must provide this interface
*/
class FAbstractSerializable {
protected:
    virtual void save(FBufferWriter&) const  = 0;
    virtual void restore(FBufferReader&) = 0;
};

#endif // FABSTRACTSERIALIZABLE_HPP
