// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
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
