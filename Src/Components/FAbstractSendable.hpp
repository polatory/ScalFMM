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


