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



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractSendable
* Please read the license
*/
class FAbstractSendable {
protected:
    /** Empty Destructor */
    virtual ~FAbstractSendable(){}

    //static const int SerializedSizeUp = sizeof(?);
    virtual void serializeUp(void* const buffer) const  = 0;
    virtual void deserializeUp(const void* const buffer) = 0;

    //static const int SerializedSizeDown = sizeof(?);
    virtual void serializeDown(void* const buffer) const = 0;
    virtual void deserializeDown(const void* const buffer) = 0;
};

#endif //FABSTRACTSENDABLE_HPP


