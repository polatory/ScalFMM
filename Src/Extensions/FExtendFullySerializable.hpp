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
#ifndef FEXTENDFULLYSENDABLE_HPP
#define FEXTENDFULLYSENDABLE_HPP

#include "../Utils/FAbstractSerializable.hpp"

#include "../Containers/FBufferWriter.hpp"
#include "../Containers/FBufferReader.hpp"

/** @author Berenger Bramas
  * This extension MUST be put as the first inherited class!
  * It copy the memory in the buffer read writer
  */
template <class ClassType>
class FExtendFullySerializable : public FAbstractSerializable {
public:
    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << (*reinterpret_cast<const ClassType*>(this));
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> (*reinterpret_cast<ClassType*>(this));
    }
};


#endif // FEXTENDFULLYSENDABLE_HPP
