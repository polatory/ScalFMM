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
#ifndef FROTATIONCELL_HPP
#define FROTATIONCELL_HPP

#include "../../Components/FAbstractSerializable.hpp"
#include "../../Components/FAbstractSendable.hpp"
#include "../../Utils/FComplexe.hpp"
#include "../../Utils/FMemUtils.hpp"

#include "../../Extensions/FExtendCellType.hpp"

#include "../../Components/FBasicCell.hpp"

#include "../../Containers/FBufferWriter.hpp"
#include "../../Containers/FBufferReader.hpp"

/** This class is a cell used for the rotation based kernel
  * The size of the multipole and local vector are based on a template
  * User should choose this parameter P carrefuly to match with the
  * P of the kernel.
  *
  * Multipole/Local vectors contain value as:
  * {0,0}{1,0}{1,1}...{P,P-1}{P,P}
  * So the size of such vector can be obtained by a suite:
  * (n+1)*n/2 => (P+2)*(P+1)/2
  */
template <int P>
class FRotationCell : public FBasicCell {
protected:
    //< Size of multipole vector
    static const int MultipoleSize = ((P+2)*(P+1))/2; // Artimethique suite (n+1)*n/2
    //< Size of local vector
    static const int LocalSize = ((P+2)*(P+1))/2;     // Artimethique suite (n+1)*n/2

    //< Multipole vector (static memory)
    FComplexe multipole_exp[MultipoleSize]; //< For multipole extenssion
    //< Local vector (static memory)
    FComplexe local_exp[LocalSize];         //< For local extenssion

public:
    /** Default constructor
      * Put 0 in vectors
      */
    FRotationCell(){
    }

    /** Copy constructor
      * Copy the value in the vectors
      */
    FRotationCell(const FRotationCell& other){
        (*this) = other;
    }

    /** Default destructor */
    virtual ~FRotationCell(){
    }

    /** Copy operator
      * copies only the value in the vectors
      */
    FRotationCell& operator=(const FRotationCell& other) {
        FMemUtils::copyall(multipole_exp, other.multipole_exp, MultipoleSize);
        FMemUtils::copyall(local_exp, other.local_exp, LocalSize);
        return *this;
    }

    /** Get Multipole array */
    const FComplexe* getMultipole() const {
        return multipole_exp;
    }
    /** Get Local array */
    const FComplexe* getLocal() const {
        return local_exp;
    }

    /** Get Multipole array */
    FComplexe* getMultipole() {
        return multipole_exp;
    }
    /** Get Local array */
    FComplexe* getLocal() {
        return local_exp;
    }

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    void serializeUp(FBufferWriter& buffer) const{
        buffer.write(multipole_exp, MultipoleSize);
    }
    void deserializeUp(FBufferReader& buffer){
        buffer.fillArray(multipole_exp, MultipoleSize);
    }

    void serializeDown(FBufferWriter& buffer) const{
        buffer.write(local_exp, LocalSize);
    }
    void deserializeDown(FBufferReader& buffer){
        buffer.fillArray(local_exp, LocalSize);
    }

    ///////////////////////////////////////////////////////
    // to extend Serializable
    ///////////////////////////////////////////////////////
    void save(FBufferWriter& buffer) const{
        FBasicCell::save(buffer);
        buffer.write(multipole_exp, MultipoleSize);
        buffer.write(local_exp, LocalSize);
    }
    void restore(FBufferReader& buffer){
        FBasicCell::restore(buffer);
        buffer.fillArray(multipole_exp, MultipoleSize);
        buffer.fillArray(local_exp, LocalSize);
    }
};

#endif // FROTATIONCELL_HPP
