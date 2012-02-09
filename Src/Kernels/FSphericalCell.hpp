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
#ifndef FSPHERICALCELL_HPP
#define FSPHERICALCELL_HPP


#include "../Utils/FAbstractSendable.hpp"
#include "../Utils/FComplexe.hpp"
#include "../Utils/FMemUtils.hpp"
#include "../Extensions/FExtendCellType.hpp"
#include "../Components/FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
*/
class FSphericalCell : public FBasicCell {
protected:
    static int DevP;
    static int LocalSize;
    static int PoleSize;
    static bool UseBlas;

    FComplexe* multipole_exp; //< For multipole extenssion
    FComplexe* local_exp;     //< For local extenssion

public:
    static void Init(const int inDevP, const bool inUseBlas = false){
        DevP  = inDevP;
        const int ExpP  = int((inDevP+1) * (inDevP+2) * 0.5);
        const int NExpP = (inDevP+1) * (inDevP+1);

        LocalSize = ExpP;
        if(inUseBlas) {
            PoleSize = NExpP;
        }
        else{
            PoleSize = ExpP;
        }
    }

    static int GetLocalSize(){
        return LocalSize;
    }

    static int GetPoleSize(){
        return PoleSize;
    }

    /** Default constructor */
    FSphericalCell()
        : multipole_exp(0), local_exp(0){
        multipole_exp = new FComplexe[PoleSize];
        local_exp = new FComplexe[LocalSize];
    }

    /** Constructor */
    FSphericalCell(const FSphericalCell& other)
        : multipole_exp(0), local_exp(0){
        multipole_exp = new FComplexe[PoleSize];
        local_exp = new FComplexe[LocalSize];
        (*this) = other;
    }

    /** Default destructor */
    virtual ~FSphericalCell(){
        delete[] multipole_exp;
        delete[] local_exp;
    }

    /** Copy constructor */
    FSphericalCell& operator=(const FSphericalCell& other) {
        FMemUtils::copyall(multipole_exp, other.multipole_exp, PoleSize);
        FMemUtils::copyall(local_exp, other.local_exp, LocalSize);
        return *this;
    }

    /** Get Multipole */
    const FComplexe* getMultipole() const {
        return multipole_exp;
    }
    /** Get Local */
    const FComplexe* getLocal() const {
        return local_exp;
    }

    /** Get Multipole */
    FComplexe* getMultipole() {
        return multipole_exp;
    }
    /** Get Local */
    FComplexe* getLocal() {
        return local_exp;
    }
};

int FSphericalCell::DevP(-1);
int FSphericalCell::LocalSize(-1);
int FSphericalCell::PoleSize(-1);


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
*/
class FTypedSphericalCell : public FSphericalCell, public FExtendCellType {
public:
};


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
*/
class FSendableSphericalCell : public FSphericalCell , public FAbstractSendable {
public:
    static int SerializedSizeDown;
    static int SerializedSizeUp;

    static void Init(const int inDevP){
        FSphericalCell::Init(inDevP);
        SerializedSizeDown = PoleSize * sizeof(FComplexe) + sizeof(FBasicCell);
        SerializedSizeUp = LocalSize * sizeof(FComplexe) + sizeof(FBasicCell);
    }

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    void serializeUp(void* const buffer) const {
        memcpy(buffer, (FBasicCell*)this, sizeof(FBasicCell));
        memcpy((char*)(buffer) + sizeof(FBasicCell), multipole_exp, sizeof(FComplexe)*PoleSize );
    }
    void deserializeUp(const void* const buffer){
        memcpy((FBasicCell*)this, buffer, sizeof(FBasicCell));
        memcpy(multipole_exp, (char*)(buffer) + sizeof(FBasicCell), sizeof(FComplexe)*PoleSize );
    }

    void serializeDown(void* const buffer) const {
        memcpy(buffer, (FBasicCell*)this, sizeof(FBasicCell));
        memcpy((char*)(buffer) + sizeof(FBasicCell), local_exp, sizeof(FComplexe)*LocalSize );
    }
    void deserializeDown(const void* const buffer){
        memcpy((FBasicCell*)this, buffer, sizeof(FBasicCell));
        memcpy(local_exp, (char*)(buffer) + sizeof(FBasicCell), sizeof(FComplexe)*LocalSize );
    }
};

int FSendableSphericalCell::SerializedSizeDown(-1);
int FSendableSphericalCell::SerializedSizeUp(-1);

#endif //FSPHERICALCELL_HPP


