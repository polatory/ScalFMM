#ifndef FFMBCOMPONENTS_HPP
#define FFMBCOMPONENTS_HPP

#include "../Extensions/FExtendForces.hpp"
#include "../Extensions/FExtendPotential.hpp"
#include "../Extensions/FExtendParticleType.hpp"
#include "../Extensions/FExtendCellType.hpp"

#include "../Components/FFmaParticle.hpp"
#include "../Components/FBasicCell.hpp"
#include "../Components/FSimpleLeaf.hpp"
#include "../Components/FBasicKernels.hpp"

#include "../Utils/FAbstractSendable.hpp"

#include "FExtendFmbCell.hpp"


class FmbParticle : public FExtendForces, public FFmaParticle, public FExtendPotential {
public:
};


class FmbCell : public FBasicCell, public FExtendFmbCell {
public:
};



class FmbTypedParticle : public FmbParticle, public FExtendParticleType {
public:
};


class FmbTypedCell : public FmbCell, public FExtendCellType {
public:
};



class FmbSendableCell : public FmbCell , public FAbstractSendable {
public:
    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    static const int SerializedSizeUp = sizeof(FComplexe)*MultipoleSize + sizeof(FBasicCell);
    void serializeUp(void* const buffer) const {
        memcpy(buffer, (FBasicCell*)this, sizeof(FBasicCell));
        memcpy((char*)(buffer) + sizeof(FBasicCell), multipole_exp, sizeof(FComplexe)*MultipoleSize );
    }
    void deserializeUp(const void* const buffer){
        memcpy((FBasicCell*)this, buffer, sizeof(FBasicCell));
        memcpy(multipole_exp, (char*)(buffer) + sizeof(FBasicCell), sizeof(FComplexe)*MultipoleSize );
    }

    static const int SerializedSizeDown = sizeof(FComplexe)*MultipoleSize + sizeof(FBasicCell);
    void serializeDown(void* const buffer) const {
        memcpy(buffer, (FBasicCell*)this, sizeof(FBasicCell));
        memcpy((char*)(buffer) + sizeof(FBasicCell), local_exp, sizeof(FComplexe)*MultipoleSize );
    }
    void deserializeDown(const void* const buffer){
        memcpy((FBasicCell*)this, buffer, sizeof(FBasicCell));
        memcpy(local_exp, (char*)(buffer) + sizeof(FBasicCell), sizeof(FComplexe)*MultipoleSize );
    }
};



#endif // FFMBCOMPONENTS_HPP
