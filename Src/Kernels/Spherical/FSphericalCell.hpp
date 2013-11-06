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
#ifndef FSPHERICALCELL_HPP
#define FSPHERICALCELL_HPP

#include "../../Utils/FComplexe.hpp"
#include "../../Utils/FMemUtils.hpp"

#include "../../Extensions/FExtendCellType.hpp"

#include "../../Components/FBasicCell.hpp"

#include "../../Containers/FBufferWriter.hpp"
#include "../../Containers/FBufferReader.hpp"

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

    /** Make it like the begining */
    void resetToInitialState(){
        for(int idx = 0 ; idx < PoleSize ; ++idx){
            multipole_exp[idx].setRealImag(FReal(0.0), FReal(0.0));
        }
        for(int idx = 0 ; idx < LocalSize ; ++idx){
            local_exp[idx].setRealImag(FReal(0.0), FReal(0.0));
        }
    }

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    void serializeUp(FBufferWriter& buffer) const{
        buffer.write(multipole_exp, PoleSize);
    }
    void deserializeUp(FBufferReader& buffer){
        buffer.fillArray(multipole_exp, PoleSize);
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
        buffer.write(multipole_exp, PoleSize);
        buffer.write(local_exp, LocalSize);
    }
    void restore(FBufferReader& buffer){
        FBasicCell::restore(buffer);
        buffer.fillArray(multipole_exp, PoleSize);
        buffer.fillArray(local_exp, LocalSize);
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
    void save(FBufferWriter& buffer) const{
        FSphericalCell::save(buffer);
        FExtendCellType::save(buffer);
    }
    void restore(FBufferReader& buffer){
        FSphericalCell::restore(buffer);
        FExtendCellType::restore(buffer);
    }
    void resetToInitialState(){
        FSphericalCell::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }
};



#endif //FSPHERICALCELL_HPP


