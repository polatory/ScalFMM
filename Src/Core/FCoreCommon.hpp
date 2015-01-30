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
#ifndef FCORECOMMON_HPP
#define FCORECOMMON_HPP

#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssert.hpp"

/**
 * @brief The FFmmOperations enum
 * To chose which operation has to be performed.
 */
enum FFmmOperations {
    FFmmP2P = (1 << 0),
    FFmmP2M = (1 << 1),
    FFmmM2M = (1 << 2),
    FFmmM2L = (1 << 3),
    FFmmL2L = (1 << 4),
    FFmmL2P = (1 << 5),

    FFmmNearField = FFmmP2P,
    FFmmFarField  = (FFmmP2M|FFmmM2M|FFmmM2L|FFmmL2L|FFmmL2P),

    FFmmNearAndFarFields = (FFmmNearField|FFmmFarField)
};

/**
 * @brief The FAbstractAlgorithm class
 * Is an abstract algorithm to be able to use the FAlgorithmBuilder
 * and execute from an abastrct pointer
 */
class FAbstractAlgorithm {
protected:
    //< Where to start the work
    int upperWorkingLevel;
    //< Where to end the work (exclusif)
    int lowerWorkingLevel;
    //< Height of the tree
    int nbLevelsInTree;

    void setNbLevelsInTree(const int inNbLevelsInTree){
        nbLevelsInTree = inNbLevelsInTree;
        lowerWorkingLevel = nbLevelsInTree;
    }

    void validateLevels() const {
        FAssertLF(FAbstractAlgorithm::upperWorkingLevel <= FAbstractAlgorithm::lowerWorkingLevel);
        FAssertLF(2 <= FAbstractAlgorithm::upperWorkingLevel);
    }

    virtual void executeCore(const unsigned operationsToProceed) = 0;

public:
    FAbstractAlgorithm()
        : upperWorkingLevel(2), lowerWorkingLevel(0), nbLevelsInTree(-1){
    }

    virtual ~FAbstractAlgorithm(){
    }

    /** Execute all the fmm but for given levels*/
    virtual void execute(const int inUpperWorkingLevel, const int inLowerWorkingLevel) final {
        upperWorkingLevel = inUpperWorkingLevel;
        lowerWorkingLevel = inLowerWorkingLevel;
        validateLevels();
        executeCore(FFmmNearAndFarFields);
    }

    /** Execute all the fmm */
    virtual void execute() final {
        upperWorkingLevel = 2;
        lowerWorkingLevel = nbLevelsInTree;
        validateLevels();
        executeCore(FFmmNearAndFarFields);
    }

    /** Execute only some FMM operation for given levels */
    virtual void execute(const unsigned operationsToProceed, const int inUpperWorkingLevel, const int inLowerWorkingLevel) final {
        upperWorkingLevel = inUpperWorkingLevel;
        lowerWorkingLevel = inLowerWorkingLevel;
        validateLevels();
        executeCore(operationsToProceed);
    }

    /** Execute only some steps */
    virtual void execute(const unsigned operationsToProceed) final {
        upperWorkingLevel = 2;
        lowerWorkingLevel = nbLevelsInTree;
        validateLevels();
        executeCore(operationsToProceed);
    }
};




#endif // FCORECOMMON_HPP
