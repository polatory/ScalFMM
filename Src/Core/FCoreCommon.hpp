#ifndef FCORECOMMON_HPP
#define FCORECOMMON_HPP

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
public:
    virtual ~FAbstractAlgorithm(){
    }

    /** Execute all the fmm */
    void execute(){
        execute(FFmmNearAndFarFields);
    }

    /** Execute only some steps */
    virtual void execute(const unsigned operationsToProceed) = 0;
};




#endif // FCORECOMMON_HPP
