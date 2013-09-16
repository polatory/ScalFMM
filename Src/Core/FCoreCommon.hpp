#ifndef FCORECOMMON_HPP
#define FCORECOMMON_HPP

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

#endif // FCORECOMMON_HPP
