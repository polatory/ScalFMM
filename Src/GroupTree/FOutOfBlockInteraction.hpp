#ifndef FOUTOFBLOCKINTERACTION_HPP
#define FOUTOFBLOCKINTERACTION_HPP

#include "../Utils/FGlobal.hpp"

#include "FStarPUDefaultAlign.hpp"

struct  alignas(FStarPUDefaultAlign::StructAlign) OutOfBlockInteraction{
    MortonIndex outIndex;
    MortonIndex insideIndex;
    int outPosition;
    // To sort
    bool operator <=(const OutOfBlockInteraction& other) const{
        return outIndex <= other.outIndex;
    }
};

#endif // FOUTOFBLOCKINTERACTION_HPP

