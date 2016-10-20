// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#ifndef FGLOBALPERIODIC_HPP
#define FGLOBALPERIODIC_HPP

///////////////////////////////////////////////////////
// Periodic condition definition
///////////////////////////////////////////////////////

/**
 * @brief The PeriodicCondition enum
 * To be able to chose the direction of the periodicity.
 */
enum PeriodicCondition {
    DirNone     = 0,

    DirPlusX    = 1 << 0,
    DirMinusX   = 1 << 1,
    DirPlusY    = 1 << 2,
    DirMinusY   = 1 << 3,
    DirPlusZ    = 1 << 4,
    DirMinusZ   = 1 << 5,

    DirX        = (DirPlusX | DirMinusX),
    DirY        = (DirPlusY | DirMinusY),
    DirZ        = (DirPlusZ | DirMinusZ),

    AllDirs     = (DirX | DirY | DirZ)
};

/**
 * @brief TestPeriodicCondition
 * @param conditions
 * @param testConditions
 * @return true if the direction is in the condition
 */
inline bool TestPeriodicCondition(const int conditions, const PeriodicCondition testConditions) {
    return (conditions & testConditions) == testConditions;
}

#endif // FGLOBALPERIODIC_HPP
