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
#ifndef FGLOBALPERIODIC_HPP
#define FGLOBALPERIODIC_HPP

///////////////////////////////////////////////////////
// Periodic condition definition
///////////////////////////////////////////////////////

enum PeriodicCondition {
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

bool testPeriodicCondition(const int conditions, const PeriodicCondition testConditions) {
    return (conditions & testConditions) == testConditions;
}

#endif // FGLOBALPERIODIC_HPP
