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
