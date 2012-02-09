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
#include "FTrace.hpp"
#ifdef SCALFMM_USE_TRACE

#ifndef SCALFMM_USE_ITAC
int FTrace::Deep = 0;
FTic FTrace::TimeSinceBegining;
#endif //SCALFMM_USE_ITAC

#endif // SCALFMM_USE_TRACE

