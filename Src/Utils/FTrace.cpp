#include "FTrace.hpp"
#ifdef SCALFMM_USE_TRACE

#ifndef SCALFMM_USE_ITAC
int FTrace::Deep = 0;
FTic FTrace::TimeSinceBegining;
#endif //SCALFMM_USE_ITAC

#endif // SCALFMM_USE_TRACE

