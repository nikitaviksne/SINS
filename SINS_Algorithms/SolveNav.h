#ifndef SOLVENAV_H
#define SOLVENAV_H
#include "std.h"

void SolveNav(Ldoub* Wp, Ldoub* Ab, Ldoub* Cbn, Ldoub* Wo, Ldoub* Ao, Ldoub* V, Ldoub* Coordinates, Ldoub* CoordError, Ldoub* Err_V, Ldoub* omo, Ldoub h, Ldoub& Rphi, Ldoub& Rlambda, const Ldoub U, Ldoub R, Ldoub e, Ldoub H0, Ldoub* V0, Ldoub k1, Ldoub* ErrVins, bool allowCorr = false);

#endif