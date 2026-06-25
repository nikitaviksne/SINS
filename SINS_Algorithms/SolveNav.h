#ifndef SOLVENAV_H
#define SOLVENAV_H
#include "std.h"
#include "ap.h"

extern const Ldoub U;
extern const Ldoub R;
extern const Ldoub g;
extern const Ldoub a;
extern const Ldoub b;
extern const Ldoub e;

void SolveNav(Ldoub* Wp, Ldoub* Ab, Ldoub* Cbn, Ldoub* Wo, Ldoub* Ao, Ldoub* V, Ldoub* Coordinates, Ldoub* CoordError, Ldoub* Err_V, Ldoub* omo, Ldoub h, Ldoub& Rphi, Ldoub& Rlambda, Ldoub H0, Ldoub* V0, Ldoub k1, alglib::real_1d_array ErrVins, bool allowCorr = false);

#endif