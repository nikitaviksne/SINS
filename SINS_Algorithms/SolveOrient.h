#ifndef SOLVEORIENT_H
#define SOLVEORIENT_H

#include "quaternions.h"
#include "std.h"
#include "ap.h"

extern const Ldoub U;
extern const Ldoub R;
extern const Ldoub g;
extern const Ldoub a;
extern const Ldoub b;
extern const Ldoub e;


void SolveOrient(Ldoub* Omb, quaternion* Qf, Ldoub* Cbn, Ldoub* Orientation, Ldoub* Coordinates, Ldoub* Omo, Ldoub* omo, Ldoub* V, Ldoub phi0,Ldoub Rphi, Ldoub Rlambda, int freq, Ldoub h, int iter, bool& sw, Ldoub k2, alglib::real_1d_array errPhi, bool allowCorr = false);
#endif