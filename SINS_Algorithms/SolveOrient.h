#ifndef SOLVEORIENT_H
#define SOLVEORIENT_H

#include "quaternions.h"
#include "std.h"

void SolveOrient(Ldoub* Omb, quaternion* Qf, Ldoub* Cbn, Ldoub* Orientation, Ldoub* Coordinates, Ldoub* Omo, Ldoub* omo, Ldoub* V, Ldoub phi0,Ldoub Rphi, Ldoub Rlambda, int freq, Ldoub h, const Ldoub U, int iter, bool& sw, Ldoub k2, Ldoub* DVerr, bool allowCorr = false);
#endif