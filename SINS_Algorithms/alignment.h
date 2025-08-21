#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include <std.h>
#include "quaternions.h"

void alignment(/*входные*/Ldoub *Ab, Ldoub *Omb, Ldoub *MeanAb, Ldoub *MeanOmb, Ldoub *StdAb, Ldoub *StdOmb, int iter, const Ldoub g, const Ldoub U, Ldoub phi0, /*выходные*/Ldoub *Cbn, quaternion *Qf);
#endif
