#ifndef MATRIXH
#define MATRIXH
#include "std.h"
#include "ap.h"

void MulMatrD(Ldoub *A, Ldoub *B, Ldoub *C, int size1, int size2, int size3);
void Transpose(Ldoub *, int);
void Transpose2M(Ldoub *, Ldoub *, int);
void print2dMatr(Ldoub*, int, int);
void print2dMatr(const alglib::real_1d_array& , int , int );

bool invers(int size,float *a);
bool inversD(int size,Ldoub *a);
void Determinant_3(Ldoub* A, Ldoub& res);
void Normalization(Ldoub* Cbn, bool sw);
void Ortogonalization(Ldoub *Cbn, bool sw);
#endif
