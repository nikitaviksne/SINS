#ifndef MATRIXH
#define MATRIXH

void MulMatrD(Ldoub *A, Ldoub *B, Ldoub *C, int size1, int size2, int size3);
void Transpose(Ldoub *, int);
void Transpose2M(Ldoub *, Ldoub *, int);
void print2dMatr(Ldoub*, int, int);

bool invers(int size,float *a);
bool inversD(int size,Ldoub *a);
#endif
