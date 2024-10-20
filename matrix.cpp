#include <std.h>
#include <cmath>
void MulMatrD(double *A, double *B, double *C, int size1, int size2, int size3)
{
	/*
	A[size1*size2] @ B[size2*size3] = C[size1*size3]
	*/
	for (int i=0; i<size1; ++i) // по строкам A
		for(int j=0; j<size3; ++j) // по столбцам B
		{
			double sum=0;
			for (int k=0; k<size2; ++k) // по столбцам A (по строкам B)
			{
#if 0
				double a = A[index(size2, i, k)];
				double b = B[index(size3, k, j)];
#endif
				sum += (double) A[index(size2, i, k)] * (B[index(size3, k, j)]);
			}
			C[index(size3, i, j)] = (double) sum;
		}
}

void MatrOB(double H, double R, double P, double* C, int size)
{
	C[0] = (double) cos(R) * cos(H) + sin(R)*sin(H)* sin(P);
	C[1] = (double) sin(R) * cos(H)*sin(P) - cos(R)*sin(H);
	C[2] = (double) -sin(R) * cos(P);
	C[3] = (double) cos(P) * sin(H);
	C[4] = (double) cos(H) * cos(P);
	C[5] = (double) sin(P);
	C[6] = (double) sin(R) * cos(H) - cos(R) * sin(H) * sin(P);
	C[7] = (double) -sin(R) * sin(H) - cos(R) * cos(H) * sin(P);
	C[8] = (double) cos(R) * cos(P);
}