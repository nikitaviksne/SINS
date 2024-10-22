#include <std.h>
#include <cmath>
#include "mathematics.h"
#include <stdio.h>

void MulMatrD(double *A, double *B, double *C, int size1, int size2, int size3)
{
	/*
	 * A[size1*size2] @ B[size2*size3] = C[size1*size3]
	 */
	for (int i=0; i<size1; ++i) // по строкам A
		for(int j=0; j<size3; ++j) // по столбцам B
		{
			double sum=0;
			double errProd=0;
			double errSum=0;
			for (int k=0; k<size2; ++k) // по столбцам A (по строкам B)
			{
				#if 1
				double tempRes=0;
				// умножение
				TwoProduct(A[index(size2, i, k)], B[index(size3, k, j)], tempRes, errProd);
				// Сложение
				TwoSum(sum, tempRes, sum, errSum, false);
				// Компенсация накопленной ошибки
				TwoSum(sum, errSum, sum, errSum, false);
				#endif
				//sum += (double) A[index(size2, i, k)] * (B[index(size3, k, j)]);
			}
			C[index(size3, i, j)] = sum;
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

void print2dMatr(double* A, int size)
{
	for(int iii=0; iii<size; ++iii)
	{
		for(int jjj=0; jjj<size; ++jjj)
			printf("%.20f\t", A[index(size, iii,jjj)]);
		printf("\n");
	}
}
