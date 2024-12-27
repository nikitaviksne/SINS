#include<std.h>
#include "mathematics.h"
#include <stdio.h>

void MulMatrD(Ldoub *A, Ldoub *B, Ldoub *C, int size1, int size2, int size3)
{
	/*
	A[size1*size2] @ B[size2*size3] = C[size1*size3]
	*/
	for (int i=0; i<size1; ++i) // по строкам A
		for(int j=0; j<size3; ++j) // по столбцам B
		{
			Ldoub sum=0;
			Ldoub errProd=0;
			Ldoub errSum=0;
			for (int k=0; k<size2; ++k) // по столбцам A (по строкам B)
			{
#if 1
				Ldoub tempRes=0;
				// умножение
				TwoProduct(A[index(size2, i, k)], B[index(size3, k, j)], tempRes, errProd);
				// Сложение
				TwoSum(sum, tempRes, sum, errSum, false);
				// Компенсация накопленной ошибки
				TwoSum(sum, errSum, sum, errSum, false);
#endif
				//sum += (Ldoub) A[index(size2, i, k)] * (B[index(size3, k, j)]);
			}
			C[index(size3, i, j)] = sum;
		}
}

void Transpose(Ldoub *A, int size1)
{
	Ldoub buf;
	for(int i=0; i<size1; ++i)
		for(int j=0; j<size1; ++j)
			if (i!=j)
			{
				buf = A[index(size1, i,j)];
				A[index(size1, i,j)] = A[index(size1, j,i)];
				A[index(size1, j,i)] = buf;
			}
}

void Transpose2M(Ldoub *A, Ldoub *B, int size1)
{
	for(int i=0; i<size1; ++i)
		for(int j=0; j<size1; ++j)
				B[index(size1, i,j)] = A[index(size1, j,i)];
}

void print2dMatr(Ldoub* A, int rows, int cols)
{
	for(int iii=0; iii<rows; ++iii)
	{
		for(int jjj=0; jjj<cols; ++jjj)
			printf("%.20f\t", A[index(cols, iii, jjj)]);
		printf("\n");
	}
}