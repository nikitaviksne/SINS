#include "std.h"
#include "mathematics_kalman.h"
#include "ap.h"
#include "linalg.h"
//Умножение одномерных массивов как матриц
#if 1
template <class T> void matMul(int aRows, int aCols, int bCols, T A, T  B, T& C)
{
	for(int iii=0; iii<aRows; ++iii) //по строкам массива А
	{
		for (int kkk=0; kkk<bCols; kkk++)// по столбцам массива B
		{
			double sum (0);
			for (int jjj=0; jjj<aCols; ++jjj) //По столбцам массива B
			{
				 sum = A[index_3(aCols, iii, jjj)] * B[index_3(bCols, jjj, kkk)];
			}
			C[index_3(bCols, iii, kkk)] = sum ;
		}
	}

}
#endif

#if 1
//Перегруженная функция умножение одномерных массивов alglib как матриц
void matMul(int aRows, int aCols, int bCols, alglib::real_1d_array A, alglib::real_1d_array B, alglib::real_1d_array& C)
{
	for(int iii=0; iii<aRows; ++iii) //по строкам массива А
	{
		for (int kkk=0; kkk<bCols; kkk++)// по столбцам массива B
		{
			double sum (0);
			for (int jjj=0; jjj<aCols; ++jjj) //По столбцам массива A и строкам массива B
			{
				 sum += A[index_3(aCols, iii, jjj)] * B[index_3(bCols, jjj, kkk)];
			}
			C[index_3(bCols, iii, kkk)] = sum ;
		}
	}

}
#endif

void transpose(int size, alglib::real_1d_array B, alglib::real_1d_array& B_t) //для квадратных матриц
{
	for (int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			B_t[index_3(size, i,j)] = B[index_3(size, j,i)];
}

void transpose(int Brows, int Bcols, alglib::real_1d_array B, alglib::real_1d_array& B_t) //для прямоугольных матриц
{
	for (int i=0; i<Brows; i++)
		for(int j=0; j<Bcols; j++)
			B_t[index_3(Brows, j,i)] = B[index_3(Bcols, i,j)];
}

#if 0
void matMul(int aRows)
{
	aRows++;
}
#endif
