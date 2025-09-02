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

#if 1
void diag(alglib::real_1d_array A, alglib::real_2d_array &A1)
{
	for(int i=0; i<A.length(); i++)
	{
		for(int j=0; j<A.length(); j++)
			A1[i][j]=0;
		A1[i][i] = A[i];
	}
		
}
#endif

void pinv(alglib::real_2d_array A, alglib::real_2d_array& A_pinv)
{
	// SVD ���������� �������
	#if 0
	printf("Function pinv\n");
	printf("Matrix A in pinv is:\n");
	print2dmatr(A);
	#endif
	
	int rows = A.rows();
	int cols = A.cols();
	alglib::real_1d_array W;
	alglib::real_2d_array U;
	alglib::real_2d_array Vt;

	W.setlength(rows);
	U.setlength(rows, cols);
	Vt.setlength(rows, cols);
	alglib::rmatrixsvd(A, rows, cols, 2, 2, 2, W, U, Vt);
#if 0
	printf("W =\n");
	for (int i=0; i<rows; i++)
	{
		printf("%e\t",W[i]);
	}
	printf("\n");
#endif
	alglib::real_2d_array S; // 2D ������� ���������� ��������
	S.setlength(rows, rows);
	diag(W, S);
#if 0
	printf("S =\n");
	print2dmatr(S);
	//  ����� ������� U
	printf("U =\n");
	print2dmatr(U);

	//  ����� ������� Vt
	printf("V^T =\n");
	print2dmatr(Vt);
#endif

	alglib::real_2d_array US;
	US.setlength(rows, cols);
	alglib::real_2d_array A1;
	A1.setlength(A.rows(), A.cols());
	rmatrixgemm(rows, cols, cols, 1, U, 0, 0, 0, S, 0, 0, 0, 0, US, 0, 0);
	rmatrixgemm(rows, cols, cols, 1, US, 0, 0, 0, Vt, 0, 0, 0, 0, A1, 0, 0);
	
#if 0
	printf("U @ S @ V^T =\n");
	print2dmatr(A1);
#endif
	#if 0
	alglib::matinvreport rep;
	rmatrixinverse(A1,cols, rep);
	printf("Report of the inverse:\t%d\n",rep);
	printf("Inverse matrix A:\n");
	print2dmatr(A1);
	#endif
	// S^{-1}
	for (int i=0; i<S.rows(); i++)
		if (S[i][i]!=0) S[i][i] = 1/S[i][i];
	// U^T
	alglib::real_2d_array Ut;
	Ut.setlength(rows, cols);
	// V
	alglib::real_2d_array V;
	V.setlength(rows, cols);
	//transpose U
	rmatrixtranspose(rows, cols, U, 0, 0, Ut, 0, 0);
	#if 0
	printf("Ut\t=\n");
	print2dmatr(Ut);
	#endif
	//transpose Vt
	rmatrixtranspose(rows, cols, Vt, 0, 0, V, 0, 0);
	#if 0
	printf("V\t=\n");
	print2dmatr(V);
	#endif
	// V*S^{-1}
	alglib::real_2d_array VS;
	VS.setlength(rows, cols);
	rmatrixgemm(rows, cols, cols, 1, V, 0, 0, 0, S, 0, 0, 0, 0, VS, 0, 0);
	// pinv of A
	
	rmatrixgemm(rows, cols, cols, 1, VS, 0, 0, 0, U, 0, 0, 1, 0, A_pinv, 0, 0);
	#if 0
	printf("Pseudo inverse of A\t=\n");
	print2dmatr(A_pinv);
	printf("Exit pinv\n");
	#endif
}