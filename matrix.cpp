#include "std.h"
#include "mathematics.h"
#include <stdio.h>

void MulMatrD(Ldoub *A, Ldoub *B, Ldoub *C, int size1, int size2, int size3)
{
	/*
	A[size1*size2] @ B[size2*size3] = C[size1*size3]
	*/
	for (int i=0; i<size1; ++i) // �� ������� A
		for(int j=0; j<size3; ++j) // �� �������� B
		{
			Ldoub sum=0;
			Ldoub errProd=0;
			Ldoub errSum=0;
			for (int k=0; k<size2; ++k) // �� �������� A (�� ������� B)
			{
#if 1
				Ldoub tempRes=0;
				// ���������
				TwoProduct(A[index_3(size2, i, k)], B[index_3(size3, k, j)], tempRes, errProd);
				// ��������
				TwoSum(sum, tempRes, sum, errSum, false);
				// ����������� ����������� ������
				TwoSum(sum, errSum, sum, errSum, false);
#endif
				//sum += (Ldoub) A[index_3(size2, i, k)] * (B[index_3(size3, k, j)]);
			}
			C[index_3(size3, i, j)] = sum;
		}
}

void Transpose(Ldoub *A, int size1)
{
	Ldoub buf;
	for(int i=0; i<size1; ++i)
		for(int j=0; j<size1; ++j)
			if (i!=j)
			{
				buf = A[index_3(size1, i,j)];
				A[index_3(size1, i,j)] = A[index_3(size1, j,i)];
				A[index_3(size1, j,i)] = buf;
			}
}

void Transpose2M(Ldoub *A, Ldoub *B, int size1)
{
	for(int i=0; i<size1; ++i)
		for(int j=0; j<size1; ++j)
				B[index_3(size1, i,j)] = A[index_3(size1, j,i)];
}

void print2dMatr(Ldoub* A, int rows, int cols)
{
	for(int iii=0; iii<rows; ++iii)
	{
		for(int jjj=0; jjj<cols; ++jjj)
			printf("%.20f\t", A[index_3(cols, iii, jjj)]);
		printf("\n");
	}
}

bool invers(int size,float *a)
{
/*Взял из ArmPsi*/
   float p,y,*a1,*b1,*c1;
   int i,j,k;
   for (k=0;k<size;k++){
	  if(a[0]<=0){
		 return false;
	  }
	  p=1/a[0];
	  a1=a+size*(size-1);
	  b1=a+1;
	  for(i=1;i<size;i++)*a1++=*b1++;
	  for (i=0;i<size-1;i++){
		 y=-a[index_3(size,size-1,i)]*p;
		 a[index_3(size,i,size-1)]=y;
		 a1=a+size*i+i;
		 b1=a+size*(i+1)+(i+1);
		 c1=a+size*(size-1)+i;
		 for (j=i;j<size-1;j++)*a1++=(*b1++)+(*c1++)*y;
	  }
	  a[index_3(size,size-1,size-1)]=-p;
   }
   for (i=0;i<size;i++){
	  a1=a+size*i+i;
	  b1=a1;
	  for (j=i;j<size;j++){
		 *a1=-(*a1);
		 *b1=*a1++;
		 b1+=size;
	  }
   }
   return true;
}

bool inversD(int size,Ldoub *a)
{
/*Взял из ArmPsi*/
   Ldoub p,y,*a1,*b1,*c1;
   int i,j,k;
   for (k=0;k<size;k++){
      if(a[0]<=0){
         return false;
      }
      p = 1/a[0]; //p = D_DIV(Constant_1,a[0]);
      a1 = a+size*(size-1);
      b1 = a+1;
      for(i=1;i<size;i++)*a1++=*b1++;
      for (i=0;i<size-1;i++){
         y= - a[index_3(size,size-1,i)]*p;        // y=D_MINUS(D_MUL(a[index_3(size,size-1,i)],p));
         a[index_3(size,i,size-1)]=y;
         a1=a+size*i+i;
         b1=a+size*(i+1)+(i+1);
         c1=a+size*(size-1)+i;
         for (j=i;j<size-1;j++)
		*a1 ++= (*b1++) + (*c1++) * y;
//*a1++=D_ADD((*b1++),D_MUL((*c1++),y));
      }
      a[index_3(size,size-1,size-1)]=-(p);
   }
   for (i=0;i<size;i++){
      a1=a+size*i+i;
      b1=a1;
      for (j=i;j<size;j++){
         *a1= - *a1;
         *b1 = *a1++;
         b1 += size;
      }
   }
   return true;
}

void Determinant_3(Ldoub* A, Ldoub& res)
{ //Определитель матрицы 3*3
	//Методом "звездочки"
	res = A[0] * A[index_3(3, 1, 1)]* A[index_3(3,2,2)] + A[index_3(3, 0, 1)] * A[index_3(3, 1, 0)] * A[index_3(3, 2, 0)] + A[index_3(3, 1, 0)] * A[index_3(3, 2, 1)]* A[index_3(3, 0, 2)] - (A[index_3(3, 2, 0)] * A[index_3(3, 1, 1)] * A[index_3(3, 0, 2)] + A[index_3(3, 1, 0)] * A[index_3(3, 0, 1)] * A[index_3(3, 2, 2)] + A[index_3(3, 2, 1)] * A[index_3(3, 1, 2)] * A[index_3(3, 0, 0)]);
}

void Normalization(Ldoub* Cbn, bool sw)
{
	Ldoub norm[1];
	//контроль масштаба
	if (sw)//Строки
		for (int iii=0; iii<3; ++iii)
			{
				Ldoub string[3] = {Cbn[index_3(3,iii,0)], Cbn[index_3(3,iii,1)], Cbn[index_3(3,iii,2)]};
				MulMatrD(string, string, norm, 1, 3, 1);
				norm[0] = 1 - norm[0];
				for (int jjj=0; jjj < 3; ++jjj)
					Cbn[index_3(3,iii,jjj)] = Cbn[index_3(3,iii,jjj)]  - 0.5*norm[0]*Cbn[index_3(3,iii,jjj)];
			}
	else//Столбцы
		for (int iii=0; iii<3; ++iii)
			{
				Ldoub string[3] = {Cbn[index_3(3,0,iii)], Cbn[index_3(3,1,iii)], Cbn[index_3(3,2,iii)]};
				MulMatrD(string, string, norm, 1, 3, 1);
				norm[0] = 1 - norm[0];
				for (int jjj=0; jjj < 3; ++jjj)
					Cbn[index_3(3,jjj,iii)] = Cbn[index_3(3,jjj,iii)]  - 0.5*norm[0]*Cbn[index_3(3,jjj,iii)];
			}
}

void Ortogonalization(Ldoub *Cbn, bool sw)
{
	// ортогонализация
	if (sw)//по строкам
		for(int iii=0; iii < 3; ++iii)
			for(int jjj=0; jjj < 3; ++jjj)
			{
				if(iii==jjj) continue;
				Ldoub string[3] = {Cbn[index_3(3,iii,0)], Cbn[index_3(3,iii,1)], Cbn[index_3(3,iii,2)]};
				Ldoub column[3] = {Cbn[index_3(3,jjj,0)], Cbn[index_3(3,jjj,1)], Cbn[index_3(3,jjj,2)]};
				Ldoub norm[1];
				MulMatrD(string, column, norm, 1,3,1);
				for(int kkk=0; kkk<3; ++kkk)
				{
					Cbn[index_3(3,iii,kkk)] = Cbn[index_3(3,iii,kkk)] - 0.5*norm[0]*column[kkk];// Cbn[index_3(3,jjj,kkk )];
					Cbn[index_3(3,jjj,kkk)] = Cbn[index_3(3,jjj,kkk)] - 0.5*norm[0]*string[kkk];// Cbn[index_3(3,iii,kkk)];
				}
			}
	else//по столбцам
		for(int iii=0; iii < 3; ++iii)
			for(int jjj=0; jjj < 3; ++jjj)
			{
				if(iii==jjj) continue;
				Ldoub string[3] = {Cbn[index_3(3,0,iii)], Cbn[index_3(3,1,iii)], Cbn[index_3(3,2,iii)]};
				Ldoub column[3] = {Cbn[index_3(3,0,jjj)], Cbn[index_3(3,1,jjj)], Cbn[index_3(3,2,jjj)]};
				Ldoub norm[1];
				MulMatrD(string, column,norm, 1,3,1);
				for(int kkk=0; kkk<3; ++kkk)
				{
					Cbn[index_3(3,kkk,iii)] = Cbn[index_3(3,kkk,iii)] - 0.5*norm[0]*column[iii];//Cbn[index_3(3,kkk,jjj )];
					Cbn[index_3(3,kkk,jjj)] = Cbn[index_3(3,kkk,jjj)] - 0.5*norm[0]*string[iii];//Cbn[index_3(3,kkk,iii)];
				}
			}

}