#include "ap.h"
#include "linalg.h"
#include "AdaptiveRKalman.h"
#include "mathematics_kalman.h"
#include "std.h"
#include <stdio.h>

// #define dev

AdaptiveRKalman::AdaptiveRKalman(int dimx, int dimz, int dimq, Ldoub h)
{
	setDimX(dimx);
	setDimZ(dimz);
	setDimQ(dimq); // установка размера матрицы q
	this->h = h;

	Papr.setlength(dim_x * dim_x);
	A.setlength(dim_x * dim_x);
	A2.setlength(dim_x * dim_x);
	A3.setlength(dim_x * dim_x);
	Phi.setlength(dim_x * dim_x);
	#if 0
	Phi[index_3(dim_x, 0, 0)] = 1;	Phi[index_3(dim_x, 0, 1)] = 1;
	Phi[index_3(dim_x, 1, 0)] = 0;	Phi[index_3(dim_x, 1, 1)] = 1;
	#endif

	Papst.setlength(dim_x * dim_x);

	x.setlength(dim_x * 1);
	x_1.setlength(dim_x * 1);
	G.setlength(dim_x*dim_q); // dim_x * dim_q
	GQGt.setlength(dim_x*dim_x); // dim_x * dim_x
	Q.setlength(dim_q * dim_q);
	C.setlength(dim_z * dim_z);
	R.setlength(dim_z * dim_z);
	K.setlength(dim_x * dim_z);
	H.setlength(dim_z * dim_x);
	v.setlength(dim_z * 1);

	I.setlength(dim_x * dim_x);
	for(int ii=0; ii<dim_x; ++ii)
	{
		for(int jj=0; jj<dim_x; ++jj)
			I[index_3(dim_x, ii,jj)] = 0;
		I[index_3(dim_x, ii, ii)] = 1;
	}

	z.setlength(dim_z * 1);

	//Инициализация ковариационных матриц
	//Априорная ошибка оценивания
	for(int iii=0; iii<getDimX(); ++iii)
		for(int jjj=0; jjj<getDimX(); ++jjj)
		{
			Papr[index_3(dim_x, iii, jjj)] = 0; //априорная ошибка оценивания

			Papst[index_3(dim_x, iii, jjj)] = 0; //апостериорная ошибка оценивания

		}

	

			
#if 0
	Papr[index_3(dim_x, 0, 0)] = 0;	Papr[index_3(dim_x, 0, 1)] = 0;
	Papr[index_3(dim_x, 1, 0)] = 0;	Papr[index_3(dim_x, 1, 1)] = 0;

	//Апостериорная ошибка оценивания
	Papst[index_3(dim_x, 0, 0)] = 0;	Papst[index_3(dim_x, 0, 1)] = 0;
	Papst[index_3(dim_x, 1, 0)] = 0;	Papst[index_3(dim_x, 1, 1)] = 0;

	//Входного шума
	Q[index_3(dim_x, 0, 0)] = q1;	Q[index_3(dim_x, 0, 1)] = 0;
	Q[index_3(dim_x, 1, 0)] = 0;	Q[index_3(dim_x, 1, 1)] = q2;

	//Измерительного шума
	R[index_3(dim_x, 0, 0)] = r1;	R[index_3(dim_x, 0, 1)] = 0;
	R[index_3(dim_x, 1, 0)] = 0;	R[index_3(dim_x, 1, 1)] = q2;
#endif
	init = false;
}
void AdaptiveRKalman::Init(Ldoub* initVal, Ldoub* q, Ldoub* r, Ldoub* h/*матрица измерений*/)
{
	Init(initVal, q, h);
}
void AdaptiveRKalman::Init(Ldoub* initVal, Ldoub* q, Ldoub* h/*матрица измерений*/)
{//Инициализация значениями


	//Инициализация начальными значениями
	for (int iii=0; iii<getDimX(); ++iii)
	{
		x[iii] = initVal[iii];
	}
	for(int iii=0; iii < getDimQ(); ++iii)
		for(int jjj=0; jjj < getDimQ(); ++jjj)
			Q[index_3(getDimQ(), iii, jjj)] = q[index_3(getDimQ(), iii, jjj)];

	for (int iii=0; iii<getDimZ(); ++iii)
		for(int jjj=0; jjj<getDimX(); ++jjj)
			H[index_3(dim_x, iii, jjj)] = h[index_3(dim_x, iii, jjj)];

	iter = 0;
	//Флаг инициализации
	init=true;
}
void AdaptiveRKalman::Predict()
{
	matMul(getDimX(), getDimX(), getDimX(), A, A, A2);
	matMul(getDimX(), getDimX(), getDimX(), A2, A, A3);
	for (int iii =0; iii < getDimX()*getDimX(); ++iii) //вычисляю матрицу перехода Phi
		Phi[iii] = I[iii] + A[iii] * h + A2[iii] * pow(h,2) / 2. + A3[iii] * pow(h, 3) / 6.; //не забываем умножить на такт интегрирования
	matMul(getDimX(), getDimX(), 1, Phi, x, x_1); //предсказываем вектор состояния x_1  = Phi @ x

	/*предсказываем априорную ошибку оценивания (используя апостериорную и ковариацию входного шума)*/

	alglib::real_1d_array Phi_Papst;
	Phi_Papst.setlength(getDimX() * getDimX());//временная переменная Phi @ Papst
	matMul(getDimX(), getDimX(), getDimX(), Phi, Papst, Phi_Papst);

	alglib::real_1d_array Phi_Papst_Phi_t;
	Phi_Papst_Phi_t.setlength(getDimX() * getDimX());//временная переменная Phi @ Papst @ Phi.T

	alglib::real_1d_array Phi_t;
	Phi_t.setlength(getDimX() * getDimX()); //временная переменная Phi.T
	transpose(getDimX(), Phi, Phi_t);

	matMul(getDimX(), getDimX(), getDimX(), Phi_Papst, Phi_t, Phi_Papst_Phi_t);

#if 0
	printf("Матрица G\n");
	Print2dMatr(G, getDimX(), getDimQ());
#endif

	//Вспомогательная матрица GQ
	alglib::real_1d_array GQ;
	GQ.setlength(getDimX()*getDimQ());
	matMul(getDimX(), getDimQ(), getDimQ(), G, Q, GQ);

#if 0
	printf("Матрица GQ\n");
	Print2dMatr(GQ, getDimX(), getDimQ());
#endif
	// делаем вспомогательную матрицу Gt = G^T
	alglib::real_1d_array Gt; 
	Gt.setlength(getDimQ()*getDimX()); 
	transpose(getDimX(), getDimQ(), G, Gt);

#if 0
	printf("Матрица G^t\n");
	Print2dMatr(Gt, getDimQ(), getDimX());
#endif
	//Вычисляем матрицу GQGt
	matMul(getDimX(), getDimQ(), getDimX(), GQ, Gt, GQGt);

#if 0
	printf("Матрица QQG^t\n");
	Print2dMatr(GQGt, getDimX(), getDimX());
#endif

	for (int iii=0; iii<getDimX()*getDimX(); ++iii)
		Papr[iii] = Phi_Papst_Phi_t[iii] + GQGt[iii];

}

void AdaptiveRKalman::Update(Ldoub* zin/*измерения обычные C-массивы*/)
{
	for(int i=0; i< getDimZ(); i++)
		this->z[i] = zin[i];
	//вычисляем различия между измерениями и предсказанием
	alglib::real_1d_array Hx_1;
	Hx_1.setlength(getDimZ() * 1);
	matMul(getDimZ(), getDimX(), 1, H, x_1, Hx_1);
	for (int iii=0; iii<getDimZ(); ++iii)
		v[iii] = this->z[iii] - Hx_1[iii];
#ifdef dev
	printf("v:\n");
	Print2dMatr(v, getDimZ(), 1);
#endif
	//Вычисляем v*v.T
	alglib::real_1d_array vvt;
	vvt.setlength(getDimZ() * getDimZ());
	for(int i=0; i<getDimZ(); i++)
	{
		for(int j=0; j<getDimZ(); j++) //все элементы по 0
		{
			vvt[index_3(getDimZ(), i, j)] = 0;
		}
		vvt[index_3(getDimZ(), i, i)] = pow(v[i],2);//главная диагональ с дисперсиями
	}
	// matMul(getDimZ(), 1, getDimZ(), v, v, vvt);

	//Вычислям C
	for(int iii=0; iii<getDimZ(); iii++) 
			for(int jjj=0; jjj<getDimZ(); ++jjj)
			{
				Ldoub buf;
				buf = iter / (iter + 1.) * (C[index_3(getDimZ(), iii, jjj)]) + 1. / (iter + 1.) *  (vvt[index_3(getDimZ(), iii, jjj)]); // для стационарной системы
				C[index_3(getDimZ(), iii, jjj)] = buf;
			}
#ifdef dev
	printf("C:\n");
	Print2dMatr(C, getDimZ(), getDimZ());
#endif
	
	alglib::real_1d_array H_t;
	H_t.setlength(getDimX() * getDimZ());
	transpose(getDimZ(), getDimX(), H, H_t);

	alglib::real_1d_array Papr_H_t; //потом еще понадобится
	Papr_H_t.setlength(getDimX() * getDimZ());//временная переменная H @ Papr
	matMul(getDimX(), getDimX(), getDimZ(), Papr, H_t, Papr_H_t);

#ifdef dev
// #if 1
	printf("Papr_H_t:\n");
	Print2dMatr(Papr_H_t, getDimX(), getDimZ());
#endif

	alglib::real_1d_array H_Papr_H_t;
	H_Papr_H_t.setlength(getDimZ() * getDimZ());//временная переменная H @ Papr @ H.T
	matMul(getDimZ(), getDimX(), getDimZ(), H, Papr_H_t, H_Papr_H_t);

#ifdef dev
	printf("H_Papr_H_t:\n");
	Print2dMatr(H_Papr_H_t, getDimZ(), getDimZ());
#endif

	//Вычисляем  R
	for(int iii=0; iii<getDimZ()*getDimZ(); ++iii)
	{
		Ldoub temp = C[iii] - H_Papr_H_t[iii];		
		R[iii] = temp;
	}
	
#if 1
	for (int iii=0; iii<getDimZ()/**getDimZ()*/; ++iii)
		if (R[index_3(getDimZ(), iii, iii)] < 0) //проверка на отрицательность ---диагональных элементов C---
		{
			for(int kkk=0; kkk<getDimZ()/**getDimZ()*/; ++kkk) //все диагональные элементы устанавливаем в 0
				R[index_3(getDimZ(), kkk, kkk)] = 0;
			break; //и выходим из внешнего цикла
		}
#else //значения из обычного фильтра Калмана
Ldoub freq = 100.0;
	R[0] = pow(0.2*sqrt(freq), 2); R[1] = 0; R[2] = 0; R[3] = 0;
	R[4] = 0; R[5] = pow(0.2*sqrt(freq), 2); R[6] = 0; R[7] = 0;  
	R[8] = 0; R[9] =  0; R[10] = pow(0.05*sqrt(freq), 2); R[11] = 0;
	R[12] = 0; R[13] = 0; R[14] = 0; R[15] = pow(0.05*sqrt(freq), 2);
#endif 
#ifdef dev
	printf("R:\n");
	Print2dMatr(R, getDimZ(), getDimZ());
#endif
	
	alglib::real_1d_array H_Papr_H_tR;
	H_Papr_H_tR.setlength(getDimZ() * getDimZ());
	for (int iii=0; iii<getDimZ()*getDimZ(); ++iii)
		H_Papr_H_tR[iii] = H_Papr_H_t[iii] + R[iii];

#ifdef dev
// #if 1
	printf("H_Papr_H_tR:\n");
	Print2dMatr(H_Papr_H_tR, getDimZ(), getDimZ());
#endif	

	alglib::real_2d_array invHPR; //обратная матрица к H_Papr_H_tR; 2D потому что либо библиотечная функция так считает, либо сингулярка, которая тоже использует библиотечные функции умножения
	alglib::real_2d_array HPR; //обратная матрица к H_Papr_H_tR; 2D потому что либо библиотечная функция так считает, либо сингулярка, которая тоже использует библиотечные функции умножения
	invHPR.setlength(getDimZ(), getDimZ());
	HPR.setlength(getDimZ(), getDimZ());
	//делаем копию
	for (int iii=0; iii<getDimZ(); ++iii)
		for (int jjj=0; jjj<getDimZ(); ++jjj)
			HPR[iii][jjj] = H_Papr_H_tR[index_3(getDimZ(), iii, jjj)];

	#if 0
	try
	{
		alglib::matinvreport rep; // по примеру
		alglib::rmatrixinverse(invHPR, rep); //перезапишется матрица invHPR
	}
	catch(alglib::ap_error alglib_exception)
	{
		printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
		return ;
	}
	#else
	pinv(HPR, invHPR);
	#endif
	//Переписываем 2D в 1D
	alglib::real_1d_array inv_H_P_H_t;
	inv_H_P_H_t.setlength(getDimZ()*getDimZ());
	for(int iii=0; iii< getDimZ(); ++iii)
		for(int jjj=0; jjj< getDimZ(); ++jjj)
			inv_H_P_H_t[index_3(getDimZ(), iii, jjj)] = invHPR[iii][jjj];
	
#ifdef dev
// #if 1
	printf("inv_H_P_H_t after inverse:\n");
	Print2dMatr(inv_H_P_H_t, getDimZ(), getDimZ());
#endif
	
	//Наконец-то вычислили коэффициент усиления
	matMul(getDimX(), getDimZ(), getDimZ(), Papr_H_t, inv_H_P_H_t, K);

#ifdef dev
	printf("K:\n");
	Print2dMatr(K, getDimX(), getDimZ());
#endif

	alglib::real_1d_array Kv;
	Kv.setlength(getDimX() * 1);
	matMul(getDimX(), getDimZ(), 1, K, v, Kv);

#ifdef dev
	printf("Kv:\n");
	Print2dMatr(Kv, getDimX(), 1);
#endif

	//поправляем наше предсказание (экстраполяцию)
	for (int iii=0; iii< getDimX(); ++iii)
		x[iii] = x_1[iii] + Kv[iii];

	//Подготовка к вычислению апостериорной ошибки оценивания
	alglib::real_1d_array KH;
	KH.setlength(getDimX() * getDimX());
	matMul(getDimX(), getDimZ(), getDimX(), K, H, KH);

	alglib::real_1d_array IKH;
	IKH.setlength(getDimX()* getDimX());
	for(int iii=0; iii<getDimX()*getDimX(); ++iii)
		IKH[iii] = I[iii] - KH[iii];

	//Вычисляем апостериорную ошибку оценивания
	matMul(getDimX(), getDimX(), getDimX(), IKH, Papr, Papst);
	++iter;
}

int AdaptiveRKalman::getDimX() //функция для получения private размерности
{
	return dim_x;
}
int AdaptiveRKalman::getDimZ() //функция для получения private размерности
{
	return dim_z;
}
int AdaptiveRKalman::getDimQ() //функция для получения private размерности
{
	return dim_q;
}
void AdaptiveRKalman::setDimX(int val) //функция для установки private размерности
{
	dim_x = val;
}
void AdaptiveRKalman::setDimZ(int val) //функция для установки private размерности
{
	dim_z = val;
}

void AdaptiveRKalman::setDimQ(int val) //функция для установки private размерности
{
	dim_q = val;
}

void AdaptiveRKalman::Print2dMatr(alglib::real_1d_array A, int dim1, int dim2)//Функция дл вывода на печать матриц
{
	for (int iii=0; iii<dim1; ++iii)
	{
		for (int jjj=0; jjj<dim2; ++jjj) printf("%e ", A[index_3(dim2, iii, jjj)]);
		printf("\n");
	}
}

void AdaptiveRKalman::Reset()//Сброс фильтра Калмана
{
	for (int i =0; i < dim_x; i++)
	{
		for (int j =0; j < dim_x; j++)
		{
		}
		x[i] = 0;
	}
}