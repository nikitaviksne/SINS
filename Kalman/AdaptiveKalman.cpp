#include "ap.h"
#include "linalg.h"
#include "AdaptiveKalman.h"
#include "mathematics_kalman.h"
#include "std.h"

AdaptiveKalman::AdaptiveKalman(int dimx, int dimz)
{
	setDimX(dimx);
	setDimZ(dimz);

	Papr.setlength(dim_x * dim_x);
	A.setlength(dim_x * dim_x);
	Phi.setlength(dim_x * dim_x);
	#if 0
	Phi[index_3(dim_x, 0, 0)] = 1;	Phi[index_3(dim_x, 0, 1)] = 1;
	Phi[index_3(dim_x, 1, 0)] = 0;	Phi[index_3(dim_x, 1, 1)] = 1;
	#endif

	Papst.setlength(dim_x * dim_x);

	x.setlength(dim_x * 1);
	x_1.setlength(dim_x * 1);
	Q.setlength(dim_x * dim_x);
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
}
void AdaptiveKalman::Init(Ldoub* initVal, Ldoub* q, Ldoub* h)
{//Инициализация значениями


	//Инициализация начальными значениями
	for (int iii=0; iii<getDimX(); ++iii)
	{
		x[iii] = initVal[iii];
	}
	for(int iii=0; iii<getDimX(); ++iii)
		for(int jjj=0; jjj<getDimX(); ++jjj)
			Q[index_3(dim_x, iii, jjj)] = q[index_3(dim_x, iii, jjj)];
	for (int iii=0; iii<getDimZ(); ++iii)
		for(int jjj=0; jjj<getDimX(); ++jjj)
			H[index_3(dim_x, iii, jjj)] = h[index_3(dim_x, iii, jjj)];

	iter = 0;
	//Флаг инициализации
	init=true;
}
void AdaptiveKalman::Predict()
{
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

	//Вот тут возможно должна быть матрица входного шума G, но я тогда не понимаю, получается, что тогда измерительный шум один и тот-же?
	for (int iii=0; iii<getDimX()*getDimX(); ++iii)
		Papr[iii] = Phi_Papst_Phi_t[iii] + Q[iii];

}

void AdaptiveKalman::Update(Ldoub* zin/*измерения обычные C-массивы*/)
{
	for(int i=0; i< getDimZ(); i++)
		this->z[i] = zin[i];
	//вычисляем различия между измерениями и предсказанием
	alglib::real_1d_array Hx_1;
	Hx_1.setlength(getDimZ() * 1);
	matMul(getDimZ(), getDimX(), 1, H, x_1, Hx_1);
	for (int iii=0; iii<getDimZ(); ++iii)
		v[iii] = this->z[iii] - Hx_1[iii];

	//Вычисляем v*v.T
	alglib::real_1d_array vvt;
	vvt.setlength(getDimZ() * getDimZ());
	matMul(getDimZ(), 1, getDimZ(), v, v, vvt);

	//Вычислям C
	for(int iii=0; iii<getDimZ(); iii++) 
			for(int jjj=0; jjj<getDimZ(); ++jjj)
			{
				Ldoub buf;
				buf = iter / (iter + 1.) * (C[index_3(getDimZ(), iii, jjj)]) + 1. / (iter + 1.) *  (vvt[index_3(getDimZ(), iii, jjj)]); // для стационарной системы
				C[index_3(getDimZ(), iii, jjj)] = buf;
			}

	alglib::real_1d_array H_t;
	H_t.setlength(getDimX() * getDimZ());
	transpose(getDimZ(), getDimX(), H, H_t);

	alglib::real_1d_array Papr_H_t; //потом еще понадобится
	Papr_H_t.setlength(getDimX() * getDimZ());//временная переменная H @ Papr
	matMul(getDimX(), getDimX(), getDimZ(), Papr, H_t, Papr_H_t);

	alglib::real_1d_array H_Papr_H_t;
	H_Papr_H_t.setlength(getDimZ() * getDimZ());//временная переменная H @ Papr @ H.T
	matMul(getDimZ(), getDimX(), getDimZ(), H, Papr_H_t, H_Papr_H_t);

	//Вычисляем  R
	for(int iii=0; iii<getDimZ()*getDimZ(); ++iii)
	{
		Ldoub temp = C[iii] - H_Papr_H_t[iii];		
		R[iii] = temp;
	}
	
	for (int iii=0; iii<getDimZ(); ++iii)
		if (R[index_3(getDimZ(), iii, iii)] < 0) //проверка на отрицательность диагональных элементов C
		{
			for(int jjj=0; jjj<getDimZ(); ++jjj) //все диагональные элементы устанавливаем в 0
				R[index_3(getDimZ(), jjj, jjj)] = 0;
			break; //и выходим из внешнего цикла
		}
	

	alglib::real_1d_array H_Papr_H_tR;
	H_Papr_H_tR.setlength(getDimZ() * getDimZ());
	for (int iii=0; iii<getDimZ()*getDimZ(); ++iii)
		H_Papr_H_tR[iii] = H_Papr_H_t[iii] + R[iii];

	alglib::real_2d_array invHPR; //обратная матрица к H_Papr_H_tR; 2D потому что библиотечная функция
	invHPR.setlength(getDimZ(), getDimZ());
	//делаем копию
	for (int iii=0; iii<getDimZ(); ++iii)
		for (int jjj=0; jjj<getDimZ(); ++jjj)
			invHPR[iii][jjj] = H_Papr_H_tR[index_3(getDimZ(), iii, jjj)];
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
	//Переписываем 2D в 1D
	alglib::real_1d_array inv_H_P_H_t;
	inv_H_P_H_t.setlength(getDimZ()*getDimZ());
	for(int iii=0; iii< getDimZ(); ++iii)
		for(int jjj=0; jjj< getDimZ(); ++jjj)
			inv_H_P_H_t[index_3(getDimZ(), iii, jjj)] = invHPR[iii][jjj];
	//Наконец-то вычислили коэффициент усиления
	matMul(getDimX(), getDimZ(), getDimZ(), Papr_H_t, inv_H_P_H_t, K);

	alglib::real_1d_array Kv;
	Kv.setlength(getDimX() * 1);
	matMul(getDimX(), getDimZ(), 1, K, v, Kv);
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

int AdaptiveKalman::getDimX() //функция для получения private размерности
{
	return dim_x;
}
int AdaptiveKalman::getDimZ() //функция для получения private размерности
{
	return dim_z;
}
void AdaptiveKalman::setDimX(int val) //функция для установки private размерности
{
	dim_x = val;
}
void AdaptiveKalman::setDimZ(int val) //функция для установки private размерности
{
	dim_z = val;
}
