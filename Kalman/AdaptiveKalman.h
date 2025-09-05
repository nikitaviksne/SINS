#ifndef ADAPTIVEKALMAN_H
#define ADAPTIVEKALMAN_H


#include "ap.h"
#include "linalg.h"
#include "std.h"

class AdaptiveKalman
{
public:
	AdaptiveKalman(int, int, Ldoub h = 1/*такт в секундах, по умолчанию*/);
	void Init(Ldoub*, Ldoub*, Ldoub* h);
	void Predict(); // функция предсказания по модели
	void Update(Ldoub* );// функция обновления, т.е. вычисление оценки вектора по измерениям

	alglib::real_1d_array Phi;
	alglib::real_1d_array A;
	alglib::real_1d_array A2; //матрица A в квадрате
	alglib::real_1d_array A3; //матрица A в кубе
	alglib::real_1d_array H; //матрица измерений
	alglib::real_1d_array G; // матрица входного шума
	alglib::real_1d_array Q; //ковариационная матрица входных значений (модели)
	alglib::real_1d_array R; // ковариационная матрица вектора измерений (входная)
	alglib::real_1d_array C; // ковариационная матрица вектора измерений (входная)
	alglib::real_1d_array Papr;//Априорная ошибка оценивания
	alglib::real_1d_array Papst;//Апостериорная ошибка оценивания
	alglib::real_1d_array K; //Коэффициент усиления
	alglib::real_1d_array x;//вектор состояний
	alglib::real_1d_array x_1; // предсказание вектора состояния
	alglib::real_1d_array z; //вектор измерений
	alglib::real_1d_array v; //вектор невязки

	alglib::real_1d_array I; //единичная матрица, чтобы задать один раз и все
	bool init; //Флаг инициализации	
	int iter; //счетчик тактов (для вычисления C_k)

	int getDimX(); //функция для получения private размерности
	int getDimZ(); //функция для получения private размерности
	void setDimX(int); //функция для установки private размерности
	void setDimZ(int); //функция для установки private размерности
	void Print2dMatr(alglib::real_1d_array A, int dim1, int dim2);//Функция дл вывода на печать матриц
	Ldoub h;
private:
	int dim_x; //размер вектора состояния
	int dim_z; //размер вектора измерений

};


#endif //ADAPTIVEKALMAN_H
