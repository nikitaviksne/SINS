#ifndef MATHEMATICS_H
#define MATHEMATICS_H

#include "ap.h"
#include "linalg.h"

//шаблон функции для умножения любых типов одномерных массивов (и alglib и обычных C-массивов)
template <class T> void matMul(int aRows, int aCols, int bCols, T A, T B, T& C);

//Перегрузка функции для умножения одномерных массивов alglib как матрицы
void matMul(int aRows, int aCols, int bCols, alglib::real_1d_array A, alglib::real_1d_array B, alglib::real_1d_array& C);
//void matMul(int);

void transpose(int, alglib::real_1d_array, alglib::real_1d_array& ); //для квадратных матриц

void transpose(int, int, alglib::real_1d_array, alglib::real_1d_array& ); //для прямоугольных матриц

void pinv(alglib::real_2d_array , alglib::real_2d_array& );

void diag(alglib::real_1d_array, alglib::real_2d_array&);

#endif //MATHEMATICS_H
