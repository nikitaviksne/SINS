#include <stdio.h>
#include "std.h"
#include "matrix.h"
#include "mathematics.h"
#include <cmath>
#include <fstream> //файловые потоки для чтения бинарного файла библиотека C++
//#include <QDataStream>
//#include <QFile>
#include "alignment.h"
#include "SolveOrient.h"
#include "SolveNav.h"
#include "AdaptiveRKalman.h" //Для адаптивного фильтра Калмана
#include "UsualKalman.h" //Для обычного фильтра Калмана
#include "mathematics_kalman.h"
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include "quaternions.h"

#define dev

#ifdef dev
#include <stdio.h>
#endif

void MatrOB(Ldoub H, Ldoub R, Ldoub P, Ldoub* C, int size)
{
	C[0] = (Ldoub) cos(R) * cos(H) + sin(R)*sin(H)* sin(P);
	C[1] = (Ldoub) sin(R) * cos(H)*sin(P) - cos(R)*sin(H);
	C[2] = (Ldoub) -sin(R) * cos(P);
	C[3] = (Ldoub) cos(P) * sin(H);
	C[4] = (Ldoub) cos(H) * cos(P);
	C[5] = (Ldoub) sin(P);
	C[6] = (Ldoub) sin(R) * cos(H) - cos(R) * sin(H) * sin(P);
	C[7] = (Ldoub) -sin(R) * sin(H) - cos(R) * cos(H) * sin(P);
	C[8] = (Ldoub) cos(R) * cos(P);
}

bool ReadFile(std::ifstream &is, bool AllowBiasAcc, bool AllowBiasGyr, bool AllowRandAcc, bool AllowRandGyr, bool AllowRandVgps, Ldoub* Ab, Ldoub* Omb, Ldoub* Vgps)
{//Чтение файла с помощью fstream
	// Чтение из файла ускорений и угловых скоростей
	if (is.read((char *) &Ab[0], sizeof(double))) //пока считывается первое измерение (ускорение по оси Xb)
	{
		//Заканчиваю считывать показания акселерометров
		for (int jjj=1; jjj<3; ++jjj)
			is.read((char *) &Ab[jjj], sizeof(double));

		//Считываю показания ДУС
		for (int jjj=0; jjj<3; ++jjj)
			is.read((char *) &Omb[jjj], sizeof(double));

		Ldoub BiasAb[3] = {0}; // Постоянные погрешности акселерометров
		Ldoub BiasOmb[3] = {0}; // Постоянные погрешности гироскопов
		Ldoub RandAb[3] = {0}; // Случайные погрешности акселерометров
		Ldoub RandOmb[3] = {0};// Случайные погрешности гироскопов

		//Считываю показания постоянных смещений нуля акселерометров
		for (int jjj=0; jjj<3; ++jjj)
			is.read((char *) &BiasAb[jjj], sizeof(double));

		//Считываю показания постоянных смещений нуля ДУС
		for (int jjj=0; jjj<3; ++jjj)
			is.read((char *) &BiasOmb[jjj], sizeof(double));

		//Считываю показания случайных смещений нуля акселерометров
		for (int jjj=0; jjj<3; ++jjj)
			is.read((char *) &RandAb[jjj], sizeof(double));

		//Считываю показания случайных смещений нуля ДУС
		for (int jjj=0; jjj<3; ++jjj)
			is.read((char *) &RandOmb[jjj], sizeof(double));

		//добавление дрейфов к показаниям
		for (int ii=0; ii<3; ++ii)
		{
			Ab[ii] += AllowBiasAcc*BiasAb[ii] + AllowRandAcc*RandAb[ii];
			Omb[ii] += AllowBiasGyr*BiasOmb[ii] + AllowRandGyr*RandOmb[ii];
		}

		//Считываем идеальные скорости от GPS
		for (int iii =0; iii<2; ++iii)
			is.read((char *) &Vgps[iii], sizeof(double));

		// Считываем шум скорости от GPS (с нулевым средним)
		Ldoub RandomVgps[2] = {0};
		for (int iii =0; iii<2; ++iii)
		{
			is.read((char *) &RandomVgps[iii], sizeof(double));
			//Добавление шума к показаниям по скорости
			Vgps[iii] += (Ldoub) AllowRandVgps*RandomVgps[iii];
		}
	}
	else
	{//формирую признак конца файла
		return false;
	}
	return true;
}

int Readfile(FILE* is, int numRes /*должное количество считанных данных из файла*/, bool AllowBiasAcc, bool AllowBiasGyr, bool AllowRandAcc, bool AllowRandGyr, bool AllowRandVgps, bool AllowRandCoogps, Ldoub* Ab, Ldoub* Omb, Ldoub* Vgps, Ldoub* CooGps, Ldoub* V0 /*всегда идеальные скорости (без погрешностей GPS)*/)
{
	Ldoub BiasAb[3] = {0}; // Постоянные погрешности акселерометров
	Ldoub BiasOmb[3] = {0}; // Постоянные погрешности гироскопов
	Ldoub RandAb[3] = {0}; // Случайные погрешности акселерометров
	Ldoub RandOmb[3] = {0};// Случайные погрешности гироскопов
	Ldoub RandomVgps[2] = {0}; //Случайные погрешности по скоростям GPS
	Ldoub RandomCoogps[2] = {0}; //Случайные погрешности по координатам GPS

	int res = fscanf(is, "%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;%Lf;\n", &Ab[0], &Ab[1], &Ab[2], &Omb[0], &Omb[1], &Omb[2], &BiasAb[0], &BiasAb[1], &BiasAb[2], &BiasOmb[0], &BiasOmb[1], &BiasOmb[2], &RandAb[0], &RandAb[1], &RandAb[2], &RandOmb[0], &RandOmb[1], 	&RandOmb[2], &V0[0], &V0[1], &RandomVgps[0], &RandomVgps[1], &CooGps[0], &CooGps[1], &RandomCoogps[0], &RandomCoogps[1]);
	if (res == numRes)
	{
		//добавление дрейфов к показаниям инерциальных датчиков
		for (int ii=0; ii<3; ++ii)
		{
			Ab[ii] += AllowBiasAcc*BiasAb[ii] + AllowRandAcc*RandAb[ii];
		}
		Omb[0] += 1.*(AllowBiasGyr*BiasOmb[0] + AllowRandGyr*RandOmb[0]);
		Omb[1] += 1.*(AllowBiasGyr*BiasOmb[1] + AllowRandGyr*RandOmb[1]);
		Omb[2] += 0.*(AllowBiasGyr*BiasOmb[2] + AllowRandGyr*RandOmb[2]);
		//добавление шума к показаниям СНС
		for(int iii=0; iii<2; ++iii)
		{
			Vgps[iii] = V0[iii] + (Ldoub) AllowRandVgps*RandomVgps[iii]; //Скорости
			CooGps[iii] += (Ldoub) AllowRandCoogps*RandomCoogps[iii]; //Скорости
		}

		return res;
	}
	else return -1; //означает что что-то не так
}

int Readfile(FILE* is, Ldoub* Ab, Ldoub* Omb) // Чтение сырых данных (от микромеханики)
{
	int timestamp;
	int res = fscanf(is, "%d;%lf;%lf;%lf;%lf;%lf;%lf\n", &timestamp, &Ab[0], &Ab[1], &Ab[2], &Omb[0], &Omb[1], &Omb[2]);
	if (res == 6)
	{
		return res;
	}
	else return -1; //означает что что-то не так
}

#if 0 //QT
void ReadFile(QDataStream &in, bool AllowBiasAcc, bool AllowBiasGyr, bool AllowRandAcc, bool AllowRandGyr, bool AllowRandVgps, Ldoub* Ab, Ldoub* Omb, Ldoub* Vgps)
{

	//Считываю показания акселерометров
	for (int jjj=0; jjj<3; ++jjj)
		in >> Ab[jjj];

	//Считываю показания ДУС
	for (int jjj=0; jjj<3; ++jjj)
	{
		in >> Omb[jjj];
		//Omb[jjj] *= (-1.);
	}
	Ldoub BiasAb[3] = {0}; // Постоянные погрешности акселерометров
	Ldoub BiasOmb[3] = {0}; // Постоянные погрешности гироскопов
	Ldoub RandAb[3] = {0}; // Случайные погрешности акселерометров
	Ldoub RandOmb[3] = {0};// Случайные погрешности гироскопов
	//Считываю показания постоянных смещений нуля акселерометров
	for (int jjj=0; jjj<3; ++jjj)
		in >> BiasAb[jjj];

	//Считываю показания постоянных смещений нуля ДУС
	for (int jjj=0; jjj<3; ++jjj)
		in >> BiasOmb[jjj];

	//Считываю показания случайных смещений нуля акселерометров
	for (int jjj=0; jjj<3; ++jjj)
		in >> RandAb[jjj];

	//Считываю показания случайных смещений нуля ДУС
	for (int jjj=0; jjj<3; ++jjj)
		in >> RandOmb[jjj];
	//добавление дрейфов к показаниям
	for (int ii=0; ii<3; ++ii)
	{
		Ab[ii] += (Ldoub) AllowBiasAcc*BiasAb[ii] + (Ldoub) AllowRandAcc*RandAb[ii];
		Omb[ii] +=(Ldoub) AllowBiasGyr*BiasOmb[ii] + (Ldoub) AllowRandGyr*RandOmb[ii];
	}

	//Считываем идеальные скорости от GPS
	for (int iii =0; iii<2; ++iii)
		in >> Vgps[iii];

	// Считываем шум скорости от GPS (с нулевым средним)
	Ldoub RandomVgps[2] = {0};
	for (int iii =0; iii<2; ++iii)
	{
		in >> RandomVgps[iii];
		//Добавление шума к показаниям по скорости
		Vgps[iii] += (Ldoub) AllowRandVgps*RandomVgps[iii];
	}

}
#endif

#if 0
void GeneratedSens(Ldoub *Ab, Ldoub *Omb, Ldoub Vabs, Ldoub H0, int cur_time, int t_alignment, Ldoub U, Ldoub g, Ldoub* Cnb)
{
	/*Генерирование (моделирование) показаний ч.э*/
	Ldoub Omo[3] = {0};
	Ldoub Ao[3] = {0};
	if (cur_time <= t_alignment)
	{
		Ao[0] = 0.0; Ao[1] = 0.0; Ao[2] = g;
		Omo[0] = 0;
		Omo[1] = (Ldoub) U*cos(phi0);
		Omo[2] = (Ldoub) U*sin(phi0);
	}
	else
	{
		/*
		Omo[0] = (Ldoub) -Vabs*cos(H0) / (Rphi + 0);
		Omo[1] = (Ldoub) Vabs*sin(H0) / ((Rlambda + 0)) + (Ldoub) U*cos(phi0);
		Omo[2] = (Ldoub) Vabs*sin(H0) * tan(phi0) / (Rlambda + 0) + (Ldoub) U*sin(phi0);
		*/
		Omo[0] = 0;
		Omo[1] = (Ldoub) U*cos(phi0);
		Omo[2] = (Ldoub) U*sin(phi0);
		Ao[0] = 0.0;// + (Ldoub) ( Omo[1]*0 -(Ldoub) Omo[2]*Vabs*cos(H0) + (Ldoub) U*cos(phi0)*0 - (Ldoub) U*sin(phi0)*Vabs*cos(H0));
		Ao[1] = 0.0;// +(Ldoub) (-Omo[0]*0 +(Ldoub) Omo[2]*Vabs*sin(H0) + (Ldoub) U*sin(phi0)*Vabs*sin(H0));
		Ao[2] = g;
		phi0 += (Ldoub) Vabs*cos(H0)/(R + 0)/freq1;
	}
	MulMatrD(Cnb, Ao, Ab, 3, 3, 1); // проекция ускорений на связанные оси
	MulMatrD(Cnb, Omo, Omb, 3, 3, 1); // проекция угловых скоростей на связанные оси
}
#endif
int main(int argc, char *argv[])
{
	//инициализация необходимых переменных и констант
	const Ldoub g = 9.81;
	const Ldoub a = 6378245.;
	const Ldoub b = 6356856.;
	const Ldoub e = sqrt(1 - b*b/a/a);
	const Ldoub R = 6400e3;
	Ldoub nu2 = sqrt(g/R); // квадрат частоты Шуллера
	const Ldoub pi = 3.141592653589793;
	const Ldoub U = 7.27220521664304e-05;
	const Ldoub rad2deg = 180./M_PI; // из градусов в час в радианы в секунду
	const Ldoub deg2rad = 1./rad2deg;
	Ldoub freq1 = 400.; // частота свехбыстрого цикла
	Ldoub h1 (1./freq1); // период дискретизации свехбыстрого цикла
	Ldoub freq = freq1/4.;//100; // частота измерений с инерциальных датчиков
	Ldoub h (1./freq); //период дискретизации
	int t_nav = 180*60; // время работы нав алгоритма в секундах
	int t_alignment = (int) 5*60*freq; // время выставки в тактах
	Ldoub Rlambda;
	Ldoub Rphi;
	Ldoub phi0 = (Ldoub) 55*deg2rad;// и для моделирования
	int cur_time = 0; // текущий такт!! измерения
	// Для моделирования показаний Ч.Э.
	Ldoub H0 = (Ldoub) (strtod(argv[1], NULL))*deg2rad;
	Ldoub P0 = (Ldoub) (strtod(argv[2], NULL))*deg2rad;
	Ldoub R0 = (Ldoub) (strtod(argv[3], NULL))*deg2rad;
	Ldoub Vabs = (Ldoub) strtod(argv[4], NULL);

	//printf("H0 = %.4f P0 = %.4f R0 = %.4f Vabs = %.4f\n", H0, P0, R0, Vabs);

	Ldoub Vgps[3] = {0};
	Ldoub CooGps[3] = {0};
	Ldoub Cnb[9];
	MatrOB(H0, R0, P0, Cnb, 3); // матрица перехода из опорной в связанную
	//Необходимое для выставки
	Ldoub lambda0 = (Ldoub) 33*deg2rad;
	Ldoub DeltaHeading = 0, DeltaRoll = 0, DeltaPitch = 0;// ошибки выставки по курсу, крену и тангажу соответственно
	Ldoub Heading = 0, Roll = 0, Pitch = 0;
	bool AlignmentContinue = true; // для начала выставки
	// Необходимые массивы для решение навигационной задачи
	Ldoub Ab[3]={0}; // Ускорения в связанных осях
	Ldoub Ao[3] = {0}; // Ускорения в опорныых осях
	Ldoub Wo[3] = {0}; //Приращение скоростей в опорных осях
	Ldoub Omb[3] = {0}; // Угловые скорости в связанных осях
	//Ldoub Oms[3] = {0}; // угловые скорости от линейного движения + Земля
	Ldoub MeanAb[3] = {0};
	Ldoub MeanOmb[3] = {0};
	Ldoub StdAb[3] = {0};
	Ldoub StdOmb[3] = {0};
	Ldoub Cbn[9] = {0};
	Ldoub Cib[9] = {0}; //матрица перехода из инерциальной системы в связанную. Начальное значение равно транспонированной матрицы на момент окончания выставки
	quaternion Qf; //медленный (филнальный) кватернион. Быстрый в функции решение задачи ориентации
	Ldoub Cin[9] = {1.,0,0,0,1.,0,0,0,1.}; //матрица перехода из инерцальной в опорную. Начальное знвчение -- единичная Cin(0)=E
	Ldoub V0[3] = {(Ldoub) Vabs*sin(H0), (Ldoub) Vabs*cos(H0), 0}; // линейные скорости E; N; Up
	Ldoub Err_V[3] = {0}; // ошибки по скоростям
	Ldoub ErrVins[4] = {0}; // разница ошибок между ИНС и GPS (для корректирующих поправок) (широта, долгота, Vx, Vy)
	// массивы для выходных значений
	Ldoub V[3] = {(Ldoub) V0[0], (Ldoub) V0[1], 0}; // линейные скорости E; N; Up (на всякий случай, пока резерв)
	Ldoub Coordinates[3] = {phi0, lambda0, 0}; // Географические кординаты: широта, долгота и высота
	Ldoub CoordError[3] = {0}; // Ошибки в м (dE, dN)
	Ldoub Orientation[3] = {0}; // Углы ориентации

	Rlambda = (Ldoub) R/sqrt(1.-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
	Rphi = (Ldoub) R*(1. - e*e)/(sqrt(1.-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1.-e*e*sin(Coordinates[0])*sin(Coordinates[0])));

	Ldoub Omo[3] = {-V0[1]/(Rphi + Coordinates[2]), V0[0]/(Rlambda + Coordinates[2]), V0[0]/(Rlambda + Coordinates[2])*tan(Coordinates[0])}; // Угловые скорости в опорных осях
	Ldoub omo[3] = { Omo[0], Omo[1] +  U*cos(Coordinates[0]), Omo[2] + U*sin(Coordinates[0]) }; //переносные (вроде даже абсолютные) угловые скорости, вчисляются в решении задачи ориентации (SolveOrient)

	bool AllowBiasAcc = (bool) (strtod(argv[6], NULL));
	bool AllowBiasGyr = (bool) (strtod(argv[7], NULL));
	bool AllowRandAcc = (bool) (strtod(argv[8], NULL));
	bool AllowRandGyr = (bool) (strtod(argv[9], NULL));
	bool AllowRandVgps = (bool) (strtod(argv[10], NULL));
	bool AllowRandCoogps = (bool) (strtod(argv[11], NULL));
	//Создаем квазикоординаты
	Ldoub alpha[12] = {0}; //малые приращения углов 3 показания на 4 тактах (матрица 3*4)
	Ldoub w[12] = {0}; // малые приращения скоростей (матрица 3*4)

	bool derectNorm (true); //направление ортогонализации и нормализации

	Ldoub kcor1 (0.006138995628986877); // метод наименьшей дисперсии () программой на python
	Ldoub kcor2 (11.54438988378206); // метод наименьшей дисперсии () программой на python
	bool allowCorr (false); //разрешение на коррекцию
	// Чтение из файла ускорений и угловых скоростей


#if 0
	QFile file(argv[1]);
	file.open(QIODevice::ReadOnly);
	QDataStream in(&file);
	in.setByteOrder(QDataStream::LittleEndian);
#endif

#if 1
	int numRes (26); //количество переменных, котоорое должно быть считано из файла в безошибочном случае
	FILE* file=fopen(argv[5], "rt");
	fscanf(file, "%*s;");//чтение строки заголовка, она не нужна, выбрасываем
	int res = Readfile(file, numRes, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, AllowRandCoogps, Ab, Omb,Vgps, CooGps, V0 /*Всегда идеальные скорости*/);
#else
	std::ifstream in(argv[5], std::ios::binary);// !in.eof() //условие цикла while для бинарного файла
	#define FILE_STREAM
#endif
	Ldoub resultV[2] = {0.};
	Ldoub Verr1[2] = {0.}; // ошибки интегрирования ускорений
	Ldoub Verr2[2] = {0.}; // ошибки накопления скоростей
	///*
	const int dim_state (10 + 2); //размер вектора состояния (+ координаты GPS)
	const int dim_sense (3 + 2); //размер вектора измерения (+ координаты GPS)
	const int dim_input_noise (3); //размер матрицы входных шумов Q
	//*/
	UsualKalman filter(dim_state, dim_sense, dim_input_noise, h/*шаг дискретизации*/); //Создаю объект фильтра Калмана с матрицей размера 6*6 и измерениями 2*1 (вертикальную скорость не учитываю)
	Ldoub x0[dim_state] = {0}; //Начальные оценочные значения дрейфов
	Ldoub H[dim_sense*dim_state] = {0}; //матрица наблюдения
	#if 1 // вектор состояния 6
	H[index_3(dim_state, 0, 0)] = 1.; //varphi широта
	H[index_3(dim_state, 1, 1)] = 1.; //lambda долгота
	H[index_3(dim_state, 2, 2)] = 1.; //Ve
	H[index_3(dim_state, 3, 3)] = 1.; //Vn
	H[index_3(dim_state, 4, 4)] = 1.; //Vup
	print2dMatr(H, dim_sense, dim_state);

	//Матрца ковариации входных шумов (модели)
	Ldoub q[dim_input_noise * dim_input_noise] = {0};
	q[index_3(dim_input_noise, (dim_input_noise - 3), (dim_input_noise - 3))] = 0*pow(0.02*deg2rad/3600., 4); //delta_omega_x
	q[index_3(dim_input_noise, (dim_input_noise - 2), (dim_input_noise - 2))] = 0*pow(0.02*deg2rad/3600., 4); //delta_omega_y
	q[index_3(dim_input_noise, (dim_input_noise - 1), (dim_input_noise - 1))] = 0*pow(0.02*deg2rad/3600., 2); //delta_omega_z
	
	// Заполнение всей матрицы G нулями
	for (int iii=0; iii<filter.getDimX(); iii++)
			for(int jjj=0; jjj<filter.getDimQ(); jjj++)
				filter.G[index_3(filter.getDimQ(), iii, jjj)] = 0;
	//printf("Для сравнения index_3(dim_state-2) = %d, index_3(7) = %d, dim_state - 2 (%d-2) = %d\n", index_3(dim_state, dim_state-2, dim_state - 2), index_3(dim_state, 7, 7), dim_state, dim_state-2);
	//printf("Матрица q\n");
	//print2dMatr(q, dim_state, dim_state);
	#else // вектор состояния 3
	H[index_3(dim_state, 0, 0)] = 1;
	//Матрца ковариации входных шумов (модели)
	Ldoub q[dim_state*dim_state] = {0};
	q[index_3(dim_state, 0, 0)] = 1e-19;
	#endif
	Ldoub r[dim_sense*dim_sense] = {pow(0.2, 2), 0, 0, 0, 0,
									0, pow(0.2, 2), 0, 0, 0, 
									0, 0, pow(0.05, 2), 0, 0,
									0, 0, 0, pow(0.05, 2), 0,
									0, 0, 0, 0, pow(0.05, 2)};

#if 0
	//Для лучшей обусловленности матрицы HPH_t увеличиваю начальные значения априорной ошибки
	for (int iii=0; iii<filter.getDimX()*filter.getDimX(); ++iii)
	{
		filter.Papr[index_3(filter.getDimX(), iii, iii)] = 1;
	}
#endif

	while(res == numRes ) // Пока возможно чтение из файла. Для бинарного файла -- !in.eof() //условие цикла while для бинарного файла //  Readfile(file, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps) == 22 --- условие цикла While для текстового файла
	{ // цилк while для частоты 100 Гц, поэтому нельзя сделать общее чтение файла для выставки и навигации

		// этап выставки
		if ((cur_time <= t_alignment) && AlignmentContinue )
		{
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
			#ifdef FILE_STREAM //если определен файловый поток, то чтение из бинарника, читаем тут
			ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps); //чтение из бинарного файла данных используемых для выставки
			#else
			res = Readfile(file, numRes, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, AllowRandCoogps, Ab, Omb,Vgps, CooGps, V0);
			#endif
			
			alignment(Ab, Omb, MeanAb, MeanOmb, StdAb, StdOmb, cur_time, g, U, phi0, Cbn, &Qf);
		
			++cur_time; // для 400 Гц
			/*
			printf("Mean Omb\n");
			print2dMatr(MeanOmb, 1, 3);
			*/
			continue;
		}
		else
		{
			if (AlignmentContinue)
			{
				for (int i=0; i<3; i++)
				{
					for (int j=0; j<3; j++)
					{
						Cbn[index_3(3, i, j)] = 0;
					}
					Cbn[index_3(3, i,i)] = 1;
				}
				printf("Cbn перед переводом из кватерниона\n");
				printf("[");
				for (int iii=0; iii<9; ++iii)
					printf("%.10Lf, ", Cbn[iii]);
				printf("]\n");
				Qf.print();
				Quat2Matr(Qf, Cbn);//проверка правильности работы функции

				printf("Cbn после перевода из кватерниона\n");
				printf("[");
				for (int iii=0; iii<9; ++iii)
					printf("%.10Lf, ", Cbn[iii]);
				printf("]\n");

				printf("Mean Omb not in Mean\n");
				print2dMatr(MeanOmb, 1, 3);
				Ldoub c0 = (Ldoub) sqrt(Cbn[index_3(3, 2, 0)]* Cbn[index_3(3, 2, 0)] + Cbn[index_3(3, 2, 2)]*Cbn[index_3(3, 2, 2)]);
				// Вычисление углов ориентации
				Heading = (Ldoub) atan2(Cbn[index_3(3, 0, 1)], Cbn[index_3(3, 1, 1)]);
				Roll = (Ldoub) - atan2(Cbn[index_3(3, 2, 0)], Cbn[index_3(3, 2, 2)]);
				Pitch = (Ldoub) atan2(Cbn[index_3(3, 2, 1)], c0);
				// вычисление ошибок выставки
				DeltaRoll = (StdAb[0] * MeanAb[2] - StdAb[2] * MeanAb[0])/(MeanAb[2]*MeanAb[2] + MeanAb[0]*MeanAb[0]);
				DeltaPitch = (StdAb[1])/sqrt(g*g - MeanAb[1]*MeanAb[1]);
				Ldoub DeltaWn[3];
				MulMatrD(Cbn, StdOmb, DeltaWn,3,3,1); // проекция дрейфов гироскопов на географические оси
				Ldoub DeltaAn[3];
				MulMatrD(Cbn, StdAb, DeltaAn,3,3,1); // проекция дрейфов акселерометров на географические оси
				DeltaHeading = - DeltaWn[0]/(U*cos(phi0)) + DeltaAn[0]/g*tan(phi0) - StdAb[2]/2/g*sin(2*Heading);
				printf("Alignment\n");
				printf("Heading (degrees) = %.30f Error = %.Lf\n", Heading*rad2deg, DeltaHeading*rad2deg);
				printf("roll (degrees) = %.30f Error = %.Lf\n", Roll*rad2deg, DeltaRoll*rad2deg);
				printf("pitch (degrees) = %.30f Error = %.Lf\n", Pitch*rad2deg, DeltaPitch*rad2deg);
				AlignmentContinue = false;
				Transpose2M(Cbn, Cib, 3); // начальное значение Cib; Cib(0)
				//Transpose(Cib,3);
				//printf("pow(5,2) = %.Lf\n", pow(5,2));
				cur_time = 0; // для 100 Гц

				Ldoub gravity[3] = {0, 0, g}; //ускорения силы тяжести в проекциях на оси опорной системы координат (географической)
				Ldoub ErrAcc[3] = {0}; //Ошибки акселерометров (нужны для невыставки)
				MulMatrD(Cib, gravity, ErrAcc,3,3,1); //здесь ErrAcc как временная матрица, а Cib=(Cbn)^t в начальный момент времени
				for (int iii=0; iii<3; ++iii)
					ErrAcc[iii] -= MeanAb[iii]; //теперь ErrAcc есть ошибки акселерометров
#if 0
				//начальные значения ошибок ориентации (для вектора состояния)
				x0[4] = -MeanAb[1]/g;//1e-4/g;
				x0[5] = MeanAb[0]/g;//-1e-4/g;
				x0[6] = 2.42407e-07 / (U*cos(phi0));//-(MeanOmb[0])/(U*cos(phi0)) + (ErrAcc[0]/g)*tan(phi0)//2.42407e-07 / (U*cos(phi0));//-(MeanOmb[0])/(U*cos(phi0));//azimuth misalignment
#endif
			}

		}
#if 1
		Ldoub MeanAlpha[3] = {0};//осредненные малые приращения углов (псевдокоординаты)
		Ldoub MeanW[3] = {0}; //осредненные малые приращения скоростей

		Ldoub Wp[3] = {0}; //проинтегрированные малые приращения. Начальные значения обнуляются на каждом такте быстрого цикла (с частотой 100 Гц)

		//Включение и выключение коррекции
	#if 0
		if (cur_time == int(60*60/h) )
			allowCorr = true;
		if (cur_time == int(70*60/h) )
			allowCorr = false;
	#endif
		res = Readfile(file, numRes, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, AllowRandCoogps, Ab, Omb,Vgps, CooGps, V0);
		/*Далее идет 100 Гц такт*/
		for (int iii=0; iii<3; iii++)
			Wp[iii] = Ab[iii]*h; //вычисляем малые приращения скорости вместо метода Рунге-Кутты, таим вот кустарным способом)
		/*Решение задачи навигации*/
		SolveNav(Wp, Ab, Cbn, Wo, Ao, V,  Coordinates, CoordError, Err_V, omo, h, Rphi, Rlambda, U, R, e, H0, V0, kcor1, ErrVins, allowCorr); //Err_V уже в этой функции вычисляется, поэтому я могу это значение использовать для коррекции
		//Решение задачи ориентации
		SolveOrient(Omb, &Qf, Cbn, Orientation, Coordinates, Omo, omo, V, phi0, Rphi, Rlambda, freq, h, U, cur_time, derectNorm, kcor2, ErrVins, allowCorr); //из одноименного заголовочного файла
		breakpiont: cur_time >= (284996 - 30000);
		#if 0
		V0[0] = Vabs*sin(Orientation[0]); //Не нужно, если из файла читаются идеальные скорости от GPS
		V0[1] = Vabs*cos(Orientation[0]); //Не нужно, если из файла читаются идеальные скорости от GPS
		#endif
		// инкремент тактов
		++cur_time;
		#if 0
		for (int iii=0; iii<3; ++iii)
		{
			for (int jjj=0; jjj<3; jjj++)
				printf("%.10Lf\t", Cbn[index_3(3, iii, jjj)]);
			printf("\n");
		}
		#endif
		/*
		По идее, куда-то сюда можно засунуть оценивание по Калману скоростей дрейфов гироскопов
		*/
	#if 1
	#if 1 //Модель для вектора состояния N>=6
		//Каждый такт пересчитываем матрицу A у фильтра Калмана
		//Delta dot varphi (широта)
		filter.A[index_3(dim_state, 0, 0)] = (Ldoub) 0; //varphi
		filter.A[index_3(dim_state, 0, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 0, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 0, 3)] = (Ldoub) 0; //Vox
		filter.A[index_3(dim_state, 0, 4)] = (Ldoub) 1./(R + Coordinates[2]); //Voy
		filter.A[index_3(dim_state, 0, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 0, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 0, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 0, 8)] = (Ldoub) 0;//U*cos(Coordinates[0]) + V[0]/(R*pow(cos(Coordinates[0]), 2)); //Phi_oz
		filter.A[index_3(dim_state, 0, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 0, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 0, 11)] = (Ldoub) 0; //d_omega_z
		//Delta dot lambda (долгота)
		filter.A[index_3(dim_state, 1, 0)] = (Ldoub) tan(Coordinates[0])*V[0]/((R + Coordinates[2])*cos(Coordinates[0])); //varphi
		filter.A[index_3(dim_state, 1, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 1, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 1, 3)] = (Ldoub) 1./((R + Coordinates[2])*cos(Coordinates[0])); //Vox
		filter.A[index_3(dim_state, 1, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 1, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 1, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 1, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 1, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 1, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 1, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 1, 11)] = (Ldoub) 0; //d_omega_z
		//Delta dot h (долгота)
		filter.A[index_3(dim_state, 2, 0)] = (Ldoub) 0; //varphi
		filter.A[index_3(dim_state, 2, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 2, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 2, 3)] = (Ldoub) 0; //Vox
		filter.A[index_3(dim_state, 2, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 2, 5)] = (Ldoub) 1; //Voz
		filter.A[index_3(dim_state, 2, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 2, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 2, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 2, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 2, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 2, 11)] = (Ldoub) 0; //d_omega_z
		// Delta dot V_ox
		filter.A[index_3(dim_state, 3, 0)] = (Ldoub) (V[0]/((R + Coordinates[2])*pow(cos(Coordinates[0]),2)) + (Ldoub) 2*U*cos(Coordinates[0]) )*V[1] + 2 * U * sin(Coordinates[0]) * V[2]; //varphi
		filter.A[index_3(dim_state, 3, 1)] = (Ldoub) 0; //lambda 
		filter.A[index_3(dim_state, 3, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 3, 3)] = (Ldoub) V[1]/(R+Coordinates[2])*tan(Coordinates[0]) - (Ldoub) V[2]/(R+Coordinates[2]); //Vox
		filter.A[index_3(dim_state, 3, 4)] = (Ldoub) V[0]/(R+Coordinates[2])*tan(Coordinates[0]) + (Ldoub) 2*U*sin(Coordinates[0]); //Voy
		filter.A[index_3(dim_state, 3, 5)] = (Ldoub) -(V[0]/(R+Coordinates[2]) + 2 * U * cos(Coordinates[0]) ); // Voz
		filter.A[index_3(dim_state, 3, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 3, 7)] = (Ldoub) -Ao[2]; //Phi_oy
		filter.A[index_3(dim_state, 3, 8)] = (Ldoub) Ao[1]; //Phi_oz
		filter.A[index_3(dim_state, 3, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 3, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 3, 11)] = (Ldoub) 0; //d_omega_z
		//Delta dot V_oy
		filter.A[index_3(dim_state, 4, 0)] = (Ldoub) -(V[0]/((R + Coordinates[2])*pow(cos(Coordinates[0]),2)) + (Ldoub) 2*U*cos(Coordinates[0]) )*V[0]; //varphi
		filter.A[index_3(dim_state, 4, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 4, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 4, 3)] = (Ldoub) -2.*(V[0]/(R + Coordinates[2])*tan(Coordinates[0]) + (Ldoub) U*sin(Coordinates[0])); //Vox
		filter.A[index_3(dim_state, 4, 4)] = (Ldoub) 0-V[2]*1/(R + Coordinates[2]); //Voy
		filter.A[index_3(dim_state, 4, 5)] = (Ldoub) - V[1]/(R+Coordinates[2]); //Voz
		filter.A[index_3(dim_state, 4, 6)] = (Ldoub) Ao[2]; //Phi_ox
		filter.A[index_3(dim_state, 4, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 4, 8)] = (Ldoub) -Ao[0]; //Phi_oz
		filter.A[index_3(dim_state, 4, 9)] = (Ldoub) 0;//d_omega_x
		filter.A[index_3(dim_state, 4, 10)] = (Ldoub) 0;//d_omega_y
		filter.A[index_3(dim_state, 4, 11)] = (Ldoub) 0;//d_omega_z
		//Delta dot V_oz
		filter.A[index_3(dim_state, 5, 0)] = (Ldoub) -2. * U * sin(Coordinates[0]) * V[0]; //varphi
		filter.A[index_3(dim_state, 5, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 5, 2)] = (Ldoub) 2*nu2; //height
		filter.A[index_3(dim_state, 5, 3)] = (Ldoub) 2.*(V[0]/(R + Coordinates[2]) + (Ldoub) U*cos(Coordinates[0])); //Vox
		filter.A[index_3(dim_state, 5, 4)] = (Ldoub) 2 * V[1]/(R + Coordinates[2]); //Voy
		filter.A[index_3(dim_state, 5, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 5, 6)] = (Ldoub) -Ao[1]; //Phi_ox
		filter.A[index_3(dim_state, 5, 7)] = (Ldoub) Ao[0]; //Phi_oy
		filter.A[index_3(dim_state, 5, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 5, 9)] = (Ldoub) 0;//d_omega_x
		filter.A[index_3(dim_state, 5, 10)] = (Ldoub) 0;//d_omega_y
		filter.A[index_3(dim_state, 5, 11)] = (Ldoub) 0;//d_omega_z
		//Phi_ox
		filter.A[index_3(dim_state, 6, 0)] = (Ldoub) 0; //varphi
		filter.A[index_3(dim_state, 6, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 6, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 6, 3)] = (Ldoub) 0; //Vox
		filter.A[index_3(dim_state, 6, 4)] = (Ldoub) -1./(R + Coordinates[2]); //Voy
		filter.A[index_3(dim_state, 6, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 6, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 6, 7)] = (Ldoub) omo[2]; //Phi_oy
		filter.A[index_3(dim_state, 6, 8)] = (Ldoub) -omo[1]; //Phi_oz
		filter.A[index_3(dim_state, 6, 9)] = (Ldoub) (-Cbn[index_3(3, 0, 0)]);//d_omega_x
		filter.A[index_3(dim_state, 6, 10)] = (Ldoub) (-Cbn[index_3(3, 0, 1)]);//d_omega_y
		filter.A[index_3(dim_state, 6, 11)] = (Ldoub) (-Cbn[index_3(3, 0, 2)]);//d_omega_z
		//Phi_oy
		filter.A[index_3(dim_state, 7, 0)] = (Ldoub) -U * sin(Coordinates[0]);//; //varphi
		filter.A[index_3(dim_state, 7, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 7, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 7, 3)] = (Ldoub) 1./(R+Coordinates[2]); //Vox
		filter.A[index_3(dim_state, 7, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 7, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 7, 6)] = (Ldoub) - omo[2]; //Phi_ox
		filter.A[index_3(dim_state, 7, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 7, 8)] = (Ldoub) omo[0]; //Phi_oz
		filter.A[index_3(dim_state, 7, 9)] = (Ldoub) (-Cbn[index_3(3, 1, 0)]); //d_omega_x
		filter.A[index_3(dim_state, 7, 10)] = (Ldoub) (-Cbn[index_3(3, 1, 1)]); //d_omega_y
		filter.A[index_3(dim_state, 7, 11)] = (Ldoub) (-Cbn[index_3(3, 1, 2)]); //d_omega_z
		#if 1
		//Phi_oz
		filter.A[index_3(dim_state, 8, 0)] = (Ldoub) U * cos(Coordinates[0]) + V[0]/(R*pow(cos(Coordinates[0]),2));//; //varphi
		filter.A[index_3(dim_state, 8, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 8, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 8, 3)] = (Ldoub) tan(Coordinates[0])/(R+Coordinates[2]); //Vox
		filter.A[index_3(dim_state, 8, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 8, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 8, 6)] = (Ldoub) omo[1]; //Phi_ox
		filter.A[index_3(dim_state, 8, 7)] = (Ldoub) -omo[0]; //Phi_oy
		filter.A[index_3(dim_state, 8, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 8, 9)] = (Ldoub) (-Cbn[index_3(3, 2, 0)]); //d_omega_x
		filter.A[index_3(dim_state, 8, 10)] = (Ldoub) (-Cbn[index_3(3, 2, 1)]); //d_omega_y
		filter.A[index_3(dim_state, 8, 11)] = (Ldoub) (-Cbn[index_3(3, 2, 2)]); //d_omega_z
		#endif
		//Delta omega_x
		filter.A[index_3(dim_state, 9, 0)] = (Ldoub) 0; //varphi
		filter.A[index_3(dim_state, 9, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 9, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 9, 3)] = (Ldoub) 0; //Vox
		filter.A[index_3(dim_state, 9, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 9, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 9, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 9, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 9, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 9, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 9, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 9, 11)] = (Ldoub) 0; //d_omega_z
		//Delta omega_y
		filter.A[index_3(dim_state, 10, 0)] = (Ldoub) 0; //varphi
		filter.A[index_3(dim_state, 10, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 10, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 10, 3)] = (Ldoub) 0; //Vox
		filter.A[index_3(dim_state, 10, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 10, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 10, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 10, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 10, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 10, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 10, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 10, 11)] = (Ldoub) 0; //d_omega_z
		//Delta omega_z
		filter.A[index_3(dim_state, 11, 0)] = (Ldoub) 0; //varphi
		filter.A[index_3(dim_state, 11, 1)] = (Ldoub) 0; //lambda
		filter.A[index_3(dim_state, 11, 2)] = (Ldoub) 0; //height
		filter.A[index_3(dim_state, 11, 3)] = (Ldoub) 0; //Vox
		filter.A[index_3(dim_state, 11, 4)] = (Ldoub) 0; //Voy
		filter.A[index_3(dim_state, 11, 5)] = (Ldoub) 0; //Voz
		filter.A[index_3(dim_state, 11, 6)] = (Ldoub) 0; //Phi_ox
		filter.A[index_3(dim_state, 11, 7)] = (Ldoub) 0; //Phi_oy
		filter.A[index_3(dim_state, 11, 8)] = (Ldoub) 0; //Phi_oz
		filter.A[index_3(dim_state, 11, 9)] = (Ldoub) 0; //d_omega_x
		filter.A[index_3(dim_state, 11, 10)] = (Ldoub) 0; //d_omega_y
		filter.A[index_3(dim_state, 11, 11)] = (Ldoub) 0; //d_omega_z
	#else //ветрок состояния 3
		filter.A[index_3(dim_state, 0, 0)] = 0; //Vn
		filter.A[index_3(dim_state, 0, 1)] = Ao[2]; //Phi_e
		filter.A[index_3(dim_state, 0, 2)] = 0; //d_omega_x
		//
		filter.A[index_3(dim_state, 1, 0)] = -1/(R+Coordinates[2]); //Vn
		filter.A[index_3(dim_state, 1, 1)] = 0; //Phi_e
		filter.A[index_3(dim_state, 1, 2)] = (Cbn[index_3(3, 0, 0)]); //d_omega_x
		//
		filter.A[index_3(dim_state, 2, 0)] = 0; //Vn
		filter.A[index_3(dim_state, 2, 1)] = 0; //Phi_e
		filter.A[index_3(dim_state, 2, 2)] = 0; //d_omega_x
		////////
	#endif
#if 1
	#if 0
		//Phi_x
		filter.G[index_3(filter.getDimQ(), (4), (0))] = (Ldoub) (-Cbn[index_3(3, 0, 0)]) * pow(h, 0);//d_omega_x
		filter.G[index_3(filter.getDimQ(), (4), (1))] = (Ldoub) (-Cbn[index_3(3, 0, 1)]) * pow(h, 0);//d_omega_y
		//Phi_y
		filter.G[index_3(filter.getDimQ(), (5), (0))] = (Ldoub) (-Cbn[index_3(3, 1, 0)]) * pow(h, 0);//d_omega_x
		filter.G[index_3(filter.getDimQ(), (5), (1))] = (Ldoub) (-Cbn[index_3(3, 1, 1)]) * pow(h, 0);//d_omega_y
		//Phi_z
		filter.G[index_3(filter.getDimQ(), (6), (0))] = (Ldoub) (-Cbn[index_3(3, 2, 0)]) * pow(h, 0);//d_omega_x
		filter.G[index_3(filter.getDimQ(), (6), (1))] = (Ldoub) (-Cbn[index_3(3, 2, 1)]) * pow(h, 0);//d_omega_y
	#endif
		// как было до этого
		int pow_h = 1;
		//Phi_x
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 3), (filter.getDimQ() - 3))] = (Ldoub) (-Cbn[index_3(3, 0, 0)]) * pow(h, pow_h);//d_omega_x
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 3), (filter.getDimQ() - 2))] = (Ldoub) (-Cbn[index_3(3, 0, 1)]) * pow(h, pow_h);//d_omega_y
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 3), (filter.getDimQ() - 1))] = (Ldoub) (-Cbn[index_3(3, 0, 2)]) * pow(h, pow_h);//d_omega_z
		//Phi_y
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 2), (filter.getDimQ() - 3))] = (Ldoub) (-Cbn[index_3(3, 1, 0)]) * pow(h, pow_h);//d_omega_x
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 2), (filter.getDimQ() - 2))] = (Ldoub) (-Cbn[index_3(3, 1, 1)]) * pow(h, pow_h);//d_omega_y
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 2), (filter.getDimQ() - 1))] = (Ldoub) (-Cbn[index_3(3, 1, 2)]) * pow(h, pow_h);//d_omega_z
		//Phi_z
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 1), (filter.getDimQ() - 3))] = (Ldoub) (-Cbn[index_3(3, 2, 0)]) * pow(h, pow_h);//d_omega_x
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 1), (filter.getDimQ() - 2))] = (Ldoub) (-Cbn[index_3(3, 2, 1)]) * pow(h, pow_h);//d_omega_y
		filter.G[index_3(filter.getDimQ(), (filter.getDimX() - 1), (filter.getDimQ() - 1))] = (Ldoub) (-Cbn[index_3(3, 2, 2)]) * pow(h, pow_h);//d_omega_z
	#endif
#if 0
		printf("G\n");
		filter.Print2dMatr(filter.G, filter.getDimX(), filter.getDimQ());

		#endif 
		if (!filter.init) //Если ранее не было инициализации, то инициализируем
		{
			filter.Init(x0, q, r, H);
		}
		else//в противном случае оцениваем
		{
			/*Тут было составление матрицы Phi, но я это перенес в AdaptiveKalman::Predict*/
		#if 1
			for (int iii=0; iii<2; ++iii)
				ErrVins[iii] = Coordinates[iii] - CooGps[iii]; //разница ошибок координат ИНС и СНС
		#endif
			for (int iii=2; iii<filter.getDimZ(); ++iii)
				ErrVins[iii] = V[iii - 2] - Vgps[iii - 2]; //разница ошибок скоростей ИНС и СНС
			filter.Predict();
			filter.Update(ErrVins);
		#if 0
			printf("v:\n");
			filter.Print2dMatr(filter.v, filter.getDimZ(), 1);
			printf("x:\n");
			filter.Print2dMatr(filter.x, filter.getDimX(), 1);
			printf("C:\n");
			filter.Print2dMatr(filter.C, filter.getDimZ(), filter.getDimZ());
			printf("K:\n");
			filter.Print2dMatr(filter.K, filter.getDimX(), filter.getDimZ());
			// printf("omega_x = %.8f \t omega_y = %.8f\n", filter.x[5], filter.x[6]);
		#endif

		#if 0
			for (int iii=0; iii<2; ++iii)
				V[iii] -= filter.x[iii]; // коррекция скоростей ИНС
		#endif
		}
	#endif
		// Запись в файл навигационного решения (скорости, координаты, углы, ошибки по скоростям, ошибки по координатам)

		static FILE* navig_res;
		if(!navig_res)
		{
			navig_res=fopen("./data/Nav_res.csv","wt");
			// Шапка
			fprintf(navig_res, "Woe;");
			fprintf(navig_res, "Won;");
			fprintf(navig_res, "Ve;");
			fprintf(navig_res, "Vn;");
			fprintf(navig_res, "Vup;");
			fprintf(navig_res, "errVe;");
			fprintf(navig_res, "errVn;");
			fprintf(navig_res, "errVup;");
			fprintf(navig_res, "Phi;");
			fprintf(navig_res, "Lambda;");
			fprintf(navig_res, "Height;");
			fprintf(navig_res, "Heading;");
			fprintf(navig_res, "Roll;");
			fprintf(navig_res, "Pitch;");
			fprintf(navig_res, "d_VE;");
			fprintf(navig_res, "d_VN;");
			fprintf(navig_res, "d_VUp;");
			fprintf(navig_res, "d_E;");
			fprintf(navig_res, "d_N;");
			fprintf(navig_res, "d_h;");
			fprintf(navig_res, "\n");
		}

		if(navig_res)
		{
			// Навигационные параметры
			// Малые приращения скоростей
			for(int i=0; i<2; ++i)
				fprintf(navig_res, "%.10Lf;", Wo[i]);
			// Скорости
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10Lf;", V[i]);
			//Ошибки по скорости
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10Lf;", Err_V[i]);
			// Координаты
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10Lf;", Coordinates[i]);
			// Углы оориентации
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10Lf;", Orientation[i]);
			//Ошибки по cкоростям в м/с
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10Lf;", Err_V[i]);
			//Ошибки по координатам в м
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10Lf;", CoordError[i]);

			fprintf(navig_res, "\n");
		}
#endif
		fflush(navig_res);

#if 0 //Запись в файл данных для оценивания дрейфов программой для дипломной работы (на Python)
		static FILE* for_Kalman;
		if(!for_Kalman)
		{
			for_Kalman=fopen("./data/For_Kalman.csv","wt"); //Timestamp,a_ll_x,a_ll_y,a_ll_z,omega_s_x,omega_s_y,omega_s_z,Ve_ins,Vn_ins,heading,roll,pitch,latitude,longitude,Ve_gps,Vn_gps
			// Шапка
			fprintf(for_Kalman, "Timestamp,");
			fprintf(for_Kalman, "a_ll_x,");
			fprintf(for_Kalman, "a_ll_y,");
			fprintf(for_Kalman, "a_ll_z,");
			fprintf(for_Kalman, "omega_s_x,");
			fprintf(for_Kalman, "omega_s_y,");
			fprintf(for_Kalman, "omega_s_z,");
			fprintf(for_Kalman, "Ve_ins,");
			fprintf(for_Kalman, "Vn_ins,");
			fprintf(for_Kalman, "heading,");
			fprintf(for_Kalman, "roll,");
			fprintf(for_Kalman, "pitch,");
			fprintf(for_Kalman, "latitude,");
			fprintf(for_Kalman, "longitude,");
			fprintf(for_Kalman, "Ve_gps,");
			fprintf(for_Kalman, "Vn_gps");
			fprintf(for_Kalman, "\n");
		}

		if(for_Kalman)
		{
			fprintf(for_Kalman, "%d,", cur_time);

			// Ускорения в опорной с.к
			for(int i=0; i<3; ++i)
				fprintf(for_Kalman, "%.10e,", Ao[i]);
			// Асболютная (?) угловая скорость опорной с.к
			for(int i=0; i<3; ++i)
				fprintf(for_Kalman, "%.10e,", omo[i]);
			// Скорости в горизонте
			for(int i=0; i<2; ++i)
				fprintf(for_Kalman, "%.10e,", V[i]);
			// Углы оориентации
			for(int i=0; i<3; ++i)
				fprintf(for_Kalman, "%.10e,", Orientation[i]);
			//Скоординаты
			for(int i=0; i<2; ++i)
				fprintf(for_Kalman, "%.10e,", Coordinates[i]);
			//Скорости от GPS
			for (int iii=0; iii<2; ++iii)
				fprintf(for_Kalman, "%.10e,", Vgps[iii]);

			fprintf(for_Kalman, "\n");
		}
		fflush(for_Kalman);
#endif

#if 1 
		//Запись в файл оцененного вектора состояния
		static FILE* estimations;
		if(!estimations)
		{
			estimations=fopen("./data/Kalman_est.csv","wt"); //Timestamp,a_ll_x,a_ll_y,a_ll_z,omega_s_x,omega_s_y,omega_s_z,Ve_ins,Vn_ins,heading,roll,pitch,latitude,lingitude,Ve_gps,Vn_gps
			// Шапка
		#if 1 //Для вектора состояния размерности 10
			fprintf(estimations, "varphi,");		// 0
			fprintf(estimations, "lambda,");		// 1
			fprintf(estimations, "height,");		// 2
			fprintf(estimations, "Ve,");			// 3
			fprintf(estimations, "Vn,");			// 4
			fprintf(estimations, "Vup,");			// 5
			fprintf(estimations, "Phi_e,");			// 6
			fprintf(estimations, "Phi_n,");			// 7
			fprintf(estimations, "Phi_z,");			// 8
			fprintf(estimations, "d_omega_x,");		// 9
			fprintf(estimations, "d_omega_y,");		// 10
			fprintf(estimations, "d_omega_z,");		// 11
			fprintf(estimations, "d_Psi,");			// 12
			fprintf(estimations, "d_Roll,");		// 13
			fprintf(estimations, "d_Pitch,");		// 14
			fprintf(estimations, "\n");
		#else
			fprintf(estimations, "Ve,");
			fprintf(estimations, "Phi_n,");
			fprintf(estimations, "d_omega_y,");
			fprintf(estimations, "d_Roll,");
			fprintf(estimations, "d_Pitch,");
			fprintf(estimations, "\n");
		#endif
		}

		if(estimations)
		{
			fprintf(estimations, "%d,", cur_time);

			// весь вектор состояния
			for(int i=0; i<dim_state; ++i)
				fprintf(estimations, "%.10e,", filter.x[i]);
			double d_psi (filter.x[8] - (filter.x[4+2] * sin(Orientation[0]) + filter.x[5+2] * cos(Orientation[0]))*tan(Orientation[2] )); //это угол курса, а нек ошибки курса
			double d_roll = -(filter.x[5+2] * cos(Orientation[0]) + filter.x[4+2] * sin(Orientation[0])) * (1./cos(Orientation[2]));//Ошибка крена по (ошибкам?) ориентации Fx, Fy
			double d_pitch = -(filter.x[4+2] * cos(Orientation[0]) - filter.x[5+2] * sin(Orientation[0]));//Ошибка крена по (ошибкам?) ориентации Fx, Fy
			fprintf(estimations, "%.10e,", d_psi);
			fprintf(estimations, "%.10e,", d_roll);
			fprintf(estimations, "%.10e,", d_pitch);
			fprintf(estimations, "\n");
		}
		
		fflush(estimations);
#endif
	} //чтение из файла while( !in.eof())
	printf("Алгоритм ИНС закончил моделирование\n");
	return 0;
}
