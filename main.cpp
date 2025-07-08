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
#include "AdaptiveKalman.h" //Для адаптивного фильтра Калмана
#include "UsualKalman.h" //Для обычного фильтра Калмана
#define _USE_MATH_DEFINES
#include <stdlib.h>

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

int Readfile(FILE* is, bool AllowBiasAcc, bool AllowBiasGyr, bool AllowRandAcc, bool AllowRandGyr, bool AllowRandVgps, Ldoub* Ab, Ldoub* Omb, Ldoub* Vgps, Ldoub* CoordGps)
{
	Ldoub BiasAb[3] = {0}; // Постоянные погрешности акселерометров
	Ldoub BiasOmb[3] = {0}; // Постоянные погрешности гироскопов
	Ldoub RandAb[3] = {0}; // Случайные погрешности акселерометров
	Ldoub RandOmb[3] = {0};// Случайные погрешности гироскопов
	Ldoub RandomVgps[2] = {0}; // Случайные погрешности скоростей GPS
	Ldoub RandomCooGPS[2] = {0}; // Случайные погрешности координат GPS

	int res = fscanf(is, "%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n", &Ab[0], &Ab[1], &Ab[2], &Omb[0], &Omb[1], &Omb[2], &BiasAb[0], &BiasAb[1], &BiasAb[2], &BiasOmb[0], &BiasOmb[1], &BiasOmb[2], &RandAb[0], &RandAb[1], &RandAb[2], &RandOmb[0], &RandOmb[1], 	&RandOmb[2], &Vgps[0], &Vgps[1], &RandomVgps[0], &RandomVgps[1], &CoordGps[0], &CoordGps[1], &RandomCooGPS[0], &RandomCooGPS[1]);
	if (res == 26)
	{
		//добавление дрейфов к показаниям инерциальных датчиков
		for (int ii=0; ii<3; ++ii)
		{
			Ab[ii] += AllowBiasAcc*BiasAb[ii] + AllowRandAcc*RandAb[ii];
		}
		Omb[0] += AllowBiasGyr*BiasOmb[0] + AllowRandGyr*RandOmb[0];
		Omb[1] += AllowBiasGyr*BiasOmb[1] + AllowRandGyr*RandOmb[1];
		Omb[2] += 0*AllowBiasGyr*BiasOmb[2] + AllowRandGyr*RandOmb[2];
		//добавление шума к показаниям СНС
		for(int iii=0; iii<2; ++iii)
		{
			Vgps[iii] += (Ldoub) AllowRandVgps*RandomVgps[iii];
			CoordGps[iii] += (Ldoub) AllowRandVgps*RandomCooGPS[iii];
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
	const Ldoub a = 6378245;
	const Ldoub b = 6356856;
	const Ldoub e = sqrt(1 - b*b/a/a);
	const Ldoub R = 6400e3;
	const Ldoub pi = 3.141592653589793;
	const Ldoub U = 7.27220521664304e-05;
	const Ldoub rad2deg = 180./M_PI; // из градусов в час в радианы в секунду
	const Ldoub deg2rad = 1./rad2deg;
	int freq1 = 400; // частота свехбыстрого цикла
	Ldoub h1 (1./freq1); // период дискретизации свехбыстрого цикла
	int freq = freq1/4;//100; // частота измерений с инерциальных датчиков
	Ldoub h (1./freq); //период дискретизации
	int t_nav = 180*60; // время работы нав алгоритма в секундах
	int t_alignment = (int) 5*60*freq1; // время выставки в тактах
	Ldoub Rlambda;
	Ldoub Rphi;
	Ldoub phi0 = (Ldoub) 55*deg2rad;// и для моделирования
	Ldoub lambda0 = (Ldoub) 33*deg2rad;
	int cur_time = 0; // текущий такт!! измерения
	// Для моделирования показаний Ч.Э.
	Ldoub H0 = (Ldoub) (strtod(argv[1], NULL))*deg2rad;
	Ldoub P0 = (Ldoub) (strtod(argv[2], NULL))*deg2rad;
	Ldoub R0 = (Ldoub) (strtod(argv[3], NULL))*deg2rad;
	Ldoub Vabs = (Ldoub) strtod(argv[4], NULL);

	//printf("H0 = %.4f P0 = %.4f R0 = %.4f Vabs = %.4f\n", H0, P0, R0, Vabs);

	Ldoub Vgps[2] = {0};
	Ldoub CoordGps[2] = {0};
	Ldoub Cnb[9];
	MatrOB(H0, R0, P0, Cnb, 3); // матрица перехода из опорной в связанную
	//Необходимое для выставки
	Ldoub DeltaHeading = 0, DeltaRoll = 0, DeltaPitch = 0;// ошибки выставки по курсу, крену и тангажу соответственно
	Ldoub Heading = 0, Roll = 0, Pitch = 0;
	bool AlignmentContinue = true; // для начала выставки
	// Необходимые массивы для решение навигационной задачи
	Ldoub Ab[3]={0}; // Ускорения в связанных осях
	Ldoub Ao[3] = {0}; // Ускорения в опорныых осях
	Ldoub Wo[3] = {0}; //Приращение скоростей в опорных осях
	Ldoub Omb[3] = {0}; // Угловые скорости в связанных осях
	Ldoub Omo[3] = {0}; // Угловые скорости в опорных осях
	//Ldoub Oms[3] = {0}; // угловые скорости от линейного движения + Земля
	Ldoub MeanAb[3] = {0};
	Ldoub MeanOmb[3] = {0};
	Ldoub StdAb[3] = {0};
	Ldoub StdOmb[3] = {0};
	Ldoub Cbn[9] = {0};
	Ldoub Cib[9] = {0}; //матрица перехода из инерциальной системы в связанную. Начальное значение равно транспонированной матрицы на момент окончания выставки
	Ldoub Cin[9] = {1.,0,0,0,1.,0,0,0,1.}; //матрица перехода из инерцальной в опорную. Начальное знвчение -- единичная Cin(0)=E
	Ldoub V0[3] = {(Ldoub) Vabs*sin(H0), (Ldoub) Vabs*cos(H0), 0}; // линейные скорости E; N; Up
	Ldoub Err_V[3] = {0}; // ошибки по скоростям
	Ldoub ErrVins[3] = {0}; // разница ошибок между ИНС и GPS (для корректирующих поправок)
	// массивы для выходных значений
	Ldoub V[3] = {(Ldoub) V0[0], (Ldoub) V0[1], 0}; // линейные скорости E; N; Up (на всякий случай, пока резерв)
	Ldoub Coordinates[3] = {phi0, lambda0, 0}; // Географические кординаты: широта, долгота и высота
	Ldoub CoordError[2] = {0}; // Ошибки в м (dE, dN)
	Ldoub Orientation[3] = {0}; // Углы ориентации

	Rlambda = (Ldoub) R/sqrt(1.-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
	Rphi = (Ldoub) R*(1. - e*e)/(sqrt(1.-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1.-e*e*sin(Coordinates[0])*sin(Coordinates[0])));

	bool AllowBiasAcc = (bool) (strtod(argv[6], NULL));
	bool AllowBiasGyr = (bool) (strtod(argv[7], NULL));
	bool AllowRandAcc = (bool) (strtod(argv[8], NULL));
	bool AllowRandGyr = (bool) (strtod(argv[9], NULL));
	bool AllowRandVgps = (bool) (strtod(argv[10], NULL));
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
	FILE* file=fopen(argv[5], "rt");
	//fscanf(file, "%*s;");//чтение строки заголовка, она не нужна, выбрасываем
	int res = Readfile(file, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps, CoordGps);
#else
	std::ifstream in(argv[5], std::ios::binary);// !in.eof() //условие цикла while для бинарного файла
	#define FILE_STREAM
#endif
	Ldoub resultV[2] = {0};
	Ldoub Verr1[2] = {0}; // ошибки интегрирования ускорений
	Ldoub Verr2[2] = {0}; // ошибки накопления скоростей
	///*
	const int dim_state (6); //размер вектора состояния
	const int dim_sense (2); //размер вектора измерения
	//*/
	AdaptiveKalman filter(dim_state, dim_sense); //Создаю объект фильтра Калмана с матрицей размера 6*6 и измерениями 2*1 (вертикальную скорость не учитываю)
	Ldoub x0[dim_state] = {0}; //Начальные оценочные значения вектора состояния
	Ldoub H[dim_sense*dim_state] = {0}; //матрица наблюдения
	#if 1 // вектор состояния 6
	H[index_3(dim_state, 0, 0)] = 1;
	H[index_3(dim_state, 1, 1)] = 1;
	//Матрца ковариации входных шумов (модели)
	Ldoub q[dim_state*dim_state] = {0};
	q[index_3(dim_state, 4, 4)] = 1e-16;
	q[index_3(dim_state, 5, 5)] = 1e-16;
	#else // вектор состояния 3
	H[index_3(dim_state, 0, 0)] = 1;
	//Матрца ковариации входных шумов (модели)
	Ldoub q[dim_state*dim_state] = {0};
	q[index_3(dim_state, 0, 0)] = 1e-19;
	#endif
	// Ldoub r[dim_sense*dim_sense] = {0.05, 0, 0, 0.05};

#if 0
	//Для лучшей обусловленности матрицы HPH_t увеличиваю начальные значения априорной ошибки
	for (int iii=0; iii<filter.getDimX()*filter.getDimX(); ++iii)
	{
		filter.Papr[index_3(filter.getDimX(), iii, iii)] = 1;
	}
#endif

	while(res == 26 ) // Пока возможно чтение из файла. Для бинарного файла -- !in.eof() //условие цикла while для бинарного файла //  Readfile(file, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps) == 22 --- условие цикла While для текстового файла
	{ // цилк while для частоты 100 Гц, поэтому нельзя сделать общее чтение файла для выставки и навигации

		// этап выставки
		if ((cur_time <= t_alignment) && AlignmentContinue )
		{
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
			#ifdef FILE_STREAM //если определен файловый поток, то чтение из бинарника, читаем тут
			ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps); //чтение из бинарного файла данных используемых для выставки
			#else
			res = Readfile(file, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps, CoordGps);
			#endif
			alignment(Ab, Omb, MeanAb, MeanOmb, StdAb, StdOmb, cur_time, g, U, phi0, Cbn);
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
				printf("Heading (degrees) = %.30f Error = %f\n", Heading*rad2deg, DeltaHeading*rad2deg);
				printf("roll (degrees) = %.30f Error = %f\n", Roll*rad2deg, DeltaRoll*rad2deg);
				printf("pitch (degrees) = %.30f Error = %f\n", Pitch*rad2deg, DeltaPitch*rad2deg);
				AlignmentContinue = false;
				Transpose2M(Cbn, Cib, 3); // начальное значение Cib; Cib(0)
				//Transpose(Cib,3);
				cur_time = 0; // для 100 Гц

				Ldoub gravity[3] = {0, 0, g}; //ускорения силы тяжести в проекциях на оси опорной системы координат (географической)
				Ldoub ErrAcc[3] = {0}; //Ошибки акселерометров (нужны для невыставки)
				MulMatrD(Cib, gravity, ErrAcc,3,3,1); //здесь ErrAcc как временная матрица, а Cib=(Cbn)^t в начальный момент времени
				for (int iii=0; iii<3; ++iii)
					ErrAcc[iii] -= MeanAb[iii]; //теперь ErrAcc есть ошибки акселерометров
#if 0
				//начальные значения ошибок ориентации (для вектора состояния)
				x0[2] = -1e-4/g;
				x0[3] = 1e-4/g;
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
		for (int in_iter=0; in_iter < 4; ++in_iter) // 4 такта, нумерация с нуля, поэтому равенство нестрогое
		{
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
			#ifdef FILE_STREAM //если определен файловый поток, то чтение из бинарника, читаем тут
			ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps); //чтение из бинарного файла данных используемых для выставки
			#else //в противном случае, читаем, что есть
			res = Readfile(file, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps, CoordGps);
			#endif
		#if 1
			//Накапливаем данные 4 тактов и заодно осредним псевдокоординаты
			for (int iii=0; iii<3; ++iii)
			{
				alpha[index_3(4, iii, in_iter)] = (Ldoub) Omb[iii]*h1;
				MeanAlpha[iii] = (Ldoub) MeanAlpha[iii] + (Omb[iii] - MeanAlpha[iii]) / (in_iter + 1);
				w[index_3(4, iii, in_iter)] = (Ldoub) Ab[iii]*h1;
				MeanW[iii] = (Ldoub) MeanW[iii] + (w[iii] - MeanW[iii]) / (in_iter + 1);
			}
			//Определение приращения скорости для каждого значения сверхбыстрого цикла
			//кососиметрическая матрица псевдокоординат
			Ldoub EigAl[9] = {0, -Omb[2]*h1, Omb[1]*h1, Omb[2]*h1, 0, -Omb[0]*h1, -Omb[1]*h1, Omb[0]*h1, 0};
			//вычисляем k1
			Ldoub k1[3] = {0};
			Ldoub al_w[3] = {0}; //вспомогательная матрица

			MulMatrD(EigAl, Wp, al_w, 3,3,1);

			for (int nnn=0; nnn<3; ++nnn)
			{
				k1[nnn] = w[index_3(4, nnn, in_iter)] - al_w[nnn];
			}
			//вычисляем k2
			Ldoub k2[3]={0.};
			Ldoub temp_w_k[3] = {0.};

			for (int i=0; i<3; i++)
				temp_w_k[i] = (Ldoub) Wp[i] + (Ldoub) h1/2.*k1[i];

			MulMatrD(EigAl, temp_w_k, al_w, 3,3,1);
			for (int nnn=0; nnn<3; ++nnn)
			{
				k2[nnn] = (Ldoub) w[index_3(4, nnn, in_iter)] - (Ldoub) al_w[nnn];
			}
			//вычисляем k3
			Ldoub k3[3]={0.};

			for (int i=0; i<3; i++) temp_w_k[i] = (Ldoub) Wp[i] + (Ldoub) h1/2.*k2[i];

			MulMatrD(EigAl, temp_w_k, al_w, 3,3,1);
			for (int nnn=0; nnn<3; ++nnn)
			{
				k3[nnn] = (Ldoub) w[index_3(4, nnn, in_iter)] - (Ldoub) al_w[nnn];
			}
			//вычисляем k4
			Ldoub k4[3]={0.};
			for (int i=0; i<3; i++)
				temp_w_k[i] = (Ldoub) Wp[i] + (Ldoub) h1*k3[i];

			MulMatrD(EigAl, temp_w_k, al_w, 3,3,1);

			for (int nnn=0; nnn<3; ++nnn)
			{
				k4[nnn] = (Ldoub) w[index_3(4, nnn, in_iter)] - (Ldoub) al_w[nnn];
			}
			//считаем приращение скорости
			for (int ii=0; ii<3; ++ii)
			{
				Wp[ii] = (Ldoub) Wp[ii] +  (Ldoub) 1./6*(k1[ii] + 2.*k2[ii] + 2.*k3[ii] + k4[ii]);
			}
		#endif
		}
		/*Далее идет 100 Гц такт*/
		//Решение задачи ориентации
		Ldoub omo[3] = {0}; //переносные угловые скорости, вчисляются в решении задачи ориентации (SolveOrient)
		SolveOrient(alpha, Cib, Cin, Cbn, Orientation, Coordinates, Omo, omo, V, phi0, Rphi, Rlambda, freq, h, U, cur_time, derectNorm, kcor2, ErrVins, allowCorr); //из одноименного заголовочного файла

		/*Решение задачи навигации*/
		SolveNav(Wp, Ab, Cbn, Wo, Ao, V,  Coordinates, CoordError, Err_V, omo, h, Rphi, Rlambda, U, R, e, H0, V0, kcor1, ErrVins, allowCorr); //Err_V уже в этой функции вычисляется, поэтому я могу это значение использовать для коррекции
		// инкремент тактов
		++cur_time;

		/*
		По идее, куда-то сюда можно засунуть оценивание по Калману скоростей дрейфов гироскопов
		*/
	#if 1
	#if 1 //Модель для вектрора состояния 6
		//Каждый такт пересчитываем матрицу A у фильтра Калмана
		// Delta dot V_ox
		filter.A[index_3(dim_state, 0, 0)] = V[1]/(R+Coordinates[2])*tan(Coordinates[0]); //Vox
		filter.A[index_3(dim_state, 0, 1)] = V[0]/(R+Coordinates[2])*tan(Coordinates[0]) + 2*U*sin(Coordinates[0]); //Voy
		filter.A[index_3(dim_state, 0, 2)] = 0; //Phi_ox
		filter.A[index_3(dim_state, 0, 3)] =  -Ao[2]; //Phi_oy
		filter.A[index_3(dim_state, 0, 4)] = 0; //d_omega_x
		filter.A[index_3(dim_state, 0, 5)] = 0; //d_omega_y
		//Delta dot V_oy
		filter.A[index_3(dim_state, 1, 0)] = -2*(V[0]/(R+Coordinates[2])*tan(Coordinates[0]) + U*sin(Coordinates[0])); //Vox
		filter.A[index_3(dim_state, 1, 1)] = 0; //Voy
		filter.A[index_3(dim_state, 1, 2)] = Ao[2]; //Phi_ox
		filter.A[index_3(dim_state, 1, 3)] = 0; //Phi_oy
		filter.A[index_3(dim_state, 1, 4)] = 0;//d_omega_x
		filter.A[index_3(dim_state, 1, 5)] = 0;//d_omega_y
		//Phi_ox
		filter.A[index_3(dim_state, 2, 0)] = 0; //Vox
		filter.A[index_3(dim_state, 2, 1)] = -1/(R+Coordinates[2]); //Voy
		filter.A[index_3(dim_state, 2, 2)] = 0; //Phi_ox
		filter.A[index_3(dim_state, 2, 3)] = omo[2]; //Phi_oy
		filter.A[index_3(dim_state, 2, 4)] = -Cbn[index_3(3, 0, 0)];//c11*omega_bx 
		filter.A[index_3(dim_state, 2, 5)] = -Cbn[index_3(3, 0, 1)];// c12 * omega_by
		//Phi_oy
		filter.A[index_3(dim_state, 3, 0)] = 1/(R+Coordinates[2]); //Vox
		filter.A[index_3(dim_state, 3, 1)] = 0; //Voy
		filter.A[index_3(dim_state, 3, 2)] = - omo[2]; //Phi_ox
		filter.A[index_3(dim_state, 3, 3)] = 0; //Phi_oy
		filter.A[index_3(dim_state, 3, 4)] = -Cbn[index_3(3, 1, 0)]; //d_omega_x
		filter.A[index_3(dim_state, 3, 5)] = -Cbn[index_3(3, 1, 1)]; //d_omega_y
		//Delta omega_x
		filter.A[index_3(dim_state, 4, 0)] = 0; //Vox
		filter.A[index_3(dim_state, 4, 1)] = 0; //Voy
		filter.A[index_3(dim_state, 4, 2)] = 0; //Phi_ox
		filter.A[index_3(dim_state, 4, 3)] = 0; //Phi_oy
		filter.A[index_3(dim_state, 4, 4)] = 0; //d_omega_x
		filter.A[index_3(dim_state, 4, 5)] = 0; //d_omega_y
		//Delta omega_y
		filter.A[index_3(dim_state, 5, 0)] = 0; //Vox
		filter.A[index_3(dim_state, 5, 1)] = 0; //Voy
		filter.A[index_3(dim_state, 5, 2)] = 0; //Phi_ox
		filter.A[index_3(dim_state, 5, 3)] = 0; //Phi_oy
		filter.A[index_3(dim_state, 5, 4)] = 0; //d_omega_x
		filter.A[index_3(dim_state, 5, 5)] = 0; //d_omega_y
	#else //ветрок состояния 3
		filter.A[index_3(dim_state, 0, 0)] = 0; //Ve
		filter.A[index_3(dim_state, 0, 1)] = -Ao[2]; //Phi_N
		filter.A[index_3(dim_state, 0, 2)] = 0; //d_omega_N
		//
		filter.A[index_3(dim_state, 1, 0)] = 1/(R+Coordinates[2]); //Ve
		filter.A[index_3(dim_state, 1, 1)] = 0; //Phi_N
		filter.A[index_3(dim_state, 1, 2)] = 1; //d_omega_N
		//
		filter.A[index_3(dim_state, 2, 0)] = 0; //Ve
		filter.A[index_3(dim_state, 2, 1)] = 0; //Phi_N
		filter.A[index_3(dim_state, 2, 2)] = 0; //d_omega_N
	#endif
		if (!filter.init) //Если ранее не было инициализации, то инициализируем
		{
			filter.Init(x0, q, /*r,*/ H);
		}
		else//в противном случае оцениваем
		{
			for (int iii =0; iii < filter.getDimX()*filter.getDimX(); ++iii) //вычисляю матрицу перехода Phi
				filter.Phi[iii] = filter.I[iii] + filter.A[iii] * h; //не забываем умножить на такт интегрирования

			for (int iii=0; iii<filter.getDimZ(); ++iii)
				ErrVins[iii] = V[iii] - Vgps[iii]; //разница ошибок ИНС и СНС
			filter.Predict();
			filter.Update(ErrVins);
			// printf("omega_x = %.8f \t omega_y = %.8f\n", filter.x[5], filter.x[6]);

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
			fprintf(navig_res, "Ve;");
			fprintf(navig_res, "Vn;");
			fprintf(navig_res, "Phi;");
			fprintf(navig_res, "Lambda;");
			fprintf(navig_res, "Height;");
			fprintf(navig_res, "Heading;");
			fprintf(navig_res, "Roll;");
			fprintf(navig_res, "Pitch;");
			fprintf(navig_res, "d_VE;");
			fprintf(navig_res, "d_VN;");
			fprintf(navig_res, "d_E;");
			fprintf(navig_res, "d_N;");
			fprintf(navig_res, "\n");
		}

		if(navig_res)
		{
			// Навигационные параметры
			// Скорости
			for(int i=0; i<2; ++i)
				fprintf(navig_res, "%.10e;", V[i]);
			// Координаты
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10e;", Coordinates[i]);
			// Углы оориентации
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%.10e;", Orientation[i]);
			//Ошибки по cкоростям в м/с
			for(int i=0; i<2; ++i)
				fprintf(navig_res, "%.10e;", Err_V[i]);
			//Ошибки по координатам в м
			for(int i=0; i<2; ++i)
				fprintf(navig_res, "%.10e;", CoordError[i]);

			fprintf(navig_res, "\n");
		}
#endif
		fflush(navig_res);

#if 1 //Запись в файл данных для оценивания дрейфов программой для дипломной работы (на Python)
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
		#if 1 //Для вектора состояния размерности 6
			fprintf(estimations, "Ve,");
			fprintf(estimations, "Vn,");
			fprintf(estimations, "Phi_e,");
			fprintf(estimations, "Phi_n,");
			fprintf(estimations, "d_omega_x,");
			fprintf(estimations, "d_omega_y,");
			fprintf(estimations, "d_Roll,");
			fprintf(estimations, "d_Pitch,");
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
			double d_roll = -(filter.x[3] * cos(/*Orientation[0]*/H0) + filter.x[2] * sin(/*Orientation[0]*/H0)) * (1./cos(/*Orientation[2]*/P0));//Ошибка крена по (ошибкам?) ориентации Fx, Fy
			double d_pitch = -(filter.x[2] * cos(/*Orientation[0]*/H0) - filter.x[3] * sin(/*Orientation[0]*/H0));//Ошибка крена по (ошибкам?) ориентации Fx, Fy
			fprintf(estimations, "%.10e,", d_roll);
			fprintf(estimations, "%.10e,", d_pitch);
			fprintf(estimations, "\n");
		}
		fflush(estimations);
#endif
	} //чтение из файла while( !in.eof())
	return 0;
}
