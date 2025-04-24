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
#define _USE_MATH_DEFINES
#include <stdlib.h>

//#define dev

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

int Readfile(FILE* is, bool AllowBiasAcc, bool AllowBiasGyr, bool AllowRandAcc, bool AllowRandGyr, bool AllowRandVgps, Ldoub* Ab, Ldoub* Omb, Ldoub* Vgps)
{
	Ldoub BiasAb[3] = {0}; // Постоянные погрешности акселерометров
	Ldoub BiasOmb[3] = {0}; // Постоянные погрешности гироскопов
	Ldoub RandAb[3] = {0}; // Случайные погрешности акселерометров
	Ldoub RandOmb[3] = {0};// Случайные погрешности гироскопов
	Ldoub RandomVgps[2] = {0};

	int res = fscanf(is, "%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n", Ab[0], Ab[1], Ab[2], Omb[0], Omb[1], Omb[2], BiasAb[0], BiasAb[1], BiasAb[2], BiasOmb[0], BiasOmb[1], BiasOmb[2], RandAb[0], RandAb[1], RandAb[2], RandOmb[0], RandOmb[1], RandOmb[2], Vgps[0], Vgps[1]);
	if (res == 22)
	{
		//добавление дрейфов к показаниям инерциальных датчиков
		for (int ii=0; ii<3; ++ii)
		{
			Ab[ii] += AllowBiasAcc*BiasAb[ii] + AllowRandAcc*RandAb[ii];
			Omb[ii] += AllowBiasGyr*BiasOmb[ii] + AllowRandGyr*RandOmb[ii];
		}
		//добавление шума к показаниям СНС
		for(int iii=0; iii<2; ++iii)
		{
			Vgps[iii] += (Ldoub) AllowRandVgps*RandomVgps[iii];
		}

		return res;
	}
	else return -1; //означает что что-то не так
}
#if 0
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
	int freq = 100; // частота измерений с инерциальных датчиков
	Ldoub h = 0.01; //период дискретизации
	int freq1 = 400; // частота свехбыстрого цикла
	Ldoub h1 (1./400); // период дискретизации свехбыстрого цикла
	int t_nav = 180*60; // время работы нав алгоритма в секундах
	int t_alignment = (int) 5*60*freq1; // время выставки в тактах
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

	Ldoub Vgps[2] = {0};
	Ldoub Cnb[9];
	MatrOB(H0, R0, P0, Cnb, 3); // матрица перехода из опорной в связанную
	//Необходимое для выставки
	Ldoub lambda0 = (Ldoub) 33*deg2rad;
	Ldoub DeltaHeading = 0, DeltaRoll = 0, DeltaPitch = 0;// ошибки выставки по курсу, крену и тангажу соответственно
	Ldoub Heading = 0, Roll = 0, Pitch = 0;
	bool AlignmentContinue = true; // для начала выставки
	// Необходимые массивы для решение навигационной задачи
	Ldoub Ab[3]={0}; // Ускорения в связанных осях
	Ldoub Ao[3] = {0}; // Ускорения в географических осях
	Ldoub Omb[3] = {0}; // Угловые скорости в связанных осях
	Ldoub Omo[3] = {0}; // Угловые скорости в географических осях
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
	// массивы для выходных значений
	Ldoub V[3] = {(Ldoub) V0[0], (Ldoub) V0[1], 0}; // линейные скорости E; N; Up
	Ldoub Coordinates[3] = {phi0, lambda0, 0}; // Географические кординаты: широта, долгота и высота
	Ldoub CoordError[2] = {0}; // Ошибки в м (dE, dN)
	Ldoub Orientation[3] = {0}; // Углы ориентации

	Rlambda = (Ldoub) R/sqrt(1.-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
	Rphi = (Ldoub) R*(1. - e*e)/(sqrt(1.-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1.-e*e*sin(Coordinates[0])*sin(Coordinates[0])));

	bool AllowBiasAcc = false;
	bool AllowBiasGyr = false;
	bool AllowRandAcc = false;
	bool AllowRandGyr = false;
	bool AllowRandVgps = false;
	//Создаем квазикоординаты
	Ldoub alpha[12] = {0}; //малые приращения углов 3 показания на 4 тактах (матрица 3*4)
	Ldoub w[12] = {0}; // малые приращения скоростей (матрица 3*4)
	
	bool derectNorm (true); //направление ортогонализации и нормализации
	// Чтение из файла ускорений и угловых скоростей
#if 0
	QFile file(argv[1]);
	file.open(QIODevice::ReadOnly);
	QDataStream in(&file);
	in.setByteOrder(QDataStream::LittleEndian);
#endif

#if 1
	FILE* file=fopen(argv[5], "rt");
	fscanf(file, "%*s;");//чтение строки заголовка, она не нужна, выбрасываем
#endif
	//std::ifstream in(argv[5], std::ios::binary); !in.eof() //условие цикла while для бинарного файла

	Ldoub resultV[2] = {0};
	Ldoub Verr1[2] = {0}; // ошибки интегрирования ускорений 
	Ldoub Verr2[2] = {0}; // ошибки накопления скоростей

	AdaptiveKalman filter(7, 3); //Создаю объект Адаптивного фильтра Калмана с матрицей размера 7*7 и измерениями 3*1 (вертикальную скорость тоже учитываю)
	Ldoub x0[2] = {0}; //Начальные оценочные значения дрейфов
	Ldoub H[7*3] = {0}; //матрица наблюдения
	H[index_3(3, 0, 0)] = 1;
	H[index_3(3, 1, 1)] = 1;
	H[index_3(3, 2, 2)] = 1;
	//Матрца ковариации входных шумов (модели)
	Ldoub q[7*7] = {0};

	
	while( Readfile(file, AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps) == 22 ) // Пока возможно чтение из файла. Для бинарного файла -- !in.eof() //условие цикла while для бинарного файла
	{
		// этап выставки
		if ((cur_time <= t_alignment) && AlignmentContinue )
		{
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
			//ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps); //чтение из бинарного файла данных используемых для выставки
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
			}

		}
#if 1
		Ldoub MeanAlpha[3] = {0};//осредненные малые приращения углов (псевдокоординаты)
		Ldoub MeanW[3] = {0}; //осредненные малые приращения скоростей

		Ldoub Wp[3] = {0}; //проинтегрированные малые приращения. Начальные значения обнуляются на каждом такте быстрого цикла (с частотой 100 Гц)

		for (int in_iter=0; in_iter < 4; ++in_iter) // 4 такта, нумерация с нуля, поэтому равенство нестрогое
		{
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
			//ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, AllowRandVgps, Ab, Omb, Vgps);//чтение из бинарного файла данных используемых для навигации
		#ifdef dev
			printf("Ab[0] = %.10f; Ab[1] = %.10f; Ab[2] = %.10f\n", Ab[0], Ab[1], Ab[2]);
		#endif //dev
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
#ifdef dev
				printf("k1[%d] = %.10f\n", ii, k1[ii]);
				printf("k2[%d] = %.10f\n", ii, k2[ii]);
				printf("k3[%d] = %.10f\n", ii, k3[ii]);
				printf("k4[%d] = %.10f\n", ii, k4[ii]);
				printf("Wp[%d] = %.10f\n", ii, Wp[ii]);
#endif
			}
		}		
		/*Далее идет 100 Гц такт*/
		//Решение задачи ориентации
		Ldoub omo[3] = {0}; //переносные угловые скорости, вчисляются в решении задачи ориентации (SolveOrient)
		SolveOrient(alpha, Cib, Cin, Cbn, Orientation, Coordinates, Omo, omo, V, phi0, Rphi, Rlambda, freq, h, U, cur_time, derectNorm); //из одноименного заголовочного файла
		
		/*Решение задачи навигации*/
		SolveNav(Wp, Ab, Cbn,  Ao, V,  Coordinates, CoordError, Err_V, omo, h, Rphi, Rlambda, U, R, e, H0, V0 );
		// инкремент тактов
		++cur_time;

		/*
		По идее, куда-то сюда можно засунуть оценивание по Калману скоростей дрейфов гироскопов
		*/
	#if 0
		//Каждый такт пересчитываем матрицу A у фильтра Калмана
		// Delta dot V_ox
		filter.A[0] = V[1]/(R+Coordinates[2])*tan(Coordinates[0]) - V[2]/(R+Coordinates[2]);
		filter.A[1] = V[0]/(R+Coordinates[2])*tan(Coordinates[0]) + 2*U*sin(Coordinates[0]);
		filter.A[2] = -(V[0]/(R+Coordinates[2]) + 2*U*cos(Coordinates[0]));
		filter.A[3] = 0;
		filter.A[4] =  -Ao[2];
		filter.A[5] = 0;
		filter.A[6] = 0;
		//Delta dot V[1]
		filter.A[7] = -2*(V[0]/(R+Coordinates[2])*tan(Coordinates[0]) + U*sin(Coordinates[0]));
		filter.A[8] = -V[2]/(R+Coordinates[2]);
		filter.A[9] = -V[1]/(R+Coordinates[2]);
		filter.A[10] = Ao[2];
		filter.A[11] = 0;
		filter.A[12] = 0;
		filter.A[13] = 0;
		// Delta dot V[2]
		filter.A[14] = 2*(V[0]/(R+Coordinates[2]) + U*cos(Coordinates[0]));
		filter.A[15] = 2*V[1]/(R+Coordinates[2]);
		filter.A[16] = 0;
		filter.A[17] = -Ao[1];
		filter.A[18] = Ao[0];
		filter.A[19] = 0;
		filter.A[20] = 0;
		//Phi_ox
		filter.A[21] = 0;
		filter.A[22] = -1/(R+Coordinates[2]);
		filter.A[23] = 0;
		filter.A[24] = 0;
		filter.A[25] = omo[2];
		filter.A[26] = -cos(Orientation[0]);
		filter.A[27] = -sin(Orientation[0]);
		//Phi_oy
		filter.A[28] = 1/(R+Coordinates[2]);
		filter.A[29] = 0;
		filter.A[30] = 0;
		filter.A[31] = - omo[2];
		filter.A[32] = 0;
		filter.A[33] = sin(Orientation[0]);
		filter.A[34] = -cos(Orientation[0]);
		//Delta omega_x
		filter.A[35] = 0;
		filter.A[36] = 0;
		filter.A[37] = 0;
		filter.A[38] = 0;
		filter.A[39] = 0;
		filter.A[40] = 0;
		filter.A[41] = 0;
		//Delta omega_y
		filter.A[42] = 0;
		filter.A[43] = 0;
		filter.A[44] = 0;
		filter.A[45] = 0;
		filter.A[46] = 0;
		filter.A[47] = 0;
		filter.A[48] = 0;

		if (!filter.init) //Если ранее не было инициализации, то инициализируем
			filter.Init(x0, q, H);
		else//в противном случае оцениваем
		{
			for (int iii =0; iii< filter.getDimX(); ++iii) //вычисляю матрицу перехода Phi
				filter.Phi[iii] = filter.I[iii] + filter.A[iii];
			
			Ldoub ErrVins[3] = {V[0] - Vgps[0], V[1] - Vgps[1], 0}; //разница ошибок БИНС и СНС
			filter.Predict();
			filter.Update(ErrVins);
			printf("omega_x = %.8f \t omega_y = %.8f\n", filter.x[0], filter.x[1]);
		}
	#endif
		// Запись в файл
		
		static FILE* navig_res;
		if(!navig_res)
		{
			if (argc == 7) //если аргументом передан файл для записи результатов, записываем в него
					navig_res=fopen(argv[6],"wt");
			else // если нет, записываем в файл по умолчанию
				navig_res=fopen("./data/Nav_res.csv","wt");
			// Шапка
			fprintf(navig_res, "Ve;");
			fprintf(navig_res, "Vn;");
			fprintf(navig_res, "Vup;");
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
			for(int i=0; i<3; ++i)
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
	}
	return 0;
}
