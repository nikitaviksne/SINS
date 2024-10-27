#if 1
#define _USE_MATH_DEFINES
#include <stdio.h>
#include "matrix.h"
#if 0
#include <math.h>
#endif
#include "std.h"
#include <cmath>
#include "mathematics.h"

int main()
{
	//инициализация необходимых переменных и констант
	const double g = 9.81;
	const double a = 6378245;
	const double b = 6356856;
	const double e=sqrt(1 - b*b/a/a);
	const double R=6400e3;
	const double pi=3.141592653589793;
	const double U = 7.27220521664304e-05;
	const double rad2deg = 180./M_PI; // из градусов в час в радианы в секунду
	const double deg2rad = 1./rad2deg;

	int freq = 100; // частота измерений с инерциальных датчиков
	double h = 0.01; //период дискретизации
	int t_nav = 90*60; // время работы нав алгоритма в секундах
	int t_alignment = 5*60*freq; // время выставки в тактах

	int cur_time = 0; // текущий такт!! измерения
	/*
	 * Начальные значения.
	 *Чтобы расчитать начальные значения скоростей, с которых проводить интегрирование
	*/
	// Для моделирования показаний Ч.Э.
	double H0 = (double) (50)*deg2rad;
	double P0 = (double) (0)*deg2rad;
	double R0 = (double) (0)*deg2rad;
	double Vabs = 300;
	double Cnb[9];
	MatrOB(H0, R0, P0, Cnb, 3); // матрица перехода из опорной в связанную
	//Необходимое для выставки
	double phi0 = (double) 55*deg2rad;
	double lambda0 = (double) 33*deg2rad;
	double DeltaHeading = 0, DeltaRoll = 0, DeltaPitch = 0;// ошибки выставки по курсу, крену и тангажу соответственно
	double Heading = 0, Roll = 0, Pitch = 0;
	bool AlignmentContinue = true; // для начала выставки
	// Необходимые массивы для решение навигационной задачи
	double Ab[3] = {0}; // Ускорения в связанных осях
	double BiasAb[3] = {0}; //Смещение нулей акселерометров
	double BiasOmb[3] = {0}; //Смещение нулей гироскопов
	double RandAb[3] = {0}; //Случайные погрешности акселерометров
	double RandOmb[3] = {0}; //Случайные погрешности гироскопов
	double Ao[3] = {0}; // Ускорения в географических осях
	double Omb[3] = {0}; // Угловые скорости в связанных осях
	double Omo[3] = {0}; // Угловые скорости в географических осях
	//double Oms[3] = {0}; // угловые скорости от линейного движения + Земля
	double MeanAb[3] = {0};
	double MeanOmb[3] = {0};
	double StdAb[3] = {0};
	double StdOmb[3] = {0};
	double Cbn[9] = {0};
	// массивы для выходных значений
	double V[3] = {(double) Vabs*sin(H0), (double) Vabs*cos(H0), 0}; // линейные скорости E; N; Up
	double Coordinates[3] = {phi0, lambda0, 0}; // Географические кординаты: широта, долгота и высота
	double Orientation[3] = {0}; // Углы ориентации
	double CoordErr[2] = {0}; // Ошибки по координатам в метрах
	double Rlambda;

	double sqrEErr {0}; //Ошибка возведения e в квадрат 
	double E2E {0};
	TwoProduct(e, e, E2E, sqrEErr);
	TwoSum(E2E, sqrEErr, E2E, sqrEErr, false);
	Rlambda = (double) R/sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
	double Rphi;
	Rphi = (double) R*(1 - e*e)/(sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1-e*e*sin(Coordinates[0])*sin(Coordinates[0])));

	bool AllowBiasAcc = false;
	bool AllowBiasGyr = false;
	bool AllowRandAcc = false;
	bool AllowRandGyr = false;
	// Чтение из файла ускорений и угловых скоростей
	FILE* file=fopen("/home/nikita_viksne/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc.csv", "rt");
	fscanf(file, "%*[^\n]"); //для файла, разделенного пробелами
	//fscanf(file, "%*s;"); // для файла, разделенного точкой с запятой

	double resultV[2] = {0};
	double Verr1[2] = {0}; // ошибки интегрирования ускорений
	double Verr2[2] = {0}; // ошибки накопления скоростей

	while(true)
	{
		// int res = fscanf(file, "%f;%f;%f;%f;%f;%f;", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2]);
		//int res = fscanf(file, "%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2], &BiasAb[0],&BiasAb[1],&BiasAb[2],&BiasOmb[0],&BiasOmb[1], &BiasOmb[2], &RandAb[0], &RandAb[1],&RandAb[2],&RandOmb[0],&RandOmb[1],&RandOmb[2]);
		int res = fscanf(file, "%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e;", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2], &BiasAb[0],&BiasAb[1],&BiasAb[2],&BiasOmb[0],&BiasOmb[1], &BiasOmb[2], &RandAb[0], &RandAb[1],&RandAb[2],&RandOmb[0],&RandOmb[1],&RandOmb[2]);
		if (res!=18)
		{
			printf("Check error at %d iteration\n", cur_time);
			break;
		}
#if 1
		/*Генерирование (моделирование) показаний ч.э*/
		if (cur_time <= t_alignment)
		{
			Ao[0] = 0.0; Ao[1] = 0.0; Ao[2] = g;
			Omo[0] = 0;
			Omo[1] = (double) U*cos(phi0);
			Omo[2] = (double) U*sin(phi0);
		}
		else
		{
			Omo[0] = (double) -Vabs*cos(H0) / (R + 0); //Rphi
			Omo[1] = (double) Vabs*sin(H0) / ((R + 0)) + (double) U*cos(phi0); // Rlambda
			Omo[2] = (double) Vabs*sin(H0) * tan(phi0) / (R + 0) + (double) U*sin(phi0);// Rlambda
			Ao[0] = 0.0;//(double) ( Omo[1]*0 -(double) Omo[2]*Vabs*cos(H0) + (double) U*cos(phi0)*0 - (double) U*sin(phi0)*Vabs*cos(H0));
			Ao[1] = 0.0;//(double) (-Omo[0]*0 +(double) Omo[2]*Vabs*sin(H0) + (double) U*sin(phi0)*Vabs*sin(H0));
			Ao[2] = g;
			phi0 += (double) Vabs*cos(H0)/(R + 0) * h;
		}
		//printf("Cur_time = %d\n", cur_time);
		//printf("Modelled Ao = [%.20f; %.20f; %.20f]\n", Ao[0], Ao[1], Ao[2]);
		MulMatrD(Cnb, Ao, Ab,3,3,1); // проекция ускорений на связанные оси
		//printf("Modelled Ab = [%.20e; %.20e; %.20e]\n", Ab[0], Ab[1], Ab[2]);
		MulMatrD(Cnb, Omo, Omb,3,3,1); // проекция угловых скоростей на связанные оси
#endif
		// этап выставки
		if (cur_time <= t_alignment)
		{
			for(int i=0; i<3; ++i)
			{
				// применяем метод Уэлфорда
				MeanAb[i] = (double) MeanAb[i] + (Ab[i] - MeanAb[i]) / (cur_time+1);
				StdAb[i] = (double) (1 - 1/(cur_time + 1))*StdAb[i] + (Ab[i] - MeanAb[i])*(Ab[i] - MeanAb[i])/(cur_time + 1);
				MeanOmb[i] = (double) MeanOmb[i] + (Omb[i] - MeanOmb[i]) / (cur_time+1);
				StdOmb[i] = (double) (1 - 1/(cur_time + 1))*StdOmb[i] + (Omb[i] - MeanOmb[i])*(Omb[i] - MeanOmb[i])/(cur_time + 1);
			}
			//Вычисление (ориентации) матрицы перехода Cbn = [c00, c01, c02, c10, c11, c12, c20, c21, c22]
			// Ищем обратную (транспонированную) матрицу
			for(int i=0; i<3; ++i)
			{
				Cbn[index(3, 2, i)] = (double) MeanAb[i] / g;
				Cbn[index(3, 1, i)] = (double) (MeanOmb[i]/ U - (double) MeanAb[i]/g*sin(phi0))/cos(phi0);
			}
			// по алгебраическому дополнению
			Cbn[index(3, 0, 0)] = (double) Cbn[index(3, 1, 1)] * Cbn[index(3, 2, 2)] - Cbn[index(3, 2, 1)] * Cbn[index(3, 1, 2)];
			Cbn[index(3, 0, 1)] = (double) - Cbn[index(3, 1, 0)] * Cbn[index(3, 2, 2)] + Cbn[index(3, 1, 2)] * Cbn[index(3, 2, 0)];
			Cbn[index(3, 0, 2)] = (double) Cbn[index(3, 1, 0)] * Cbn[index(3, 2, 1)] - Cbn[index(3, 1, 1)] * Cbn[index(3, 2, 0)];
			++cur_time;
			continue;
		}
		else
		{
			if (AlignmentContinue)
			{
				double c0 = (double) sqrt(Cbn[index(3, 2, 0)]* Cbn[index(3, 2, 0)] + Cbn[index(3, 2, 2)]*Cbn[index(3, 2, 2)]);
				// Вычисление углов ориентации
				Heading = (double) atan2(Cbn[index(3, 0, 1)], Cbn[index(3, 1, 1)]);
				Roll = (double) - atan2(Cbn[index(3, 2, 0)], Cbn[index(3, 2, 2)]);
				Pitch = (double) atan2(Cbn[index(3, 2, 1)], c0);
				// вычисление ошибок выставки
				DeltaRoll = (StdAb[0] * MeanAb[2] - StdAb[2] * MeanAb[0])/(MeanAb[2]*MeanAb[2]
				+ MeanAb[0]*MeanAb[0]);
				DeltaPitch = (StdAb[1])/sqrt(g*g - MeanAb[1]*MeanAb[1]);
				double DeltaWn[3];
				MulMatrD(Cbn, StdOmb, DeltaWn,3,3,1); // проекция дрейфов гироскопов на географические оси
				double DeltaAn[3];
				MulMatrD(Cbn, StdAb, DeltaAn,3,3,1); // проекция дрейфов акселерометров на географические оси
				DeltaHeading = - DeltaWn[0]/(U*cos(phi0)) + DeltaAn[0]/g*tan(phi0) - StdAb[2]/2/g*sin(2*Heading);
				printf("Alignment\n");
				printf("Heading (degrees) = %.30f Error = %f\n", (double) Heading*rad2deg, (double) DeltaHeading*rad2deg);
				printf("roll (degrees) = %.30f Error = %f\n", (double) Roll*rad2deg, (double) DeltaRoll*rad2deg);
				printf("pitch (degrees) = %.30f Error = %f\n", (double) Pitch*rad2deg, (double) DeltaPitch*rad2deg);
				AlignmentContinue = false;
			}
		}
#if 1
		Omo[0] = (double) -V[1]/(R + Coordinates[2]);// Rphi
		Omo[1] = (double) V[0]/(R + Coordinates[2]);// Rlambda
		Omo[2] = (double) V[0]/(R + Coordinates[2])*tan(Coordinates[0]);// Rlambda
		double omo[3] = {(double) Omo[0], (double) Omo[1] + (double) U*cos(Coordinates[0]), (double) Omo[2] + (double) U*sin(Coordinates[0])}; // абсолютные угловые скорости
		Coordinates[0] += (double) (V[1]/(R + Coordinates[2])) * h; // Rphi
		Coordinates[1] += (double) (V[0]/((R + Coordinates[2])*cos(Coordinates[0]))) * h; // Rlambda
		Coordinates[2] += ((double) V[2]) * h;
		// Решение уравнения Пуассона
		double EigWb[9] = {0, -Omb[2], Omb[1], Omb[2], 0, -Omb[0], -Omb[1], Omb[0], 0};
		double EigWo[9] = {0, -omo[2], omo[1], omo[2], 0, -omo[0], -omo[1], omo[0], 0};
		double CEigWb[9]; // первое слагаемое уравнения Пуассона
		double EigWoC[9]; // второе слагаемое уравнения Пуассона
 		MulMatrD(Cbn, EigWb, CEigWb,3,3,3);
		MulMatrD(EigWo, Cbn, EigWoC,3,3,3);

		for(int i=0; i<9; ++i)
			Cbn[i] += (CEigWb[i] - EigWoC[i]) * h;
		// вычисление углов ориентации через МНК (как в выставке)
		Orientation[0] = (double) atan2(Cbn[index(3, 0, 1)], Cbn[index(3, 1, 1)]);
		Orientation[1] = (double) - atan2(Cbn[index(3, 2, 0)], Cbn[index(3, 2, 2)]);
		double c0 = (double) sqrt(Cbn[index(3, 2, 0)]* Cbn[index(3, 2, 0)] + Cbn[index(3, 2, 2)]*Cbn[index(3, 2, 2)]);
		Orientation[2] = (double) atan2(Cbn[index(3, 2, 1)], c0);

		// решение задачи навигации
		MulMatrD(Cbn, Ab, Ao,3,3,1); // перепроектирование из связаных осей в навигационные
		//V[2] += (Ao[2] + (Omo[1] + U*cos(Coordinates[0]))*V[0] + V[1]*Omo[0] - g*(1-2*Coordinates[2]/Rphi)) * h; //Vup
#if 1
		//Кориолисовы добавки
		double aCoriolis[3] = {0};
		aCoriolis[0] =(double) ((double) omo[1]*V[2] - (double) omo[2]*V[1] + (double) U*cos(Coordinates[0])*V[2] - (double) U*sin(Coordinates[0])*V[1]);
		aCoriolis[1] =(double) ((double) -omo[0]*V[2] + (double) omo[2]*V[0] + (double) U*sin(Coordinates[0])*V[0]);
		aCoriolis[2] =(double) ((double) omo[0]*V[1] - (double) omo[1]*V[0] - (double) U*cos(Coordinates[0])*V[0]);

		double V_dot[2] = {(double) Ao[0], (double) Ao[1] };

		/*
		V[0] += (double) (Ao[0] - aCoriolis[0]) * h;// Ve
		V[1] += (double) (Ao[1] - aCoriolis[1]) * h; //Vn

		V[0] += (double) (Ao[0]) * h;// Ve
		V[1] += (double) (Ao[1]) * h; //Vn
		*/

		// умножение
		TwoProduct(Ao[0],h, V_dot[0], Verr2[0]);
		TwoProduct(Ao[1],h, V_dot[1], Verr2[1]);
		// компенсация ошибок интегрирования ускорений
		TwoSum(Verr2[0], V_dot[0], V_dot[0], Verr2[0], false);
		TwoSum(Verr2[1], V_dot[1], V_dot[1], Verr2[1], false);
		// сложение скоростей с предыдущего такта и только что проинтегрированных ускорений (приращения скоростей). Оценка погрешности этого сложения
		TwoSum(V_dot[0], V[0], V[0], Verr1[0], false);
		TwoSum(V_dot[1], V[1],  V[1], Verr1[1], false);
		//Компенсация погрешности сложения скоростей с пред. такта и приращения скоростей на тек. такте
		TwoSum(Verr1[0], V[0], V[0], Verr1[0], false);
		TwoSum(Verr1[1], V[1],  V[1], Verr1[1], false);

		CoordErr[0] +=(V[0] - (double) Vabs*sin(H0)) * h;
		CoordErr[1] +=(V[1] - (double) Vabs*cos(H0)) * h;
#endif
#if 1 // Пересчет радиусов сильно влияет на ошибки
#if 0
		double sqrSinErr {0}; //Ошибка возведения синуса в квадрат 
		double sin2sin {0};

		TwoProduct(sin(Coordinates[0]), sin(Coordinates[0]), sin2sin, sqrSinErr); 
		TwoSum(sin2sin, sqrSinErr, sin2sin, sqrSinErr, false); //sin*sin
		
		double e2Sin {0}; // e*e*sin*sin
		double e2SinErr {0};
		TwoProduct(E2E, sin2sin, e2Sin, e2SinErr); //// e*e*sin*sin
		TwoSum(e2Sin, e2SinErr, e2Sin, e2SinErr, false); //e*e*sin*sin

		e2SinErr = 0;
		double oneE2Sin {0}; // 1 - e*e*sin*sin
		TwoSum(1., -e2SinErr, oneE2Sin, e2SinErr, false);
		TwoSum(oneE2Sin, e2SinErr, oneE2Sin, e2SinErr, false);
		
		Rlambda = (double) R/sqrt(oneE2Sin);
		Rphi = (double) R*(1 - e*e)/(sqrt(oneE2Sin) * (oneE2Sin));
#endif

#if 1
		Rlambda = (double) R/sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
		Rphi = (double) R*(1 - e*e)/(sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1-e*e*sin(Coordinates[0])*sin(Coordinates[0])));
#endif
#endif

		// инкремент тактов
		++cur_time;

		// Запись в файл

		static FILE* navig_res;
		if(!navig_res)
		{
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

			//Ошибки по координатам 
			for(int iii=0; iii<2; ++iii)
				fprintf(navig_res, "%.10e;", CoordErr[iii]);
			fprintf(navig_res, "\n");
		}
#endif
	}
	printf("File end\n");
	//getc(stdin);
	return 0;
}
#endif