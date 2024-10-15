#include <stdio.h>
#include <matrix.h>
#if 0
#include <math.h>
#endif
#include <std.h>
#include <cmath>

int main()
{
	//инициализация необходимых переменных и констант
	const float g = 9.81;
	const float a = 6378245;
	const float b = 6356856;
	const float e=sqrt(1 - b*b/a/a);
	const float R=6400e3;
	const float pi=3.141592653589793;
	const float U = 7.27220521664304e-05;
	const float rad2deg = 180/pi; // из градусов в час в радианы в секунду
	const float deg2rad = 1/rad2deg;
	/*
	 * Начальные значения.
	 *Чтобы расчитать начальные значения скоростей, с которых проводить интегрирование
	 */
	float Vabs = 300; // Модуль линейной скорости
	float H0 = 50*deg2rad; // Начальный угол курса в радианах

	int freq = 100; // частота измерений с инерциальных датчиков
	int t_nav = 90*60; // время работы нав алгоритма в секундах
	int t_alignment = 5*60*freq; // время выставки в тактах

	int cur_time = 0; // текущий такт!! измерения
	//Необходимое для выставки
	float phi0 = (float) 55*deg2rad;
	float lambda0 = (float) 33*deg2rad;
	float DeltaHeading = 0, DeltaRoll = 0, DeltaPitch = 0;// ошибки выставки по курсу, крену и тангажу соответственно
	float Heading = 0, Roll = 0, Pitch = 0;
	bool AlignmentContinue = true; // для начала выставки
	// Необходимые массивы для решение навигационной задачи
	float Ab[3] = {0}; // Ускорения в связанных осях
	float BiasAb[3] = {0}; //Смещение нулей акселерометров
	float BiasOmb[3] = {0}; //Смещение нулей гироскопов
	float RandAb[3] = {0}; //Случайные погрешности акселерометров
	float RandOmb[3] = {0}; //Случайные погрешности гироскопов
	float Ao[3] = {0}; // Ускорения в географических осях
	float Omb[3] = {0}; // Угловые скорости в связанных осях
	float Omo[3] = {0}; // Угловые скорости в географических осях
	//float Oms[3] = {0}; // угловые скорости от линейного движения + Земля
	float MeanAb[3] = {0};
	float MeanOmb[3] = {0};
	float StdAb[3] = {0};
	float StdOmb[3] = {0};
	float Cbn[9] = {0};
	// массивы для выходных значений
	float V[3] = {(float) Vabs*sin(H0), (float) Vabs*cos(H0), 0}; // линейные скорости E; N; Up
	float Coordinates[3] = {phi0, lambda0, 0}; // Географические кординаты: широта, долгота и высота
	float Orientation[3] = {0}; // Углы ориентации
	float CoordErr[2] = {0}; // Ошибки по координатам в метрах
	float Rlambda;
	Rlambda = (float) R/sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
	float Rphi;
	Rphi = (float) R*(1 - e*e)/(sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1-e*e*sin(Coordinates[0])*sin(Coordinates[0])));
	// Чтение из файла ускорений и угловых скоростей
	FILE* file=fopen("/home/nikita_viksne/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc.csv", "rt");
	fscanf(file, "%*[^\n]");
	while(true)
	{
		// int res = fscanf(file, "%f;%f;%f;%f;%f;%f;", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2]);
		int res = fscanf(file, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2], &BiasAb[0], &BiasAb[1], &BiasAb[2], &BiasOmb[1], &BiasOmb[2], &BiasOmb[3], &RandAb[0], &RandAb[1], &RandAb[2], &RandOmb[0], &RandOmb[1], &RandOmb[2]);
		if (res!=18)
		{
			printf("Check error at %d iteration\n", cur_time);
			break;
		}
		// этап выставки
		if (cur_time <= t_alignment)
		{
			for(int i=0; i<3; ++i)
			{
				// применяем метод Уэлфорда
				MeanAb[i] = MeanAb[i] + (Ab[i] - MeanAb[i]) / (cur_time+1);
				StdAb[i] = (1 - 1/(cur_time + 1))*StdAb[i] + (Ab[i] - MeanAb[i])*(Ab[i] - MeanAb[i])/(cur_time + 1);
				MeanOmb[i] = MeanOmb[i] + (Omb[i] - MeanOmb[i]) / (cur_time+1);
				StdOmb[i] = (1 - 1/(cur_time + 1))*StdOmb[i] + (Omb[i] - MeanOmb[i])*(Omb[i] - MeanOmb[i])/(cur_time + 1);
			}
			//Вычисление (ориентации) матрицы перехода Cbn = [c00, c01, c02, c10, c11, c12, c20, c21, c22]
			// Ищем обратную (транспонированную) матрицу
			for(int i=0; i<3; ++i)
			{
				Cbn[index(3, 2, i)] = MeanAb[i] / g;
				Cbn[index(3, 1, i)] = (MeanOmb[i]/ U - MeanAb[i]/g*sin(phi0))/cos(phi0);
			}
			// по алгебраическому дополнению
			Cbn[index(3, 0, 0)] = Cbn[index(3, 1, 1)] * Cbn[index(3, 2, 2)] - Cbn[index(3, 2, 1)] * Cbn[index(3, 1, 2)];
			Cbn[index(3, 0, 1)] = - Cbn[index(3, 1, 0)] * Cbn[index(3, 2, 2)] + Cbn[index(3, 1, 2)] * Cbn[index(3, 2, 0)];
			Cbn[index(3, 0, 2)] = Cbn[index(3, 1, 0)] * Cbn[index(3, 2, 1)] - Cbn[index(3, 1, 1)] * Cbn[index(3, 2, 0)];
			++cur_time;
			continue;
		}
		else
		{
			if (AlignmentContinue)
			{
				float c0 = (float) sqrt(Cbn[index(3, 2, 0)]* Cbn[index(3, 2, 0)] + Cbn[index(3, 2, 2)]*Cbn[index(3, 2, 2)]);
				// Вычисление углов ориентации
				Heading = (float) atan2(Cbn[index(3, 0, 1)], Cbn[index(3, 1, 1)]);
				Roll = (float) - atan2(Cbn[index(3, 2, 0)], Cbn[index(3, 2, 2)]);
				Pitch = (float) atan2(Cbn[index(3, 2, 1)], c0);
				// вычисление ошибок выставки
				DeltaRoll = (StdAb[0] * MeanAb[2] - StdAb[2] * MeanAb[0])/(MeanAb[2]*MeanAb[2]
				+ MeanAb[0]*MeanAb[0]);
				DeltaPitch = (StdAb[1])/sqrt(g*g - MeanAb[1]*MeanAb[1]);
				float DeltaWn[3];
				MulMatrD(Cbn, StdOmb, DeltaWn,3,3,1); // проекция дрейфов гироскопов на географические оси
				float DeltaAn[3];
				MulMatrD(Cbn, StdAb, DeltaAn,3,3,1); // проекция дрейфов акселерометров на географические оси
				DeltaHeading = - DeltaWn[0]/(U*cos(phi0)) + DeltaAn[0]/g*tan(phi0) - StdAb[2]/2/g*sin(2*Heading);
				printf("Alignment\n");
				printf("Heading (degrees) = %f Error = %f\n", Heading*rad2deg, DeltaHeading*rad2deg);
				printf("roll (degrees) = %f Error = %f\n", Roll*rad2deg, DeltaRoll*rad2deg);
				printf("pitch (degrees) = %f Error = %f\n", Pitch*rad2deg, DeltaPitch*rad2deg);
				AlignmentContinue = false;
			}
		}
#if 1
		Omo[0] = (float) -V[1]/(Rphi + Coordinates[2]);
		Omo[1] = (float) V[0]/(Rlambda + Coordinates[2]);
		Omo[2] = (float) V[0]/(Rlambda + Coordinates[2])*tan(Coordinates[0]);
		float omo[3] = {Omo[0], Omo[1] + U*cos(Coordinates[0]), Omo[2] + U*sin(Coordinates[0])}; // абсолютные угловые скорости
		Coordinates[0] += (float) (V[1]/(Rphi + Coordinates[2]))/freq;
		Coordinates[1] += (float) (V[0]/((Rlambda + Coordinates[2])*cos(Coordinates[0])))/freq; //
		Coordinates[2] += ((float) V[2])/freq;
		// Решение уравнения Пуассона
		float EigWb[9] = {0, -Omb[2], Omb[1], Omb[2], 0, -Omb[0], -Omb[1], Omb[0], 0};
		float EigWo[9] = {0, -omo[2], omo[1], omo[2], 0, -omo[0], -omo[1], omo[0], 0};
		float CEigWb[9]; // первое слагаемое уравнения Пуассона
		float EigWoC[9]; // второе слагаемое уравнения Пуассона
 		MulMatrD(Cbn, EigWb, CEigWb,3,3,3);
		MulMatrD(EigWo, Cbn, EigWoC,3,3,3);

		for(int i=0; i<9; ++i)
			Cbn[i] += (CEigWb[i] - EigWoC[i])/freq;
		// вычисление углов ориентации через МНК (как в выставке)
		Orientation[0] = (float) atan2(Cbn[index(3, 0, 1)], Cbn[index(3, 1, 1)]);
		Orientation[1] = (float) - atan2(Cbn[index(3, 2, 0)], Cbn[index(3, 2, 2)]);
		float c0 = (float) sqrt(Cbn[index(3, 2, 0)]* Cbn[index(3, 2, 0)] + Cbn[index(3, 2, 2)]*Cbn[index(3, 2, 2)]);
		Orientation[2] = (float) atan2(Cbn[index(3, 2, 1)], c0);

		// решение задачи навигации
		MulMatrD(Cbn, Ab, Ao,3,3,1); // перепроектирование из связаных осей в навигационные
		//V[2] += (Ao[2] + (Omo[1] + U*cos(Coordinates[0]))*V[0] + V[1]*Omo[0] - g*(1-2*Coordinates[2]/Rphi))/freq; //Vup
#if 1
		//Кориолисовы добавки
		float aCoriolis[3] = {0};
		aCoriolis[0] = omo[1]*V[2] - omo[2]*V[1] + U*cos(Coordinates[0])*V[2] - U*sin(Coordinates[0])*V[1];
		aCoriolis[1] = -omo[0]*V[2] + omo[2]*V[0] + U*sin(Coordinates[0])*V[0];
		aCoriolis[2] = omo[0]*V[1] - omo[1]*V[0] - U*cos(Coordinates[0])*V[0];

		V[0] += (float) (Ao[0] - aCoriolis[0])/freq;// Ve
		V[1] += (float) (Ao[1] - aCoriolis[1])/freq; //Vn

		CoordErr[0] +=(V[0] - Vabs*sin(H0))/freq;
		CoordErr[1] +=(V[1] - Vabs*cos(H0))/freq;
#endif

		Rlambda = (float) R/sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
		Rphi = (float) R*(1 - e*e)/(sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1-e*e*sin(Coordinates[0])*sin(Coordinates[0])));

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
