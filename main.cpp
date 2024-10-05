#include <stdio.h>
#include <matrix.h>
#include <math.h>
#include <std.h>

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
	const double rad2deg = 180/pi; // из градусов в час в радианы в секунду
	const double deg2rad = 1/rad2deg;
	int freq = 100; // частота измерений с инерциальных датчиков
	int t_nav = 90*60; // время работы нав алгоритма в секундах
	int t_alignment = 5*60*freq; // время выставки в тактах

	int cur_time = 0; // текущий такт!! измерения
	//Необходимое для выставки
	double phi0 = 55*deg2rad;
	double lambda0 = 33*deg2rad;
	double DeltaHeading = 0, DeltaRoll = 0, DeltaPitch = 0;// ошибки выставки по курсу, крену и тангажу соответственно
	double Heading = 0, Roll = 0, Pitch = 0;
	bool AlignmentContinue = true; // для начала выставки
	// Необходимые массивы для решение навигационной задачи
	double Ab[3] = {0}; // Ускорения в связанных осях
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
	double V[3] = {100*sin(50*deg2rad), 100*cos(50*deg2rad), 0}; // линейные скорости E; N; Up
	double Coordinates[3] = {phi0, lambda0, 0}; // Географические кординаты: широта, долгота и высота
	double Orientation[3] = {0}; // Углы ориентации
	double Rlambda;
	Rlambda = R/sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
	double Rphi;
	Rphi = R*(1 - e*e)/(sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1-e*e*sin(Coordinates[0])*sin(Coordinates[0])));
	// Чтение из файла ускорений и угловых скоростей
	FILE* file=fopen("/home/nikita_viksne/Документы/Python/Modelling_sensetive_elements/Data_files/data_acc.csv", "rt");
	while(true)
	{
		// int res = fscanf(file, "%f;%f;%f;%f;%f;%f;", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2]);
		int res = fscanf(file, "%lf%lf%lf%lf%lf%lf", &Ab[0],&Ab[1],&Ab[2],&Omb[0],&Omb[1],&Omb[2]);
		if (res!=6)
		{
			printf("Check error of End of file\n");
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
				double c0 = sqrt(Cbn[index(3, 2, 0)]* Cbn[index(3, 2, 0)] + Cbn[index(3, 2, 2)]*Cbn[index(3, 2, 2)]);
				// Вычисление углов ориентации
				Heading = atan2(Cbn[index(3, 0, 1)], Cbn[index(3, 1, 1)]);
				Roll = - atan2(Cbn[index(3, 2, 0)], Cbn[index(3, 2, 2)]);
				Pitch = atan2(Cbn[index(3, 2, 1)], c0);
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
				printf("Heading (degrees) = %f Error = %f\n", Heading*rad2deg, DeltaHeading*rad2deg);
				printf("roll (degrees) = %f Error = %f\n", Roll*rad2deg, DeltaRoll*rad2deg);
				printf("pitch (degrees) = %f Error = %f\n", Pitch*rad2deg, DeltaPitch*rad2deg);
				AlignmentContinue = false;
			}
		}
#if 1
		Omo[0] = -V[1]/(Rphi + Coordinates[2]);
		Omo[1] = V[0]/(Rlambda + Coordinates[2]) + U*cos(Coordinates[0]);
		Omo[2] = V[0]/(Rlambda + Coordinates[2])*tan(Coordinates[0]) + U*sin(Coordinates[0]);
		Coordinates[0] += (V[1]/(Rphi + Coordinates[2]))/freq;
		Coordinates[1] += (V[0]/((Rlambda + Coordinates[2])*cos(Coordinates[0])))/freq;
		Coordinates[2] += (V[2])/freq;
		// Решение уравнения Пуассона
		double EigWb[9] = {0, -Omb[2], Omb[1], Omb[2], 0, -Omb[0], -Omb[1], Omb[0], 0};
		double EigWo[9] = {0, -Omo[2], Omo[1], Omo[2], 0, -Omo[0], -Omo[1], Omo[0], 0};
		double CEigWb[9]; // первое слагаемое уравнения Пуассона
		double EigWoC[9]; // второе слагаемое уравнения Пуассона
 		MulMatrD(Cbn, EigWb, CEigWb,3,3,3);
		MulMatrD(EigWo, Cbn, EigWoC,3,3,3);

		for(int i=0; i<9; ++i)
			Cbn[i] += (CEigWb[i] - EigWoC[i])/freq;
		// вычисление углов ориентации через МНК (как в выставке)
		Orientation[0] = atan2(Cbn[index(3, 0, 1)], Cbn[index(3, 1, 1)]);
		Orientation[1] = - atan2(Cbn[index(3, 2, 0)], Cbn[index(3, 2, 2)]);
		double c0 = sqrt(Cbn[index(3, 2, 0)]* Cbn[index(3, 2, 0)] + Cbn[index(3, 2, 2)]*Cbn[index(3, 2, 2)]);
		Orientation[2] = atan2(Cbn[index(3, 2, 1)], c0);

		// решение задачи навигации
		MulMatrD(Cbn, Ab, Ao,3,3,1); // перепроектирование из связаных осей в навигационные
		//V[2] += (Ao[2] + (Omo[1] + U*cos(Coordinates[0]))*V[0] + V[1]*Omo[0] - g*(1-2*Coordinates[2]/Rphi))/freq; //Vup
		V[0] += (Ao[0] + (U*sin(Coordinates[0]) + Omo[2]) * V[1] - V[2] *(U * cos(Coordinates[0]) + Omo[1]))/freq;// Ve
		V[1] += (Ao[1] - (U*sin(Coordinates[0]) + Omo[2]) * V[0] + V[2] * Omo[0])/freq; //Vn

		Rlambda = R/sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0]));
		Rphi = R*(1 - e*e)/(sqrt(1-e*e*sin(Coordinates[0])*sin(Coordinates[0])) * (1-e*e*sin(Coordinates[0])*sin(Coordinates[0])));

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
			fprintf(navig_res, "\n");
		}

		if(navig_res)
		{
			// Навигационные параметры
			// Скорости
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%e;", V[i]);
			// Координаты
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%e;", Coordinates[i]);
			// Углы оориентации
			for(int i=0; i<3; ++i)
				fprintf(navig_res, "%e;", Orientation[i]);

			fprintf(navig_res, "\n");
		}
#endif
	}
	//getc(stdin);
	return 0;
}
