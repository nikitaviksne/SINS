#define _USE_MATH_DEFINES
#include <stdio.h>
#include "std.h"
#include "matrix.h"
#include "mathematics.h"
#include <cmath>
#include <fstream> //файловые потоки для чтения бинарного файла
#include <QDataStream>
#include <QFile>

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

void ReadFile(QDataStream &in, bool AllowBiasAcc, bool AllowBiasGyr, bool AllowRandAcc, bool AllowRandGyr, Ldoub* Ab, Ldoub* Omb)
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

}
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
int main()
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
	Ldoub H0 = (Ldoub) (0.)*deg2rad;
	Ldoub P0 = (Ldoub) (0.)*deg2rad;
	Ldoub R0 = (Ldoub) (0.)*deg2rad;
	Ldoub Vabs = 30;
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
	//Создаем квазикоординаты
	Ldoub alpha[12] = {0}; //малые приращения углов 3 показания на 4 тактах (матрица 3*4)
	Ldoub w[12] = {0}; // малые приращения скоростей (матрица 3*4)
	
	// Чтение из файла ускорений и угловых скоростей
#if 1
	QFile file("/home/nikita_viksne/Modelling_sensetive_elements/Data_files/data_acc_veloc_30_heading_0_freq_400.csv");
	file.open(QIODevice::ReadOnly);
	QDataStream in(&file);
	in.setByteOrder(QDataStream::LittleEndian);
#endif

#if 0
	FILE* file=fopen("C:/Users/Viksne_NA/Documents/Python/data_files/data_acc.csv", "rt");
	fscanf(file, "%*s;");
#endif
	Ldoub resultV[2] = {0};
	Ldoub Verr1[2] = {0}; // ошибки интегрирования ускорений 
	Ldoub Verr2[2] = {0}; // ошибки накопления скоростей
	
	while( !( in.atEnd() ))// && ((cur_time <= (int) 30*60*freq ) ||  AlignmentContinue ))	
	{
		// этап выставки
		if ((cur_time <= t_alignment) && AlignmentContinue )
		{
			#if 1
			ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, Ab, Omb); //чтение из файла
			#endif
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
			for(int i=0; i<3; ++i)
			{
				// применяем метод Уэлфорда
				MeanAb[i] = (Ldoub) MeanAb[i] + (Ab[i] - MeanAb[i]) / (cur_time + 1); //(cur_time * MeanAb[i] + Ab[i])/(cur_time + 1);
				StdAb[i] = (Ldoub) (1 - 1/(cur_time + 1))*StdAb[i] + (Ab[i] - MeanAb[i])*(Ab[i] - MeanAb[i])/(cur_time + 1);
				MeanOmb[i] = (Ldoub) MeanOmb[i] + (Omb[i] - MeanOmb[i]) / (cur_time + 1);//(cur_time * MeanOmb[i] + Omb[i])/(cur_time + 1);
				StdOmb[i] = (Ldoub) (1 - 1/(cur_time + 1))*StdOmb[i] + (Omb[i] - MeanOmb[i])*(Omb[i] - MeanOmb[i])/(cur_time + 1);
			}
			//Вычисление (ориентации) матрицы перехода Cbn = [c00, c01, c02, c10, c11, c12, c20, c21, c22]
			// Ищем обратную (транспонированную) матрицу
			for(int i=0; i<3; ++i)
			{
				Cbn[index_3(3, 2, i)] = (Ldoub) MeanAb[i] / g;
				Cbn[index_3(3, 1, i)] = (Ldoub) (MeanOmb[i]/ U - MeanAb[i]/g*sin(phi0))/cos(phi0);
			}
			// по алгебраическому дополнению
			Ldoub res1=0, err1=0, res2=0, err2=0;
			TwoProduct(Cbn[index_3(3, 1, 1)], Cbn[index_3(3, 2, 2)], res1, err1);
			TwoProduct(Cbn[index_3(3, 2, 1)], Cbn[index_3(3, 1, 2)], res2, err2);
			Cbn[index_3(3, 0, 0)] = (Ldoub) res1 + err1 + res2 + err2;
			res1=0; err1=0; res2=0; err2=0;
			TwoProduct( Cbn[index_3(3, 1, 2)], Cbn[index_3(3, 2, 0)], res1, err1);
			TwoProduct(Cbn[index_3(3, 1, 0)], Cbn[index_3(3, 2, 2)], res2, err2);
			Cbn[index_3(3, 0, 1)] = (Ldoub) res1 + err1 - (res2 + err2);
			res1=0; err1=0; res2=0; err2=0;
			TwoProduct(Cbn[index_3(3, 1, 0)], Cbn[index_3(3, 2, 1)], res1, err1);
			TwoProduct(Cbn[index_3(3, 1, 1)], Cbn[index_3(3, 2, 0)], res2, err2);
			Cbn[index_3(3, 0, 2)] = (Ldoub) res1 + err1 + res2 + err2;
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

		for (int in_iter=0; in_iter < 4; ++in_iter) // 4 такта, нумерация с нуля, поэтому равентсов нестрогое
		{
			#if 1
			ReadFile(in,  AllowBiasAcc, AllowBiasGyr, AllowRandAcc, AllowRandGyr, Ab, Omb); //чтение из файла
			#endif
			//GeneratedSens(Ab, Omb, Vabs, H0, cur_time, t_alignment, U, g, Cnb);
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
			Ldoub al_w[3] = {0};
			
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
#if 0
				printf("k1[%d] = %.10f\n", ii, k1[ii]);
				printf("k2[%d] = %.10f\n", ii, k2[ii]);
				printf("k3[%d] = %.10f\n", ii, k3[ii]);
				printf("k4[%d] = %.10f\n", ii, k4[ii]);
				printf("Wp[%d] = %.10f", ii, Wp[ii]);
#endif
			}
		}		
		/*Далее идет 100 Гц такт*/
		//Решение задачи ориентации
		Ldoub Thet4[3] = {0}; //Вектор Эйлера
		for (int mmm=0; mmm<3; ++mmm)
			for (int kkk=0; kkk<4; ++kkk)
		{
			Thet4[mmm] += alpha[index_3(4, mmm, kkk)];
		}
		//Для вектора Эйлера необходимо векторное умножение
		Ldoub al_1_2[3] = {0};
		Ldoub al_3_4[3] = {0};
		for (int iii=0; iii<3; ++iii)
			for (int jjj=0; jjj<2; ++jjj)
			{
				al_1_2[iii] += alpha[index_3(4, iii, jjj)];
				al_3_4[iii] += alpha[index_3(4, iii, 3-jjj)];
			}
		Ldoub al_1_2_eig[9] = {0, -al_1_2[2], al_1_2[1], al_1_2[2], 0, -al_1_2[0], -al_1_2[1], al_1_2[0], 0};
		Ldoub temp_res[3] = {0};
		MulMatrD(al_1_2_eig, al_3_4, temp_res, 3, 3, 1);
		
		for (int mmm=0; mmm<3; ++mmm)
		{
			Thet4[mmm] += 2./3*temp_res[mmm]; // после этого вектор Эйлера 
		}
		Ldoub EigThet[9] = {0, -Thet4[2], Thet4[1], Thet4[2], 0, -Thet4[0], -Thet4[1], Thet4[0], 0}; // кососиметрическая матрица верктора Эйлера
		Ldoub EigThet2[9] = {0}; //квадрат кососиметрической матрицы вектора Эйлера
		MulMatrD(EigThet, EigThet, EigThet2,3,3,3);
		Ldoub absThet2 = pow(Thet4[0],2) + pow(Thet4[1],2) + pow(Thet4[2],2); // Квадрат модуля вектора Эйлера
		//absThet2 = 3.970459053e-13;
		
		Ldoub dCbb[9] = {0}; //матрица перехода их связанной в связанную за 1 такт (4*h4)
		for (int iii=0; iii<3; ++iii)
			for(int jjj=0; jjj<3; ++jjj)
			{
				if (iii==jjj)
					dCbb[index_3(3,iii,jjj)] = 1 - (1 - absThet2/6.)*EigThet[index_3(3,iii,jjj)] + (0.5 - absThet2/24.)*EigThet2[index_3(3,iii,jjj)];
				else
					dCbb[index_3(3,iii,jjj)] = 0 - (1 - absThet2/6.)*EigThet[index_3(3,iii,jjj)] + (0.5 - absThet2/24.)*EigThet2[index_3(3,iii,jjj)];
			}
		Ldoub tempCib[9] = {0};
		MulMatrD(dCbb, Cib, tempCib, 3,3,3);
		//переприсваивание Cib = tempCib идет ниже, вместе с Cin
		
		//Вычисление переносных, относительных и абсолютных угловых скорокстей опопрной системы координат
		Omo[0] = (Ldoub) -V[1]/(Rphi + Coordinates[2]);
		Omo[1] = (Ldoub) V[0]/(Rlambda + Coordinates[2]);
		Omo[2] = (Ldoub) V[0]/(Rlambda + Coordinates[2])*tan(phi0); // phi0
		Ldoub omo[3] = {(Ldoub) Omo[0], (Ldoub) Omo[1] + (Ldoub) U*cos(Coordinates[0]), (Ldoub) Omo[2] + (Ldoub) U*sin(Coordinates[0])};// phi0
		Coordinates[0] += (Ldoub) (V[1]/(Rphi + Coordinates[2]))/freq;
		Coordinates[1] += (Ldoub) (V[0]/((Rlambda + Coordinates[2])*cos(Coordinates[0])))/freq;
		Coordinates[2] += (Ldoub) (V[2])/freq;
		//Ldoub EigWb[9] = {0, -Omb[2], Omb[1], Omb[2], 0, -Omb[0], -Omb[1], Omb[0], 0};
		Ldoub EigWo[9] = {0, -omo[2], omo[1], omo[2], 0, -omo[0], -omo[1], omo[0], 0};
		//Ldoub EigWo[9] = {0, -Omo[2], Omo[1], Omo[2], 0, -Omo[0], -Omo[1], Omo[0], 0};
		Ldoub EigWo2[9]={0};//квадрат кососиметрической матрицы абсолютной голвой скорости опорной с.к
		
		MulMatrD(EigWo, EigWo, EigWo2, 3, 3, 3);

		Ldoub tempCin[9] = {0};
		Ldoub dCnn[9] = {0};
		for (int iii=0; iii<3; ++iii)
			for(int jjj=0; jjj<3; ++jjj)
			{
				if (iii==jjj)
					dCnn[index_3(3,iii,jjj)] = 1 - h*EigWo[index_3(3,iii,jjj)] + pow(h,2)*EigWo2[index_3(3,iii,jjj)]/2. ;
				else
					dCnn[index_3(3,iii,jjj)] = 0 - h*EigWo[index_3(3,iii,jjj)] + pow(h,2)*EigWo2[index_3(3,iii,jjj)]/2. ;
			}
		//Cin = dCnn * Cin	
		MulMatrD(dCnn, Cin, tempCin, 3, 3, 3);
		for (int ii=0; ii<9; ii++) 
		{
			Cib[ii] = tempCib[ii];
			Cin[ii] = tempCin[ii];
		}
#if 0
		//Шаманим с матрицей body
		MatrOB(Thet4[2], Thet4[1], Thet4[0], Cib, 3);
#endif
#if 0
		//Шаманим с матрицей o, она же n
		MatrOB(Omo[2], Omo[1], Omo[0], Cin, 3); // если Omo, то ошибки по углам на уровне 1е-11
#endif
		Ldoub Cbi[9] = {0};
		Transpose2M(Cib, Cbi, 3);//транспонированная матрица Cbi
		
		// Решение уравнения Пуассона
		MulMatrD(Cin, Cbi, Cbn, 3,3,3);
#if 0
		for (int iii=0; iii<3; ++iii)
			for(int jjj=0; jjj<3; ++jjj)
				if (iii==jjj) Cbn[index_3(3, iii,jjj)] = 1;
				else Cbn[index_3(3, iii,jjj)] = 0;
#endif
		/*Процедура нормирования и ортогонализации*/
#if 0
		//контроль масштаба
		//Строки
		for (int iii=0; iii<3; ++iii)
			{
				Ldoub string[3] = {Cbn[index_3(3,iii,0)], Cbn[index_3(3,iii,1)], Cbn[index_3(3,iii,2)]};
				Ldoub norm[1];
				MulMatrD(string, string, norm, 1,3,1);
				norm[0] = 1 - norm[0];
				for (int jjj=0; jjj < 3; ++jjj)
					Cbn[index_3(3,iii,jjj)] = Cbn[index_3(3,iii,jjj)]  - 0.5*norm[0]*Cbn[index_3(3,iii,jjj)];
			}
		//Столбцы
		for (int iii=0; iii<3; ++iii)
			{
				Ldoub string[3] = {Cbn[index_3(3,0,iii)], Cbn[index_3(3,1,iii)], Cbn[index_3(3,2,iii)]};
				Ldoub norm[1];
				MulMatrD(string, string, norm, 1,3,1);
				norm[0] = 1 - norm[0];
				for (int jjj=0; jjj < 3; ++jjj)
					Cbn[index_3(3,jjj,iii)] = Cbn[index_3(3,jjj,iii)]  - 0.5*norm[0]*Cbn[index_3(3,jjj,iii)];
			}
#endif
#if 0
		// ортогонализация
		//по строкам
		for(int iii=0; iii < 3; ++iii)
			for(int jjj=0; jjj < 3; ++jjj)
			{
				if(iii==jjj) continue;
				Ldoub string[3] = {Cbn[index_3(3,iii,0)], Cbn[index_3(3,iii,1)], Cbn[index_3(3,iii,2)]};
				Ldoub column[3] = {Cbn[index_3(3,jjj,0)], Cbn[index_3(3,jjj,1)], Cbn[index_3(3,jjj,2)]};
				Ldoub norm[1];
				MulMatrD(string, column,norm, 1,3,1);
				for(int kkk=0; kkk<2; ++kkk)
				{
					Cbn[index_3(3,iii,kkk)] = Cbn[index_3(3,iii,kkk)] - 0.5*norm[0]*Cbn[index_3(3,jjj,kkk )];
					Cbn[index_3(3,jjj,kkk)] = Cbn[index_3(3,jjj,kkk)] - 0.5*norm[0]*Cbn[index_3(3,iii,kkk)];
				}

			}
		//по столбцам
		for(int iii=0; iii < 3; ++iii)
			for(int jjj=0; jjj < 3; ++jjj)
			{
				if(iii==jjj) continue;
				Ldoub string[3] = {Cbn[index_3(3,0,iii)], Cbn[index_3(3,1,iii)], Cbn[index_3(3,2,iii)]};
				Ldoub column[3] = {Cbn[index_3(3,0,jjj)], Cbn[index_3(3,1,jjj)], Cbn[index_3(3,2,jjj)]};
				Ldoub norm[1];
				MulMatrD(string, column,norm, 1,3,1);
				for(int kkk=0; kkk<2; ++kkk)
				{
					Cbn[index_3(3,kkk,iii)] = Cbn[index_3(3,kkk,iii)] - 0.5*norm[0]*Cbn[index_3(3,kkk,jjj )];
					Cbn[index_3(3,kkk,jjj)] = Cbn[index_3(3,kkk,jjj)] - 0.5*norm[0]*Cbn[index_3(3,kkk,iii)];
				}

			}
#endif
		
		// вычисление углов ориентации через МНК (как в выставке)
		Orientation[0] = (Ldoub) atan2(Cbn[index_3(3, 0, 1)],Cbn[index_3(3, 1, 1)]);
		Orientation[1] =  (Ldoub) - atan2(Cbn[index_3(3, 2, 0)],Cbn[index_3(3, 2, 2)]);
		Ldoub c0 = (Ldoub) sqrt(Cbn[index_3(3, 2, 0)]* Cbn[index_3(3, 2, 0)] + Cbn[index_3(3, 2, 2)]*Cbn[index_3(3, 2, 2)]);
		Orientation[2] = (Ldoub) atan2(Cbn[index_3(3, 2, 1)],c0);
		/*Проверка ре*/
		
		
		/*Решение задачи навигации*/
		MulMatrD(Cbn, Wp, Ao, 3, 3, 1); // перепроектирование из связаных осей в навигационные. Здесь Ao -- уже не ускорения, а приращшение скоросетй
		
		//Кориолисовы добавки
		Ldoub aCoriolis[3] = {0};
		aCoriolis[0] = (Ldoub) ((Ldoub) omo[1]*V[2] - (Ldoub) omo[2]*V[1] + (Ldoub) U*cos(Coordinates[0])*V[2] - (Ldoub) U*sin(Coordinates[0])*V[1]);
		aCoriolis[1] = (Ldoub) ((Ldoub) -omo[0]*V[2] + (Ldoub) omo[2]*V[0] + (Ldoub) U*sin(Coordinates[0])*V[0]);
		aCoriolis[2] = (Ldoub) ((Ldoub) omo[0]*V[1] - (Ldoub) omo[1]*V[0] - (Ldoub) U*cos(Coordinates[0])*V[0]);
		
		V[0] = V[0] + Ao[0] - 1*aCoriolis[0]*h; // Ve 
		V[1] = V[1] + Ao[1] - 1*aCoriolis[1]*h; //Vn 
		
		//Ошибки по скоростям
		for (int iii=0; iii<2; ++iii)
			Err_V[iii] = V[iii] - V0[iii];

		//Ошибки по координатам в м
		CoordError[0] += (V[0] - (Ldoub) Vabs*sin(H0)) * h;
		CoordError[1] += (V[1] - (Ldoub) Vabs*cos(H0)) * h;
	
		Rlambda = (Ldoub) R/sqrt(1.-pow(e,2)*pow(sin(Coordinates[0]),2) );
		Rphi = (Ldoub) R*(1. - pow(e,2))/(sqrt(1.-pow(e,2)*pow(sin(Coordinates[0]),2) ) * (1.-pow(e,2)*pow(sin(Coordinates[0]),2)  ) );
		
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
	//getc(stdin);
	return 0;
}
