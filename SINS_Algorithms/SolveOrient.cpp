#include "SolveOrient.h"
#include "std.h"
#include "matrix.h"
#include <cmath>
#include "quaternions.h"
//#define dev
#ifdef dev
#include <stdio.h>
#endif

void SolveOrient(Ldoub* Omb, quaternion* Qf, Ldoub* Cbn, Ldoub* Orientation, Ldoub* Coordinates, Ldoub* Omo, Ldoub* omo, Ldoub* V, Ldoub phi0, Ldoub Rphi, Ldoub Rlambda, int freq, Ldoub h, const Ldoub U, int iter, bool& sw, Ldoub k2, Ldoub* DVerr, bool allowCorr)
{
#if 0
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
	
	//Вычисление переносных, относительных и абсолютных угловых скоростей опопрной системы координат
	Omo[0] = (Ldoub) -V[1]/(Rphi + Coordinates[2]);
	Omo[1] = (Ldoub) V[0]/(Rlambda + Coordinates[2]);
	Omo[2] = (Ldoub) V[0]/(Rlambda + Coordinates[2])*tan(Coordinates[0]); // tan(phi0)

	omo[0] = (Ldoub) Omo[0] - allowCorr * k2/Rphi * DVerr[1]; 
    omo[1] = (Ldoub) Omo[1] + (Ldoub) U*cos(Coordinates[0]) + allowCorr * k2/Rphi * DVerr[0];
    omo[2] = Omo[2] + (Ldoub) U*sin(Coordinates[0]);
	Coordinates[0] += (Ldoub) (V[1]/(Rphi + Coordinates[2]))/freq;
	Coordinates[1] += (Ldoub) (V[0]/((Rlambda + Coordinates[2])*cos(Coordinates[0])))/freq;
	//Coordinates[2] += (Ldoub) (V[2])/freq;

	
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


#ifdef dev //dev --- develop, разработка, отладка
		
	printf("%d:\t столбцы:\n", iter);
	//Проверяем выполняются ли свойства матрицы направляющих косинусов
	//Скалярное произведение столбцов
	Ldoub Prodij (0.); //результат произведения столбцов МНК
	for (int kkk=0; kkk<3; kkk++)
	{
		Ldoub ci[3] = {Cib[index_3(3, kkk, 0)], Cib[index_3(3, kkk, 1)], Cib[index_3(3, kkk, 2)]};
		for (int iii=0; iii<3; iii++)
		{
			Ldoub cj[3];
			for (int jjj=0; jjj<3; jjj++)
			{
				cj[jjj] = Cib[index_3(3, iii, jjj)];
			}
			MulMatrD(ci, cj, &Prodij, 1,3,1);
			printf("%d * %d = %.8f\n",kkk, iii, Prodij);
		}
	}
	//Скалярное произведение строк	
	printf("%d:\t строки:\n", iter);
	Prodij=0.; //результат произведения строк МНК
	for (int kkk=0; kkk<3; kkk++)
	{
		Ldoub ci[3] = {Cib[index_3(3, 0, kkk)], Cib[index_3(3, 1, kkk)], Cib[index_3(3, 2, kkk)]};
		for (int iii=0; iii<3; iii++)
		{
			Ldoub cj[3];
			for (int jjj=0; jjj<3; jjj++)
			{
				cj[jjj] = Cib[index_3(3, jjj, iii)];
			}
			MulMatrD(ci, cj, &Prodij, 1,3,1);
			printf("%d * %d = %.8f\n",kkk, iii, Prodij);
		}
	}
#endif


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
	MulMatrD(Cin, Cbi, Cbn, 3, 3, 3);
#if 0
	for (int iii=0; iii<3; ++iii)
		for(int jjj=0; jjj<3; ++jjj)
			if (iii==jjj) Cbn[index_3(3, iii,jjj)] = 1;
			else Cbn[index_3(3, iii,jjj)] = 0;
#endif

#if 0 // коррекция матрицы перехода
	Normalization(Cbn, sw);
	Ortogonalization(Cbn, sw);
	sw = !sw; //меняем направление нормализации/ортогонализации
#endif
		
	// вычисление углов ориентации через МНК (как в выставке)
	Orientation[0] = (Ldoub) atan2(Cbn[index_3(3, 0, 1)],Cbn[index_3(3, 1, 1)]);
	Orientation[1] =  (Ldoub) - atan2(Cbn[index_3(3, 2, 0)],Cbn[index_3(3, 2, 2)]);
	Ldoub c0 = (Ldoub) sqrt(Cbn[index_3(3, 2, 0)]* Cbn[index_3(3, 2, 0)] + Cbn[index_3(3, 2, 2)]*Cbn[index_3(3, 2, 2)]);
	Orientation[2] = (Ldoub) atan2(Cbn[index_3(3, 2, 1)],c0);
	/*Проверка ре*/
#endif
	//Вычисление переносных, относительных и абсолютных угловых скоростей опопрной системы координат
	Omo[0] = (Ldoub) -V[1]/(Rphi + Coordinates[2]);
	Omo[1] = (Ldoub) V[0]/(Rlambda + Coordinates[2]);
	Omo[2] = (Ldoub) V[0]/(Rlambda + Coordinates[2])*tan(Coordinates[0]); // tan(phi0)

	omo[0] = (Ldoub) Omo[0] - allowCorr * k2/Rphi * DVerr[1]; 
    omo[1] = (Ldoub) Omo[1] + (Ldoub) U*cos(Coordinates[0]) + allowCorr * k2/Rphi * DVerr[0];
    omo[2] = Omo[2] + (Ldoub) U*sin(Coordinates[0]);
	Coordinates[0] += (Ldoub) (V[1]/(Rphi + Coordinates[2]))/freq;
	Coordinates[1] += (Ldoub) (V[0]/((Rlambda + Coordinates[2])*cos(Coordinates[0])))/freq;

	Ldoub F;
	MulMatrD(Omb, Omb, &F, 1, 3, 1); //сумма квадратов элементов Omb
	F = sqrt(F);
	quaternion DL(cos(F*h/2.), Omb[0]/F*sin(F*h/2.), Omb[1]/F*sin(F*h/2.), Omb[2]/F*sin(F*h/2.));

	quaternion Qp = (*Qf)*DL; //предварительный (быстрый) кватернион

	Ldoub Om0;
	MulMatrD(omo, omo, &Om0, 1, 3, 1); //сумма квадратов элементов Omb
	Om0 = sqrt(Om0);

	quaternion Dm(cos(Om0*h/2.), -omo[0]/Om0*sin(Om0*h/2.), -omo[1]/Om0*sin(Om0*h/2.), -omo[2]/Om0*sin(Om0*h/2.));
	quaternion tempQuat = Dm*(Qp);
	
	//Костыль выглядит так, вообще, хотелось сделать Qf = Dm * (Qp)
	Qf->w = tempQuat.w;
	Qf->x = tempQuat.x;
	Qf->y = tempQuat.y;
	Qf->z = tempQuat.z;

	Qf->normalize(1e-5);

	Quat2Matr(Qf, Cbn); //пересчет кватерниона в матрицу

	// вычисление углов ориентации через МНК (как в выставке)
	Orientation[0] = (Ldoub) atan2(Cbn[index_3(3, 0, 1)],Cbn[index_3(3, 1, 1)]);
	Orientation[1] =  (Ldoub) - atan2(Cbn[index_3(3, 2, 0)],Cbn[index_3(3, 2, 2)]);
	Ldoub c0 = (Ldoub) sqrt(Cbn[index_3(3, 2, 0)]* Cbn[index_3(3, 2, 0)] + Cbn[index_3(3, 2, 2)]*Cbn[index_3(3, 2, 2)]);
	Orientation[2] = (Ldoub) atan2(Cbn[index_3(3, 2, 1)],c0);
}