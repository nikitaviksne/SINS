#include "alignment.h"
#include "std.h"
#include "mathematics.h"
#include <cmath>
#include "quaternions.h"

void alignment(/*входные*/Ldoub *Ab, Ldoub *Omb, Ldoub *MeanAb, Ldoub *MeanOmb, Ldoub *StdAb, Ldoub *StdOmb, int iter, const Ldoub g, const Ldoub U, Ldoub phi0, /*выходные*/Ldoub *Cbn, quaternion *Qf)
{	/*процедура выставки (Б)ИНС*/
	for(int i=0; i<3; ++i)
	{
		// применяем метод Уэлфорда
			MeanAb[i] = (Ldoub) (iter * MeanAb[i] + Ab[i])/(iter + 1.); //MeanAb[i] + (Ab[i] - MeanAb[i]) / (iter + 1);
			StdAb[i] = (Ldoub) (1 - 1./(iter + 1))*StdAb[i] + (Ab[i] - MeanAb[i])*(Ab[i] - MeanAb[i])/(iter + 1.);
			MeanOmb[i] = (Ldoub) (iter * MeanOmb[i] + Omb[i])/(iter + 1.);//MeanOmb[i] + (Omb[i] - MeanOmb[i]) / (iter + 1);
			StdOmb[i] = (Ldoub) (1 - 1./(iter + 1))*StdOmb[i] + (Omb[i] - MeanOmb[i])*(Omb[i] - MeanOmb[i])/(iter + 1.);
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
	Matr2Quat(Cbn, Qf);
}
