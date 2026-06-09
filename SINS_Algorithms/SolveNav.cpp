#include "SolveNav.h"
#include "std.h"
#include "ap.h"
#include "matrix.h"
#include <cmath>

//#define dev

#ifdef dev
#include <stdio.h>
#endif

void SolveNav(Ldoub* Wp, Ldoub* Ab, Ldoub* Cbn, Ldoub* Wo, Ldoub* Ao, Ldoub* V, Ldoub* Coordinates, Ldoub* CoordError, Ldoub* Err_V, Ldoub* Omo/*Переносная скорость */, Ldoub h, Ldoub& Rphi, Ldoub& Rlambda, const Ldoub U, Ldoub R, Ldoub e, Ldoub H0, Ldoub* V0, Ldoub k1, alglib::real_1d_array ErrVins, bool allowCorr)
{
    MulMatrD(Cbn, Wp, Wo, 3, 3, 1); // перепроектирование из связаных осей в навигационные. Здесь Wo -- уже не ускорения, а приращение скоростей

	MulMatrD(Cbn, Ab, Ao, 3, 3, 1); // перепроектирование из связаных осей в навигационные. Ao ускорения в опорной с.к
	
	#ifdef dev
	Ldoub detC (0);
	Determinant_3(Cbn, detC);
	printf("Определитель = %.8f\n", detC);
	Ldoub AccO[3] = {0};// Ускорения в связанной с.к
	MulMatrD(Cbn, Ab, AccO, 3, 3, 1); // перепроектирование Ускорений из связаных осей в навигационные.
	printf("Ускорения в опорной с.к AccO[0] = %.8f AccO[1] = %.8f AccO[2] = %.8f\n", AccO[0], AccO[1], AccO[2]);
	printf("Прииращения скоростей в связанной с.к Wp[0] = %.8f Wp[1] = %.8f Wp[2] = %.8f\n", Wp[0], Wp[1], Wp[2]);
	printf("Прииращения скоростей в опорной с.к Wo[0] = %.8f Wo[1] = %.8f Wo[2] = %.8f\n", Wo[0], Wo[1], Wo[2]);
	#endif //dev
	
	//Кориолисовы добавки
	Ldoub aCoriolis[3] = {0};
	#if 1
	aCoriolis[0] = (Ldoub) ((Ldoub) Omo[1]*V[2] - (Ldoub) Omo[2]*V[1] + (Ldoub) U*cos(Coordinates[0])*V[2] - (Ldoub) U*sin(Coordinates[0])*V[1]);
	aCoriolis[1] = (Ldoub) ((Ldoub) -Omo[0]*V[2] + (Ldoub) Omo[2]*V[0] + (Ldoub) U*sin(Coordinates[0])*V[0]);
	aCoriolis[2] = (Ldoub) ((Ldoub) Omo[0]*V[1] - (Ldoub) Omo[1]*V[0] - (Ldoub) U*cos(Coordinates[0])*V[0]);
	
	/*// по Салычеву стр 70 applied navigation algorithm strapdown system navigation alogithm
	aCoriolis[0] = (Ldoub) ((Ldoub) Omo[1]*V[2] - (Ldoub) 2*U*sin(Coordinates[0])*V[1] + (Ldoub) U*cos(Coordinates[0])*V[2] - (Ldoub) U*sin(Coordinates[0])*V[1]);
	aCoriolis[1] = (Ldoub) ((Ldoub) -Omo[0]*V[2] + (Ldoub) 2*U*sin(Coordinates[0])*V[0] );
	aCoriolis[2] = (Ldoub) ((Ldoub) Omo[0]*V[1] - (Ldoub) Omo[1]*V[0] - (Ldoub) U*cos(Coordinates[0])*V[0]);*/
	#else // по соображениям теоретической механики. Ускорение Кориолиса равно удвоенному векторному произведению абсолютной угловой скорости подвижного базиса на относительнуюлинейную скорость
	aCoriolis[0] = (Ldoub) 2 * ((Ldoub) Omo[1]*V[2] - (Ldoub) Omo[2]*V[1]);
	aCoriolis[1] = (Ldoub) 2 * ((Ldoub) -Omo[0]*V[2] + (Ldoub) Omo[2]*V[0]);
	aCoriolis[2] = (Ldoub) 2 * ((Ldoub) Omo[0]*V[1] - (Ldoub) Omo[1]*V[0]);
	#endif
	
	#ifdef dev
	printf("aCoriolis[0] = %.8f aCoriolis[1] = %.8f aCoriolis[2] = %.8f\n", aCoriolis[0], aCoriolis[1], aCoriolis[2]);
	#endif //dev
	
	for (int iii=0; iii < 3; ++iii)	
	{
		V[iii] = V[iii] + Wo[iii] - (1*aCoriolis[iii] )*h; // Ve 
	}

	if (allowCorr)
	{
		for (int i =0; i< 3; ++i)
		{
			V[i] -= ErrVins[3+i]; //+ allowCorr * k1 * ErrVins[iii]*h
		}
	}

	const Ldoub g (9.81);
	// V[2] = V[2] + Wo[2] - (1*aCoriolis[2] + g)*h;
	
	
	Coordinates[0] += (Ldoub) (V[1]/(Rphi + Coordinates[2]))*h;
	Coordinates[1] += (Ldoub) (V[0]/((Rlambda + Coordinates[2])*cos(Coordinates[0])))*h;
	Coordinates[2] += (Ldoub) (V[2]) * h;

	//Ошибки по скоростям в м/с и координатам в м
	for (int iii=0; iii<3; ++iii)
    {
		Err_V[iii] = V[iii] - V0[iii];
        CoordError[iii] += (Err_V[iii]) * h;
    }


	Rlambda = R/sqrt(1.-pow(e,2)*pow(sin(Coordinates[0]),2) );
	Rphi = R*(1. - pow(e,2))/(sqrt(1. - pow(e,2)*pow(sin(Coordinates[0]),2) ) * (1.-pow(e,2)*pow(sin(Coordinates[0]),2)  ) );
		
}