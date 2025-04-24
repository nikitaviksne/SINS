#include "SolveNav.h"
#include "std.h"
#include "matrix.h"
#include <cmath>

//#define dev

#ifdef dev
#include <stdio.h>
#endif

void SolveNav(Ldoub* Wp, Ldoub* Ab, Ldoub* Cbn, Ldoub* Ao, Ldoub* V, Ldoub* Coordinates, Ldoub* CoordError, Ldoub* Err_V, Ldoub* omo, Ldoub h, Ldoub& Rphi, Ldoub& Rlambda, const Ldoub U, Ldoub R, Ldoub e, Ldoub H0, Ldoub* V0 )
{
    MulMatrD(Cbn, Wp, Ao, 3, 3, 1); // перепроектирование из связаных осей в навигационные. Здесь Ao -- уже не ускорения, а приращение скоростей
	
	#ifdef dev
	Ldoub detC (0);
	Determinant_3(Cbn, detC);
	printf("Определитель = %.8f\n", detC);
	Ldoub AccO[3] = {0};// Ускорения в связанной с.к
	MulMatrD(Cbn, Ab, AccO, 3, 3, 1); // перепроектирование Ускорений из связаных осей в навигационные.
	printf("Ускорения в опорной с.к AccO[0] = %.8f AccO[1] = %.8f AccO[2] = %.8f\n", AccO[0], AccO[1], AccO[2]);
	printf("Прииращения скоростей в связанной с.к Wp[0] = %.8f Wp[1] = %.8f Wp[2] = %.8f\n", Wp[0], Wp[1], Wp[2]);
	printf("Прииращения скоростей в опорной с.к Wo[0] = %.8f Wo[1] = %.8f Wo[2] = %.8f\n", Ao[0], Ao[1], Ao[2]);
	#endif //dev
	
	//Кориолисовы добавки
	Ldoub aCoriolis[3] = {0};
	aCoriolis[0] = (Ldoub) ((Ldoub) omo[1]*V[2] - (Ldoub) omo[2]*V[1] + (Ldoub) U*cos(Coordinates[0])*V[2] - (Ldoub) U*sin(Coordinates[0])*V[1]);
	aCoriolis[1] = (Ldoub) ((Ldoub) -omo[0]*V[2] + (Ldoub) omo[2]*V[0] + (Ldoub) U*sin(Coordinates[0])*V[0]);
	aCoriolis[2] = (Ldoub) ((Ldoub) omo[0]*V[1] - (Ldoub) omo[1]*V[0] - (Ldoub) U*cos(Coordinates[0])*V[0]);
	
	#ifdef dev
	printf("aCoriolis[0] = %.8f aCoriolis[1] = %.8f aCoriolis[2] = %.8f\n", aCoriolis[0], aCoriolis[1], aCoriolis[2]);
	#endif //dev
	V[0] = V[0] + Ao[0] - 1*aCoriolis[0]*h; // Ve 
	V[1] = V[1] + Ao[1] - 1*aCoriolis[1]*h; //Vn 
	
	//Ошибки по скоростям в м/с и координатам в м
	for (int iii=0; iii<2; ++iii)
    {
		Err_V[iii] = V[iii] - V0[iii];
        CoordError[iii] += (Err_V[iii]) * h;
    }


	Rlambda = R/sqrt(1.-pow(e,2)*pow(sin(Coordinates[0]),2) );
	Rphi = R*(1. - pow(e,2))/(sqrt(1.-pow(e,2)*pow(sin(Coordinates[0]),2) ) * (1.-pow(e,2)*pow(sin(Coordinates[0]),2)  ) );
		
}