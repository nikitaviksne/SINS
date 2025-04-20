#include "SolveNav.h"
#include "std.h"
#include "matrix.h"
#include <cmath>

void SolveNav(Ldoub* Wp, Ldoub* Cbn, Ldoub* Ao, Ldoub* V, Ldoub* Coordinates, Ldoub* CoordError, Ldoub* Err_V, Ldoub* omo, Ldoub h, Ldoub& Rphi, Ldoub& Rlambda, const Ldoub U, Ldoub R, Ldoub e, Ldoub H0, Ldoub* V0 )
{
    MulMatrD(Cbn, Wp, Ao, 3, 3, 1); // перепроектирование из связаных осей в навигационные. Здесь Ao -- уже не ускорения, а приращшение скоросетй
		
	//Кориолисовы добавки
	Ldoub aCoriolis[3] = {0};
	aCoriolis[0] = (Ldoub) ((Ldoub) omo[1]*V[2] - (Ldoub) omo[2]*V[1] + (Ldoub) U*cos(Coordinates[0])*V[2] - (Ldoub) U*sin(Coordinates[0])*V[1]);
	aCoriolis[1] = (Ldoub) ((Ldoub) -omo[0]*V[2] + (Ldoub) omo[2]*V[0] + (Ldoub) U*sin(Coordinates[0])*V[0]);
	aCoriolis[2] = (Ldoub) ((Ldoub) omo[0]*V[1] - (Ldoub) omo[1]*V[0] - (Ldoub) U*cos(Coordinates[0])*V[0]);
	
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