#ifndef QUATERNIONSH
#define QUATERNIONSH

#include "std.h"

class quaternion
{
	public:
		Ldoub w;
		Ldoub x;
		Ldoub y;
		Ldoub z;
		
		quaternion(Ldoub a, Ldoub b, Ldoub c, Ldoub d);
        quaternion();
        #if 1
		
		quaternion& operator / (const Ldoub& val);  //перегрузка оператора деления на число

		quaternion& operator * (const Ldoub& val ); //перегрузка оператора умножения на число
        
		quaternion operator * (const quaternion& quat ) const; //перегрузка оператора умножения кватернионов
		
		quaternion operator = (const quaternion& quat); //перегрузка оператора присваивания
		//quaternion& operator = (const quaternion& quat); //перегрузка оператора присваивания
		#endif
		
		void print();
		
		Ldoub norma();

		void normalize(Ldoub limit = 1e-5);
};


void MulQuat();

void Quat2Matr(quaternion Q, Ldoub* C);
void Quat2Matr(quaternion *Q, Ldoub* C);

void Matr2Quat(Ldoub* C, quaternion* Q);

#endif