#ifndef QUATERNIONS_H
#define QUATERNIONS_H

#include "std.h"

class quaternion
{
	public:
		Ldoub w;
		Ldoub x;
		Ldoub y;
		Ldoub z;
		
		quaternion(Ldoub a, Ldoub b, Ldoub c, Ldoub d);
        #if 1
		
		quaternion& operator / (const Ldoub& val);  //перегрузка оператора деления на число

		quaternion operator * (const Ldoub& val ); //перегрузка оператора умножения на число
        
		quaternion operator * (const quaternion& quat ) const; //перегрузка оператора умножения кватернионов
		
		quaternion operator = (const quaternion& quat); //перегрузка оператора присваивания
		#endif
		
		void print();
		
		Ldoub norma();
};

void MulQuat();

void Quat2Matr();

void Matr2Quat();

#endif //QUATERNIONS_H