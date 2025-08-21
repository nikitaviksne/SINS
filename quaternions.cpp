#include <cmath>
#include "include/std.h"
#include "quaternions.h"
#include <stdio.h>

quaternion::quaternion(Ldoub a, Ldoub b, Ldoub c, Ldoub d)
{
	this->w = a;
	this->x = b;
	this->y = c;
	this->z = d;
}

quaternion::quaternion()
{
	this->w = 0.;
	this->x = 0.;
	this->y = 0.;
	this->z = 0.;
}

quaternion& quaternion::operator / (const Ldoub& val)  //перегрузка оператора деления на число
{
	this->w /= val;
	this->x /= val;
	this->y /= val;
	this->z /= val;
	return *this;
}			

quaternion& quaternion::operator * (const Ldoub& val ) //перегрузка оператора умножения на число
{
	this->w *= val;
	this->x *= val;
	this->y *= val;
	this->z *= val;
	return *this;
}

quaternion quaternion::operator * (const quaternion& quat ) const //перегрузка оператора умножения кватернионов
{
	Ldoub w, x, y, z;
	w = this->w * quat.w - this->x * quat.x - this->y * quat.y - this->z * quat.z;
	x = this->w * quat.x + this->x * quat.w + this->y * quat.z - this->z * quat.y;
	y = this->w * quat.y - this->x * quat.z + this->y * quat.w + this->z * quat.x;
	z = this->w * quat.z + this->x * quat.y - this->y * quat.x + this->z * quat.w; 
	quaternion tempQuat = quaternion(w,x,y,z);
	return quaternion(w,x,y,z);
}

quaternion quaternion::operator = (const quaternion& quat) //перегрузка оператора присваивания
{
	this->w = quat.w;
	this->x = quat.x;
	this->y = quat.y;
	this->z = quat.z;
	return *this;
}

void quaternion::print()
{
	printf("[%.6f, %.6f, %.6f, %.6f]\n", this->w, this->x, this->y, this->z);
}

Ldoub quaternion::norma()
{
	return pow(this->w, 2.) + pow(this->x, 2.) + pow(this->y, 2.) + pow(this->z, 2.);
}
void quaternion::normalize(Ldoub limit)
{
	Ldoub norma (this->norma());
	if ( (1 - norma) > limit)
	{
		this->w *= ((1 + norma/2.));
		this->x *= ((1 + norma/2.));
		this->y *= ((1 + norma/2.));
		this->z *= ((1 + norma/2.));
	}
}
void Quat2Matr(quaternion Q, Ldoub* C)
{
	//из лекций Быковского
#if 0
	C[index_3(3, 0, 0)] = 2*(pow(Q.w,2) + pow(Q.x,2)) - 1;
	C[index_3(3, 0, 1)] = 2*(Q.x * Q.y - Q.w * Q.z);
	C[index_3(3, 0, 2)] = 2*(Q.x * Q.z + Q.w * Q.y);
	C[index_3(3, 1, 0)] = 2*(Q.x * Q.y + Q.w * Q.z);
	C[index_3(3, 1, 1)] = 2*(pow(Q.w,2) + pow(Q.y,2)) - 1;
	C[index_3(3, 1, 2)] = 2*(Q.y * Q.z - Q.w * Q.x);
	C[index_3(3, 2, 0)] = 2*(Q.x * Q.z - Q.w * Q.y);
	C[index_3(3, 2, 1)] = 2*(Q.y * Q.z + Q.w * Q.x);
	C[index_3(3, 2, 2)] = 2*(pow(Q.w,2) + pow(Q.z,2)) - 1;
#endif
#if 0
//из интернета
	C[index_3(3, 0, 0)] = 1 - 2*(pow(Q.y,2) + pow(Q.z,2));
	C[index_3(3, 0, 1)] = 2*(Q.x * Q.y + Q.w * Q.z);
	C[index_3(3, 0, 2)] = 2*(Q.x * Q.z - Q.w * Q.y);
	C[index_3(3, 1, 0)] = 2*(Q.x * Q.y - Q.w * Q.z);
	C[index_3(3, 1, 1)] = 1 - 2*(pow(Q.x,2) + pow(Q.z,2));
	C[index_3(3, 1, 2)] = 2*(Q.y * Q.z + Q.w * Q.x);
	C[index_3(3, 2, 0)] = 2*(Q.x * Q.z + Q.w * Q.y);
	C[index_3(3, 2, 1)] = 2*(Q.y * Q.z - Q.w * Q.x);
	C[index_3(3, 2, 2)] = 1 - 2*(pow(Q.x,2) + pow(Q.y,2));
#endif
#if 1
	C[index_3(3, 0, 0)] = pow(Q.w,2) + pow(Q.x,2) - pow(Q.y,2) - pow(Q.z,2);
    C[index_3(3, 0, 1)] = 2*(Q.x*Q.y - Q.w*Q.z);
    C[index_3(3, 0, 2)] = 2*(Q.x*Q.z + Q.w*Q.y);
    C[index_3(3, 1, 0)] = 2*(Q.x*Q.y + Q.w*Q.z);
    C[index_3(3, 1, 1)] = pow(Q.w,2) + pow(Q.y,2) - pow(Q.x,2) - pow(Q.z,2);
    C[index_3(3, 1, 2)] = 2*(Q.y*Q.z - Q.w*Q.x);
    C[index_3(3, 2, 0)] = 2*(Q.x*Q.z - Q.w*Q.y);
    C[index_3(3, 2, 1)] = 2*(Q.y*Q.z + Q.w*Q.x);
    C[index_3(3, 2, 2)] = pow(Q.w,2) + pow(Q.z,2) - pow(Q.x,2) - pow(Q.y,2);
#endif
};

void Quat2Matr(quaternion *Q, Ldoub* C)
{
	//из лекций Быковского
#if 0
	C[index_3(3, 0, 0)] = 2*(pow(Q->w,2) + pow(Q->x,2)) - 1;
	C[index_3(3, 0, 1)] = 2*(Q->x * Q->y - Q->w * Q->z);
	C[index_3(3, 0, 2)] = 2*(Q->x * Q->z + Q->w * Q->y);
	C[index_3(3, 1, 0)] = 2*(Q->x * Q->y + Q->w * Q->z);
	C[index_3(3, 1, 1)] = 2*(pow(Q->w,2) + pow(Q->y,2)) - 1;
	C[index_3(3, 1, 2)] = 2*(Q->y * Q->z - Q->w * Q->x);
	C[index_3(3, 2, 0)] = 2*(Q->x * Q->z - Q->w * Q->y);
	C[index_3(3, 2, 1)] = 2*(Q->y * Q->z + Q->w * Q->x);
	C[index_3(3, 2, 2)] = 2*(pow(Q->w,2) + pow(Q->z,2)) - 1;
#endif
#if 0
//из интернета
	C[index_3(3, 0, 0)] = 1 - 2*(pow(Q->y,2) + pow(Q->z,2));
	C[index_3(3, 0, 1)] = 2*(Q->x * Q->y + Q->w * Q->z);
	C[index_3(3, 0, 2)] = 2*(Q->x * Q->z - Q->w * Q->y);
	C[index_3(3, 1, 0)] = 2*(Q->x * Q->y - Q->w * Q->z);
	C[index_3(3, 1, 1)] = 1 - 2*(pow(Q->x,2) + pow(Q->z,2));
	C[index_3(3, 1, 2)] = 2*(Q->y * Q->z + Q->w * Q->x);
	C[index_3(3, 2, 0)] = 2*(Q->x * Q->z + Q->w * Q->y);
	C[index_3(3, 2, 1)] = 2*(Q->y * Q->z - Q->w * Q->x);
	C[index_3(3, 2, 2)] = 1 - 2*(pow(Q->x,2) + pow(Q->y,2));
#endif
#if 1
	C[index_3(3, 0, 0)] = pow(Q->w,2) + pow(Q->x,2) - pow(Q->y,2) - pow(Q->z,2);
    C[index_3(3, 0, 1)] = 2*(Q->x*Q->y - Q->w*Q->z);
    C[index_3(3, 0, 2)] = 2*(Q->x*Q->z + Q->w*Q->y);
    C[index_3(3, 1, 0)] = 2*(Q->x*Q->y + Q->w*Q->z);
    C[index_3(3, 1, 1)] = pow(Q->w,2) + pow(Q->y,2) - pow(Q->x,2) - pow(Q->z,2);
    C[index_3(3, 1, 2)] = 2*(Q->y*Q->z - Q->w*Q->x);
    C[index_3(3, 2, 0)] = 2*(Q->x*Q->z - Q->w*Q->y);
    C[index_3(3, 2, 1)] = 2*(Q->y*Q->z + Q->w*Q->x);
    C[index_3(3, 2, 2)] = pow(Q->w,2) + pow(Q->z,2) - pow(Q->x,2) - pow(Q->y,2);
#endif
};

void Matr2Quat(Ldoub* C, quaternion* Q)
{
	Q->w = 0.5*sqrt(C[index_3(3, 0, 0)] + C[index_3(3, 1, 1)] + C[index_3(3, 2, 2)] + 1);
	if (Q->w == 0)
	{
		Q->x = sqrt(0.5*(C[index_3(3, 0, 0)]+1));
		Q->y = sqrt(0.5*(C[index_3(3, 1, 1)]+1));
		Q->z = sqrt(0.5*(C[index_3(3, 2, 2)]+1));
	}
	else
	{
		Q->x = (C[index_3(3, 2, 1)] - C[index_3(3, 1, 2)])/(4.*Q->w);
		Q->y = (C[index_3(3, 0, 2)] - C[index_3(3, 2, 0)])/(4.*Q->w);
		Q->z = (C[index_3(3, 1, 0)] - C[index_3(3, 0, 1)])/(4.*Q->w);
	}
}