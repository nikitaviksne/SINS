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

quaternion& operator / (const Ldoub& val);  //перегрузка оператора деления на число
{
	this->w /= val;
	this->x /= val;
	this->y /= val;
	this->z /= val;
	return *this;
}			

quaternion operator * (const Ldoub& val ) //перегрузка оператора умножения на число
{
	this->w *= val;
	this->x *= val;
	this->y *= val;
	this->z *= val;
	return *this;
}

quaternion operator * (const quaternion& quat ) const //перегрузка оператора умножения кватернионов
{
	Ldoub w, x, y, z;
	w = this->w * quat.w - this->x * quat.x - this->y * quat.y - this->z * quat.z;
	x = this->w * quat.x + this->x * quat.w + this->y * quat.z - this->z * quat.y;
	y = this->w * quat.y - this->x * quat.z + this->y * quat.w + this->z * quat.x;
	z = this->w * quat.z + this->x * quat.y - this->y * quat.x + this->z * quat.w; 
	return quaternion(w,x,y,z);
}

quaternion operator = (const quaternion& quat) //перегрузка оператора присваивания
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
	return pow(this->w, 2) + pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2);
}


void MulQuat(Ldoub* A, Ldoub* B, Ldoub* C)
{
/*� ����������� 4 ��������, ������� ������ �������� ������ 4*/
	C[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
	C[1] = A[0]*B[1] + A[1]*B[0] + A[2]*B[3] - A[3]*B[2];
	C[2] = A[0]*B[2] - A[1]*B[3] + A[2]*B[0] + A[3]*B[1];
	C[3] = A[0]*B[3] + A[1]*B[2] - A[2]*B[1] + A[3]*B[0];
}

void Quat2Matr(Ldoub* Q, Ldoub* C)
{
	/*Q - ����������, C - ������� ���� ���������*/
	C[0] = pow(Q[0],2) + pow(Q[1],2) + pow(Q[2],2) + pow(Q[3],2); //c11
	C[1] = 2 * (Q[1]*Q[2] - Q[0]*Q[3]); //c12
	C[2] = 2 * (Q[1]*Q[3] + Q[0]*Q[2]); //c13
	C[3] = 2 * (Q[1]*Q[2] + Q[0]*Q[3]); // c21
	C[4] = pow(Q[0],2) + pow(Q[2],2) - pow(Q[1],2) - pow(Q[3],2); //c22
	C[5] = 2 * (Q[2]*Q[3] - Q[0]*Q[1]); // c23
	C[6] = 2 * (Q[1]*Q[3] - Q[0]*Q[2]); //c31
	C[7] = 2 * (Q[2]*Q[3] + Q[0]*Q[1]); //c32
	C[8] = pow(Q[0],2) + pow(Q[3],2) - pow(Q[1],2) - pow(Q[2],2);//c33

}

void Matr2Quat(Ldoub* C, Ldoub* Q)
{
	Q[0] = 0.5*sqrt(C[index_3(3, 0,0)] + C[index_3(3, 1,1)] + C[index_3(3,2,2)] + 1);
	if (Q[0] !=0)
	{
		Q[1] = (C[index_3(3,2,1)] - C[index_3(3,1,2)])/(4*Q[0]);
		Q[2] = (C[index_3(3,0,2)] - C[index_3(3,2,0)])/(4*Q[0]);
		Q[3] = (C[index_3(3,1,0)] - C[index_3(3,0,1)])/(4*Q[0]);
	}
	else
	{
		Q[1] = sqrt((C[index_3(3,0,0)] + 1)/2);
		Q[2] = sqrt((C[index_3(3,1,1)] + 1)/2);
		Q[3] = sqrt((C[index_3(3,2,2)] + 1)/2);
	}


}