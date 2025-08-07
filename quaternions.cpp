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

quaternion& quaternion::operator / (const Ldoub& val)  //перегрузка оператора деления на число
{
	this->w /= val;
	this->x /= val;
	this->y /= val;
	this->z /= val;
	return *this;
}			

quaternion quaternion::operator * (const Ldoub& val ) //перегрузка оператора умножения на число
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
	return pow(this->w, 2) + pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2);
}