#include "mathematics.h"
#include <cmath>
void TwoSum(double a, double b, double& res, double& err, bool isNull)
{
	/*a + b = res + err*/
	res = a+b;
	double b_virt = res - a;
	double a_virt = res - b_virt;
	double b_roundoff = b - b_virt;
	double a_roundoff = a - a_virt;
	double y = a_roundoff + b_roundoff;
	err = (isNull)? y : err + y;
}

void TwoProduct(double a, double b, double& res, double& err)
{
	/*a * b = res + err*/
	res = a*b;
	double a_hi, a_low, b_hi, b_low;
	Split(a, 12, a_hi, a_low);
	Split(b, 12, b_hi, b_low);
	double err1, err2, err3;
	err1 = res - (a_hi*b_hi);
	err2 = err1 - (a_low*b_hi);
	err3 = err2 - (a_hi*b_low);
	err += ((a_low*b_low) - err3);
}

void Split(double a, int s, double& a_hi, double& a_lo)
{
	double c = (pow(2., s) +  1)*a;
	double a_big = c - a;
	a_hi = c - a_big;
	a_lo = a - a_hi;
}
