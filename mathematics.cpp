#include "mathematics.h"
#include <cmath>
#include "std.h"
void TwoSum(Ldoub a, Ldoub b, Ldoub& res, Ldoub& err, bool isNull)
{
	/*a + b = res + err*/
	res = a+b;
	Ldoub b_virt = res - a;
	Ldoub a_virt = res - b_virt;
	Ldoub b_roundoff = b - b_virt;
	Ldoub a_roundoff = a - a_virt;
	Ldoub y = a_roundoff + b_roundoff;
	err = (isNull)? y : err + y;
}

void Split(Ldoub a, int s, Ldoub& a_hi, Ldoub& a_lo)
{
	Ldoub c = (pow(2., s) +  1)*a;
	Ldoub a_big = c - a;
	a_hi = c - a_big;
	a_lo = a - a_hi;
}

void TwoProduct(Ldoub a, Ldoub b, Ldoub& res, Ldoub& err)
{
	/*a * b = res + err*/
	res = a*b;
	Ldoub a_hi, a_low, b_hi, b_low;
	Split(a, 12, a_hi, a_low);
	Split(b, 12, b_hi, b_low);
	Ldoub err1, err2, err3;
	err1 = res - (a_hi*b_hi);
	err2 = err1 - (a_low*b_hi);
	err3 = err2 - (a_hi*b_low);
	err += ((a_low*b_low) - err3);
}

