/* $Id: main_c.c,v 1.2 2003-03-09 04:17:13 paklein Exp $ */
#include <stdio.h>
#include <math.h>
#include "fortran_names.h"

extern void FORTRAN_NAME(double_it)(double*, int*);
extern void FORTRAN_NAME(sqrt_f)(double*);
extern void FORTRAN_NAME(sqrt_f_c)(double*);
extern void FORTRAN_NAME(print_double)(double*);
extern void FORTRAN_NAME(print_int)(int*);

static int d_check(double a_ref, double a_test) {

	double tol = 1.0e-14;

	if (fabs(a_ref) < tol)
	{
		if (fabs(a_ref - a_test) < tol)
			return 1;
		else
			return 0;
	}
	else
	{
		if (fabs(a_ref - a_test)/a_ref < tol)
			return 1;
		else
			return 0;
	}
};

static int i_check(int a_ref, int a_test)
{
	if (a_ref == a_test)
		return 1;
	else
		return 0;
};

int main(int argc, char** argv)
{
	double a = 1.1;
	int i = 2;

	double d1, d2;
	int i1, i2;
	
	d1 = d2 = cos(0.1);
	i1 = i2 = 94117;
	
	printf("start: d = %15.12lf, i = %5d\n", d1, i1);
	
	FORTRAN_NAME(double_it)(&d1, &i1);
	d2 *= 2.0;
	i2 *= 2;
	printf("double an int and a double: ");
	if (d_check(d1,d2) && i_check(i1,i2))
		printf("PASS\n");
	else
		printf("FAIL\n");
	
	FORTRAN_NAME(sqrt_f)(&d1);
	d2 = sqrt(d2);
	printf("sqrt of double: ");
	if (d_check(d1,d2))
		printf("PASS\n");
	else
		printf("FAIL\n");

	FORTRAN_NAME(sqrt_f_c)(&d1);
	d2 = sqrt(d2);
	printf("sqrt of double in C called from Fortran: ");
	if (d_check(d1,d2))
		printf("PASS\n");
	else
		printf("FAIL\n");

	printf("end: d = %15.12lf, i = %5d\n", d1, i1);

	printf("print values from Fortran functions:\n");
	FORTRAN_NAME(print_double)(&d1);
	FORTRAN_NAME(print_int)(&i1);
		
	return 0;
}
