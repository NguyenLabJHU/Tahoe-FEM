#include <math.h>
#include "fortran_names.h"

double FORTRAN_NAME(sqrt_c)(double* a);
double FORTRAN_NAME(sqrt_c)(double* a) 

/*
double SQRT_C(double* a);
double SQRT_C(double* a)
*/
{ 
	*a = sqrt(*a);
	return *a;
}

	 
