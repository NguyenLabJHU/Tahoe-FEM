#include <math.h>
#include "fortran_names.h"

double FORTRAN_NAME(sqrt_c)(double* a);
double FORTRAN_NAME(sqrt_c)(double* a) 
{ 
	*a = sqrt(*a);
	return *a;
}

	 
