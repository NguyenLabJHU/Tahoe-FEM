#include "f2c.h"

#ifdef KR_headers
double erf();
double derf_(x) doublereal *x;
#else
extern double IMT_erf(double);				/* IMT 9Nov97 */
double derf_(doublereal *x)
#endif
{
return( IMT_erf(*x) );						/* IMT 9Nov97 */
}
