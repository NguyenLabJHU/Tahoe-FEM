#include "f2c.h"

#ifdef KR_headers
double erf();
double erf_(x) real *x;
#else
extern double IMT_erf(double);					/* IMT 9Nov97 */
double erf_(real *x)
#endif
{
return( IMT_erf(*x) );								/* IMT 9Nov97 */
}
