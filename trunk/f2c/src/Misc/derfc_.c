#include "f2c.h"

#ifdef KR_headers
extern double erfc();

double derfc_(x) doublereal *x;
#else
extern double IMT_erfc(double);				/* IMT 9Nov97 */
double derfc_(doublereal *x)
#endif
{
return( IMT_erfc(*x) );						/* IMT 9Nov97 */
}
