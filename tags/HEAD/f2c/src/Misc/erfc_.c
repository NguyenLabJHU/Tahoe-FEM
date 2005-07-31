#include "f2c.h"

#ifdef KR_headers
double erfc();
double erfc_(x) real *x;
#else
extern double IMT_erfc(double);				/* IMT 9Nov97 */
double erfc_(real *x)
#endif
{
return( IMT_erfc(*x) );						/* IMT 9Nov97 */
}
