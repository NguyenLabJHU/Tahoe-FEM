/* PAK (04/27/98) */

/* y12cck.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <time.h>

#include "y12m.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

doublereal y12cck(i__, j)
integer *i__, *j;
{
    /* Local variables */
/* original f2c code
    static integer k;
    extern integer itime_();

    ret_val = (real) itime_(&k) / (float)60.;
    return ret_val;
*/

	return (((doublereal) clock())/CLOCKS_PER_SEC);

} /* y12cck_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */

