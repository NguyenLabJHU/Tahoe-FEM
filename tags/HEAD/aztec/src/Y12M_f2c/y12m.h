/*====================================================================
 * y12m.h
 *====================================================================*/

#ifndef _Y_12_M_H_
#define _Y_12_M_H_

#include "f2c.h"

doublereal y12cck(integer *i__, integer *j);

extern int y12mae(integer *n, integer *z__, real *a, integer *snr, integer *nn, 
	integer *rnr, integer *nn1, real *pivot, integer *ha, integer *iha, 
	real *aflag, integer *iflag, real *b, integer *ifail);

extern int y12maf(integer *n, integer *z__, doublereal *a, integer *snr, 
	integer *nn, integer *rnr, integer *nn1, doublereal *pivot, 
	integer *ha, integer *iha, doublereal *aflag, integer *iflag,
	doublereal *b, integer *ifail);

extern int y12mbe(integer *n, integer *z__, real *a, integer *snr, 
	integer *nn, integer *rnr, integer *nn1, integer *ha, integer *iha, 
	real *aflag, integer *iflag, integer *ifail);

extern int y12mbf(integer *n, integer *z__, doublereal *a, integer *snr, 
	integer *nn, integer *rnr, integer *nn1, integer *ha, integer *iha, 
	doublereal *aflag, integer *iflag, integer *ifail);

extern int y12mce(integer *n, integer *z__, real *a, integer *snr, integer *nn, 
	integer *rnr, integer *nn1, real *pivot, real *b, integer *ha, integer *iha,
	real *aflag, integer *iflag, integer *ifail);

extern int y12mcf(integer *n, integer *z__, doublereal *a, integer *snr, 
	integer *nn, integer *rnr, integer *nn1, doublereal *pivot, doublereal *b,
	integer *ha, integer *iha, doublereal *aflag, integer *iflag, integer *ifail);

extern int y12mde(integer *n, real *a, integer *nn, real *b, real *pivot,
	integer *snr, integer *ha, integer *iha, integer *iflag, integer *ifail);

extern int y12mdf(integer *n, doublereal *a, integer *nn, doublereal *b, 
	doublereal *pivot, integer *snr, integer *ha, integer *iha, 
	integer *iflag, integer *ifail);

extern int y12mfe(integer *n, real *a, integer *snr, integer *nn, integer *rnr, 
	integer *nn1, real *a1, integer *sn, integer *nz, integer *ha, integer *iha,
	real *b, real *b1, real *x, real *y, real *aflag, integer *iflag, integer *ifail);

extern int y12mge(integer *n, integer *nn, real *a, integer *snr, real *w, 
	real *pivot, real *anorm, real *rcond, integer *iha, integer *ha, 
	integer *iflag, integer *ifail);

extern int y12mhe(integer *n, integer *nz, real *a, integer *snr, real *work, 
	real *anorm);

#endif  /* _Y_12_M_H_ */
