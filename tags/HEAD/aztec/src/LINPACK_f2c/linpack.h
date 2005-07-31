/*====================================================================
 * linpack.h
 *====================================================================*/

#ifndef _LINPACK_H_
#define _LINPACK_H_

#include "f2c.h"
#include "blas123.h"

extern int dgedi(doublereal *a, integer *lda, integer *n, integer *ipvt,
	doublereal *det, doublereal *work, integer *job);

extern int dgefa(doublereal *a, integer *lda, integer *n, integer *ipvt, 
	integer *info);

extern int dgetf2(integer *m, integer *n, doublereal *a, integer *lda, 
	integer *ipiv, integer *info);

extern int dgetrf(integer *m, integer *n, doublereal *a, integer *lda, 
	integer *ipiv, integer *info);

extern int dgetri(integer *n, doublereal *a, integer *lda, integer *ipiv,
	doublereal *work, integer *lwork, integer *info);
	
extern int dgetrs(char *trans, integer *n, integer *nrhs, doublereal *a,
	integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
	
extern int dlaic1(integer *job, integer *j, doublereal *x, doublereal *sest, 
	doublereal *w, doublereal *gamma, doublereal *sestpr, doublereal *s, 
	doublereal *c__);

extern doublereal dlamch(char *cmach);
extern int dlamc1(integer *beta, integer *t, logical *rnd, integer *ieee1);
extern int dlamc2(integer *beta, integer *t, logical *rnd, doublereal *eps,
	integer *emin, doublereal *rmin, integer *emax, doublereal *rmax);
extern doublereal dlamc3(doublereal *a, doublereal *b);
extern int dlamc4(integer *emin, doublereal *start, integer *base);
extern int dlamc5(integer *beta, integer *p, integer *emin, logical *ieee,
	integer *emax, doublereal *rmax);

extern int dlaswp(integer *n, doublereal *a, integer *lda, integer *k1, 
	integer *k2, integer *ipiv, integer *incx);

extern int dpotf2(char *uplo, integer *n, doublereal *a, integer *lda, 
	integer *info);

extern int dtrti2(char *uplo, char *diag, integer *n, doublereal *a, 
	integer *lda, integer *info);

extern int dtrtri(char *uplo, char *diag, integer *n, doublereal *a,
	integer *lda, integer *info);

extern int dgeco(doublereal *a, integer *lda, integer *n, integer *ipvt,
	doublereal *rcond, doublereal *z__);

extern integer ilaenv(integer *ispec, char *name__, char *opts, 
	integer *n1, integer *n2, integer *n3, integer *n4);

#endif  /* _LINPACK_H_ */
