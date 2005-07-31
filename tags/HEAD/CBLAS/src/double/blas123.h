/*====================================================================
 * blas123.h
 *====================================================================*/

/* BLAS Level 1 function prototypes */

#ifndef _BLAS_123_H_
#define _BLAS_123_H_

#include "f2c.h"

/* BLAS Level 1 routines */

extern doublereal dasum_(integer *n, doublereal *dx, integer *incx);

extern int daxpy_(integer *n, doublereal *da, doublereal *dx, integer *incx,
	doublereal *dy, integer *incy);

extern int dcopy_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy);

extern doublereal ddot_(integer *n, doublereal *dx, integer *incx,
	doublereal *dy, integer *incy);

extern doublereal dmach_(integer *job);

extern doublereal dnrm2_(integer *n, doublereal *x, integer *incx);

extern int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s);

extern int drotg_(doublereal *da, doublereal *db, doublereal *c__,
	doublereal *s);

extern int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx);
	
extern int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy);

extern integer idamax_(integer *n, doublereal *dx, integer *incx);

/* BLAS Level 2 routines */

extern int dgbmv_(char *trans, integer *m, integer *n, integer *kl, 
 	integer *ku, doublereal *alpha, doublereal *a, integer *lda, 
 	doublereal *x, integer *incx, doublereal *beta, doublereal *y, 
 	integer *incy);

extern int dgemv_(char *trans, integer *m, integer *n, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy);

extern int dger_(integer *m, integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *y, integer *incy, doublereal *a, 
	integer *lda);

extern int dsbmv_(char *uplo, integer *n, integer *k, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy);

extern int dspmv_(char *uplo, integer *n, doublereal *alpha, doublereal *ap, 
	doublereal *x, integer *incx, doublereal *beta, doublereal *y ,
	integer *incy);

extern int dspr_(char *uplo, integer *n, doublereal *alpha, doublereal *x,
	integer *incx, doublereal *ap);
	
extern int dspr2_(char *uplo, integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *y, integer *incy, doublereal *ap);

extern int dsymv_(char *uplo, integer *n, doublereal *alpha, doublereal *a,
	integer *lda, doublereal *x, integer *incx, doublereal *beta, 
	doublereal *y, integer *incy);

extern int dsyr_(char *uplo, integer *n, doublereal *alpha, doublereal *x,
	integer *incx, doublereal *a, integer *lda);

extern int dsyr2_(char *uplo, integer *n, doublereal *alpha, doublereal *x,
 	integer *incx, doublereal *y, integer *incy, doublereal *a, 
 	integer *lda);

extern int dtbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *x, integer *incx);

extern int dtbsv_(char *uplo, char *trans, char *diag, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *x, integer *incx);

extern int dtpmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *ap, 
	doublereal *x, integer *incx);
	
extern int dtpsv_(char *uplo, char *trans, char *diag, integer *n, doublereal *ap, 
	doublereal *x, integer *incx);
	
extern int dtrmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a,
	integer *lda, doublereal *x, integer *incx);
	
extern int dtrsv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a, 
	integer *lda, doublereal *x, integer *incx);
	
/* BLAS Level 3 routines */

extern int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k,
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *beta, doublereal *c__, integer *ldc);

extern int dsymm_(char *side, char *uplo, integer *m, integer *n, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, 
	doublereal *c__, integer *ldc);
	
extern int dsyr2k_(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);
	
extern int dsyrk_(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *beta, 
	doublereal *c__, integer *ldc);

extern int dtrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, 
	integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb);

extern int dtrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, 
	integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b,
	integer *ldb);

#endif  /* _BLAS_123_H_ */
