/* PAK (04/27/98) */

/* blas2.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "blas123.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int dgbmv_(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta,
	 y, incy)
char *trans;
integer *m, *n, *kl, *ku;
doublereal *alpha, *a;
integer *lda;
doublereal *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j;
    static integer ix, iy, jx, jy, kx, ky;

    static char uplo = 'z'; /* not defined, guessing means it should not be "u" or "l" */


/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGBMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n band matrix, with kl sub-diagonals and ku super-diagonals. */

/*  Parameters */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           fol X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. On exit, Y is overwritten by the updated */
/*           vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame(&uplo, "U") && ! lsame(&uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla("DSYMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame(&uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L50: */
		}
		y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		y[j] += temp1 * a[j + j * a_dim1];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[j] += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		y[jy] += temp1 * a[j + j * a_dim1];
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DSYMV . */

} /* dgbmv_ */

/* Subroutine */ int dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, 
	incy)
char *trans;
integer *m, *n;
doublereal *alpha, *a;
integer *lda;
doublereal *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer lenx, leny, i__, j;
    static integer ix, iy, jx, jy, kx, ky;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Parameters */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. 
*/
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the 
*/
/*           updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame(trans, "N") && ! lsame(trans, "T") && ! 
	    lsame(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla("DGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
*/
/*     up the start points in  X  and  Y. */

    if (lsame(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * a[i__ + j * a_dim1];
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * a[i__ + j * a_dim1];
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DGEMV . */

} /* dgemv_ */

/* Subroutine */ int dger_(m, n, alpha, x, incx, y, incy, a, lda)
integer *m, *n;
doublereal *alpha, *x;
integer *incx;
doublereal *y;
integer *incy;
doublereal *a;
integer *lda;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, ix, jy, kx;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGER   performs the rank 1 operation */

/*     A := alpha*x*y' + A, */

/*  where alpha is a scalar, x is an m element vector, y is an n element 
*/
/*  vector and A is an m by n matrix. */

/*  Parameters */
/*  ========== */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. 
*/
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the m */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. On exit, A is */
/*           overwritten by the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla("DGER  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] += x[i__] * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a[i__ + j * a_dim1] += x[ix] * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

/*     End of DGER  . */

} /* dger_ */

/* Subroutine */ int dsbmv_(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
char *uplo;
integer *n, *k;
doublereal *alpha, *a;
integer *lda;
doublereal *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j, l;
    static integer kplus1, ix, iy, jx, jy, kx, ky;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSBMV  performs the matrix-vector  operation */

/*     y := alpha*A*x + beta*y, */

/*  where alpha and beta are scalars, x and y are n element vectors and */
/*  A is an n by n symmetric band matrix, with k super-diagonals. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the band matrix A is being supplied as */
/*           follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  being supplied. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  being supplied. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry, K specifies the number of super-diagonals of the */
/*           matrix A. K must satisfy  0 .le. K. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/*           by n part of the array A must contain the upper triangular */
/*           band part of the symmetric matrix, supplied column by */
/*           column, with the leading diagonal of the matrix in row */
/*           ( k + 1 ) of the array, the first super-diagonal starting at 
*/
/*           position 2 in row k, and so on. The top left k by k triangle 
*/
/*           of the array A is not referenced. */
/*           The following program segment will transfer the upper */
/*           triangular part of a symmetric band matrix from conventional 
*/
/*           full matrix storage to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = K + 1 - J */
/*                    DO 10, I = MAX( 1, J - K ), J */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/*           by n part of the array A must contain the lower triangular */
/*           band part of the symmetric matrix, supplied column by */
/*           column, with the leading diagonal of the matrix in row 1 of 
*/
/*           the array, the first sub-diagonal starting at position 1 in 
*/
/*           row 2, and so on. The bottom right k by k triangle of the */
/*           array A is not referenced. */
/*           The following program segment will transfer the lower */
/*           triangular part of a symmetric band matrix from conventional 
*/
/*           full matrix storage to band storage: */

/*                 DO 20, J = 1, N */
/*                    M = 1 - J */
/*                    DO 10, I = J, MIN( N, J + K ) */
/*                       A( M + I, J ) = matrix( I, J ) */
/*              10    CONTINUE */
/*              20 CONTINUE */

/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           ( k + 1 ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the */
/*           vector y. On exit, Y is overwritten by the updated vector y. 
*/

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k < 0) {
	info = 3;
    } else if (*lda < *k + 1) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla("DSBMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of the array A 
*/
/*     are accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame(uplo, "U")) {

/*        Form  y  when upper triangle of A is stored. */

	kplus1 = *k + 1;
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		l = kplus1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *k;
		i__4 = j - 1;
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
/* L50: */
		}
		y[j] = y[j] + temp1 * a[kplus1 + j * a_dim1] + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		ix = kx;
		iy = ky;
		l = kplus1 - j;
/* Computing MAX */
		i__4 = 1, i__2 = j - *k;
		i__3 = j - 1;
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * a[kplus1 + j * a_dim1] + *alpha * 
			temp2;
		jx += *incx;
		jy += *incy;
		if (j > *k) {
		    kx += *incx;
		    ky += *incy;
		}
/* L80: */
	    }
	}
    } else {

/*        Form  y  when lower triangle of A is stored. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		y[j] += temp1 * a[j * a_dim1 + 1];
		l = 1 - j;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = min(i__4,i__2);
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[j] += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		y[jy] += temp1 * a[j * a_dim1 + 1];
		l = 1 - j;
		ix = jx;
		iy = jy;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = min(i__4,i__2);
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DSBMV . */

} /* dsbmv_ */

/* Subroutine */ int dspmv_(uplo, n, alpha, ap, x, incx, beta, y, incy)
char *uplo;
integer *n;
doublereal *alpha, *ap, *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j, k;
    static integer kk, ix, iy, jx, jy, kx, ky;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSPMV  performs the matrix-vector operation */

/*     y := alpha*A*x + beta*y, */

/*  where alpha and beta are scalars, x and y are n element vectors and */
/*  A is an n by n symmetric matrix, supplied in packed form. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the matrix A is supplied in the packed */
/*           array AP as follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  supplied in AP. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  supplied in AP. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/*           and a( 2, 2 ) respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/*           and a( 3, 1 ) respectively, and so on. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. On exit, Y is overwritten by the updated */
/*           vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --y;
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 6;
    } else if (*incy == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla("DSPMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of the array AP 
*/
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    kk = 1;
    if (lsame(uplo, "U")) {

/*        Form  y  when AP contains the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		k = kk;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * ap[k];
		    temp2 += ap[k] * x[i__];
		    ++k;
/* L50: */
		}
		y[j] = y[j] + temp1 * ap[kk + j - 1] + *alpha * temp2;
		kk += j;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = kk + j - 2;
		for (k = kk; k <= i__2; ++k) {
		    y[iy] += temp1 * ap[k];
		    temp2 += ap[k] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * ap[kk + j - 1] + *alpha * temp2;
		jx += *incx;
		jy += *incy;
		kk += j;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when AP contains the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		y[j] += temp1 * ap[kk];
		k = kk + 1;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * ap[k];
		    temp2 += ap[k] * x[i__];
		    ++k;
/* L90: */
		}
		y[j] += *alpha * temp2;
		kk += *n - j + 1;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		y[jy] += temp1 * ap[kk];
		ix = jx;
		iy = jy;
		i__2 = kk + *n - j;
		for (k = kk + 1; k <= i__2; ++k) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * ap[k];
		    temp2 += ap[k] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
		kk += *n - j + 1;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DSPMV . */

} /* dspmv_ */

/* Subroutine */ int dspr_(uplo, n, alpha, x, incx, ap)
char *uplo;
integer *n;
doublereal *alpha, *x;
integer *incx;
doublereal *ap;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
    static integer kk, ix, jx, kx;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSPR    performs the symmetric rank 1 operation */

/*     A := alpha*x*x' + A, */

/*  where alpha is a real scalar, x is an n element vector and A is an */
/*  n by n symmetric matrix, supplied in packed form. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the matrix A is supplied in the packed */
/*           array AP as follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  supplied in AP. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  supplied in AP. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/*           and a( 2, 2 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the upper triangular part of the */
/*           updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/*           and a( 3, 1 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the lower triangular part of the */
/*           updated matrix. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    }
    if (info != 0) {
	xerbla("DSPR  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set the start point in X if the increment is not unity. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of the array AP 
*/
/*     are accessed sequentially with one pass through AP. */

    kk = 1;
    if (lsame(uplo, "U")) {

/*        Form  A  when upper triangle is stored in AP. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.) {
		    temp = *alpha * x[j];
		    k = kk;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			ap[k] += x[i__] * temp;
			++k;
/* L10: */
		    }
		}
		kk += j;
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    ix = kx;
		    i__2 = kk + j - 1;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] += x[ix] * temp;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
		kk += j;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when lower triangle is stored in AP. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.) {
		    temp = *alpha * x[j];
		    k = kk;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			ap[k] += x[i__] * temp;
			++k;
/* L50: */
		    }
		}
		kk = kk + *n - j + 1;
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    ix = jx;
		    i__2 = kk + *n - j;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] += x[ix] * temp;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of DSPR  . */

} /* dspr_ */

/* Subroutine */ int dspr2(uplo, n, alpha, x, incx, y, incy, ap)
char *uplo;
integer *n;
doublereal *alpha, *x;
integer *incx;
doublereal *y;
integer *incy;
doublereal *ap;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j, k;
    static integer kk, ix, iy, jx, jy, kx, ky;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSPR2  performs the symmetric rank 2 operation */

/*     A := alpha*x*y' + alpha*y*x' + A, */

/*  where alpha is a scalar, x and y are n element vectors and A is an */
/*  n by n symmetric matrix, supplied in packed form. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the matrix A is supplied in the packed */
/*           array AP as follows: */

/*              UPLO = 'U' or 'u'   The upper triangular part of A is */
/*                                  supplied in AP. */

/*              UPLO = 'L' or 'l'   The lower triangular part of A is */
/*                                  supplied in AP. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/*           and a( 2, 2 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the upper triangular part of the */
/*           updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular part of the symmetric matrix */
/*           packed sequentially, column by column, so that AP( 1 ) */
/*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/*           and a( 3, 1 ) respectively, and so on. On exit, the array */
/*           AP is overwritten by the lower triangular part of the */
/*           updated matrix. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;
    --y;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla("DSPR2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
*/
/*     unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of the array AP 
*/
/*     are accessed sequentially with one pass through AP. */

    kk = 1;
    if (lsame(uplo, "U")) {

/*        Form  A  when upper triangle is stored in AP. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    k = kk;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
			++k;
/* L10: */
		    }
		}
		kk += j;
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = kk + j - 1;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
		kk += j;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when lower triangle is stored in AP. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    k = kk;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			ap[k] = ap[k] + x[i__] * temp1 + y[i__] * temp2;
			++k;
/* L50: */
		    }
		}
		kk = kk + *n - j + 1;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = kk + *n - j;
		    for (k = kk; k <= i__2; ++k) {
			ap[k] = ap[k] + x[ix] * temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of DSPR2 . */

} /* dspr2_ */

/* Subroutine */ int dsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
char *uplo;
integer *n;
doublereal *alpha, *a;
integer *lda;
doublereal *x;
integer *incx;
doublereal *beta, *y;
integer *incy;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j;
    static integer ix, iy, jx, jy, kx, ky;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSYMV  performs the matrix-vector  operation */

/*     y := alpha*A*x + beta*y, */

/*  where alpha and beta are scalars, x and y are n element vectors and */
/*  A is an n by n symmetric matrix. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the array A is to be referenced as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper 
*/
/*           triangular part of the symmetric matrix and the strictly */
/*           lower triangular part of A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower 
*/
/*           triangular part of the symmetric matrix and the strictly */
/*           upper triangular part of A is not referenced. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. On exit, Y is overwritten by the updated */
/*           vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < max(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla("DSYMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame(uplo, "U")) {

/*        Form  y  when A is stored in upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L50: */
		}
		y[j] = y[j] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		y[jy] = y[jy] + temp1 * a[j + j * a_dim1] + *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when A is stored in lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[j];
		temp2 = 0.;
		y[j] += temp1 * a[j + j * a_dim1];
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    y[i__] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		}
		y[j] += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * x[jx];
		temp2 = 0.;
		y[jy] += temp1 * a[j + j * a_dim1];
		ix = jx;
		iy = jy;
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * a[i__ + j * a_dim1];
		    temp2 += a[i__ + j * a_dim1] * x[ix];
/* L110: */
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DSYMV . */

} /* dsymv_ */

/* Subroutine */ int dsyr_(uplo, n, alpha, x, incx, a, lda)
char *uplo;
integer *n;
doublereal *alpha, *x;
integer *incx;
doublereal *a;
integer *lda;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j;
    static integer ix, jx, kx;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSYR   performs the symmetric rank 1 operation */

/*     A := alpha*x*x' + A, */

/*  where alpha is a real scalar, x is an n element vector and A is an */
/*  n by n symmetric matrix. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the array A is to be referenced as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper 
*/
/*           triangular part of the symmetric matrix and the strictly */
/*           lower triangular part of A is not referenced. On exit, the */
/*           upper triangular part of the array A is overwritten by the */
/*           upper triangular part of the updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower 
*/
/*           triangular part of the symmetric matrix and the strictly */
/*           upper triangular part of A is not referenced. On exit, the */
/*           lower triangular part of the array A is overwritten by the */
/*           lower triangular part of the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*lda < max(1,*n)) {
	info = 7;
    }
    if (info != 0) {
	xerbla("DSYR  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set the start point in X if the increment is not unity. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

    if (lsame(uplo, "U")) {

/*        Form  A  when A is stored in upper triangle. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.) {
		    temp = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[i__] * temp;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    ix = kx;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[ix] * temp;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in lower triangle. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0.) {
		    temp = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[i__] * temp;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    ix = jx;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] += x[ix] * temp;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of DSYR  . */

} /* dsyr_ */

/* Subroutine */ int dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
char *uplo;
integer *n;
doublereal *alpha, *x;
integer *incx;
doublereal *y;
integer *incy;
doublereal *a;
integer *lda;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j;
    static integer ix, iy, jx, jy, kx, ky;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSYR2  performs the symmetric rank 2 operation */

/*     A := alpha*x*y' + alpha*y*x' + A, */

/*  where alpha is a scalar, x and y are n element vectors and A is an n 
*/
/*  by n symmetric matrix. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the upper or lower */
/*           triangular part of the array A is to be referenced as */
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper 
*/
/*           triangular part of the symmetric matrix and the strictly */
/*           lower triangular part of A is not referenced. On exit, the */
/*           upper triangular part of the array A is overwritten by the */
/*           upper triangular part of the updated matrix. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower 
*/
/*           triangular part of the symmetric matrix and the strictly */
/*           upper triangular part of A is not referenced. On exit, the */
/*           lower triangular part of the array A is overwritten by the */
/*           lower triangular part of the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla("DSYR2 ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Set up the start points in X and Y if the increments are not both 
*/
/*     unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

    if (lsame(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[j] != 0. || y[j] != 0.) {
		    temp1 = *alpha * y[j];
		    temp2 = *alpha * x[j];
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 + y[i__] * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0. || y[jy] != 0.) {
		    temp1 = *alpha * y[jy];
		    temp2 = *alpha * x[jx];
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 + y[iy] * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return 0;

/*     End of DSYR2 . */

} /* dsyr2_ */

/* Subroutine */ int dtbmv_(uplo, trans, diag, n, k, a, lda, x, incx)
char *uplo, *trans, *diag;
integer *n, *k;
doublereal *a;
integer *lda;
doublereal *x;
integer *incx;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, l;
    static integer kplus1, ix, jx, kx;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < *k + 1) {
	info = 7;
    } else if (*incx == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla("DTBMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame(diag, "N");

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX   too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame(trans, "N")) {

/*         Form  x := A*x. */

	if (lsame(uplo, "U")) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			temp = x[j];
			l = kplus1 - j;
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__4 = j - 1;
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			    x[i__] += temp * a[l + i__ + j * a_dim1];
/* L10: */
			}
			if (nounit) {
			    x[j] *= a[kplus1 + j * a_dim1];
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			l = kplus1 - j;
/* Computing MAX */
			i__4 = 1, i__2 = j - *k;
			i__3 = j - 1;
			for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			    x[ix] += temp * a[l + i__ + j * a_dim1];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a[kplus1 + j * a_dim1];
			}
		    }
		    jx += *incx;
		    if (j > *k) {
			kx += *incx;
		    }
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			temp = x[j];
			l = 1 - j;
/* Computing MIN */
			i__1 = *n, i__3 = j + *k;
			i__4 = j + 1;
			for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
			    x[i__] += temp * a[l + i__ + j * a_dim1];
/* L50: */
			}
			if (nounit) {
			    x[j] *= a[j * a_dim1 + 1];
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			l = 1 - j;
/* Computing MIN */
			i__4 = *n, i__1 = j + *k;
			i__3 = j + 1;
			for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
			    x[ix] += temp * a[l + i__ + j * a_dim1];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a[j * a_dim1 + 1];
			}
		    }
		    jx -= *incx;
		    if (*n - j >= *k) {
			kx -= *incx;
		    }
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame(uplo, "U")) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    l = kplus1 - j;
		    if (nounit) {
			temp *= a[kplus1 + j * a_dim1];
		    }
/* Computing MAX */
		    i__4 = 1, i__1 = j - *k;
		    i__3 = max(i__4,i__1);
		    for (i__ = j - 1; i__ >= i__3; --i__) {
			temp += a[l + i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    kx -= *incx;
		    ix = kx;
		    l = kplus1 - j;
		    if (nounit) {
			temp *= a[kplus1 + j * a_dim1];
		    }
/* Computing MAX */
		    i__4 = 1, i__1 = j - *k;
		    i__3 = max(i__4,i__1);
		    for (i__ = j - 1; i__ >= i__3; --i__) {
			temp += a[l + i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
		    temp = x[j];
		    l = 1 - j;
		    if (nounit) {
			temp *= a[j * a_dim1 + 1];
		    }
/* Computing MIN */
		    i__1 = *n, i__2 = j + *k;
		    i__4 = min(i__1,i__2);
		    for (i__ = j + 1; i__ <= i__4; ++i__) {
			temp += a[l + i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__3 = *n;
		for (j = 1; j <= i__3; ++j) {
		    temp = x[jx];
		    kx += *incx;
		    ix = kx;
		    l = 1 - j;
		    if (nounit) {
			temp *= a[j * a_dim1 + 1];
		    }
/* Computing MIN */
		    i__1 = *n, i__2 = j + *k;
		    i__4 = min(i__1,i__2);
		    for (i__ = j + 1; i__ <= i__4; ++i__) {
			temp += a[l + i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTBMV . */

} /* dtbmv_ */

/* Subroutine */ int dtbsv_(uplo, trans, diag, n, k, a, lda, x, incx)
char *uplo, *trans, *diag;
integer *n, *k;
doublereal *a;
integer *lda;
doublereal *x;
integer *incx;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, l;
    static integer kplus1, ix, jx, kx;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten 
*/
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < *k + 1) {
	info = 7;
    } else if (*incx == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla("DTBSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame(diag, "N");

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed by sequentially with one pass through A. */

    if (lsame(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame(uplo, "U")) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			l = kplus1 - j;
			if (nounit) {
			    x[j] /= a[kplus1 + j * a_dim1];
			}
			temp = x[j];
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__1 = max(i__2,i__3);
			for (i__ = j - 1; i__ >= i__1; --i__) {
			    x[i__] -= temp * a[l + i__ + j * a_dim1];
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    kx -= *incx;
		    if (x[jx] != 0.) {
			ix = kx;
			l = kplus1 - j;
			if (nounit) {
			    x[jx] /= a[kplus1 + j * a_dim1];
			}
			temp = x[jx];
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__1 = max(i__2,i__3);
			for (i__ = j - 1; i__ >= i__1; --i__) {
			    x[ix] -= temp * a[l + i__ + j * a_dim1];
			    ix -= *incx;
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			l = 1 - j;
			if (nounit) {
			    x[j] /= a[j * a_dim1 + 1];
			}
			temp = x[j];
/* Computing MIN */
			i__3 = *n, i__4 = j + *k;
			i__2 = min(i__3,i__4);
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * a[l + i__ + j * a_dim1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    kx += *incx;
		    if (x[jx] != 0.) {
			ix = kx;
			l = 1 - j;
			if (nounit) {
			    x[jx] /= a[j * a_dim1 + 1];
			}
			temp = x[jx];
/* Computing MIN */
			i__3 = *n, i__4 = j + *k;
			i__2 = min(i__3,i__4);
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[ix] -= temp * a[l + i__ + j * a_dim1];
			    ix += *incx;
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A')*x. */

	if (lsame(uplo, "U")) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    l = kplus1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *k;
		    i__4 = j - 1;
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			temp -= a[l + i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    if (nounit) {
			temp /= a[kplus1 + j * a_dim1];
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    l = kplus1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *k;
		    i__3 = j - 1;
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			temp -= a[l + i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= a[kplus1 + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx += *incx;
		    if (j > *k) {
			kx += *incx;
		    }
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    l = 1 - j;
/* Computing MIN */
		    i__1 = *n, i__3 = j + *k;
		    i__4 = j + 1;
		    for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
			temp -= a[l + i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    if (nounit) {
			temp /= a[j * a_dim1 + 1];
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    l = 1 - j;
/* Computing MIN */
		    i__4 = *n, i__1 = j + *k;
		    i__3 = j + 1;
		    for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
			temp -= a[l + i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= a[j * a_dim1 + 1];
		    }
		    x[jx] = temp;
		    jx -= *incx;
		    if (*n - j >= *k) {
			kx -= *incx;
		    }
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTBSV . */

} /* dtbsv_ */

/* Subroutine */ int dtpmv_(uplo, trans, diag, n, ap, x, incx)
char *uplo, *trans, *diag;
integer *n;
doublereal *ap, *x;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
    static integer kk, ix, jx, kx;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTPMV  performs one of the matrix-vector operations */

/*     x := A*x,   or   x := A'*x, */

/*  where x is an n element vector and  A is an n by n unit, or non-unit, 
*/
/*  upper or lower triangular matrix, supplied in packed form. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   x := A*x. */

/*              TRANS = 'T' or 't'   x := A'*x. */

/*              TRANS = 'C' or 'c'   x := A'*x. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/*           respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/*           respectively, and so on. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of 
*/
/*           A are not referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*incx == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla("DTPMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame(diag, "N");

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

    if (lsame(trans, "N")) {

/*        Form  x:= A*x. */

	if (lsame(uplo, "U")) {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			temp = x[j];
			k = kk;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * ap[k];
			    ++k;
/* L10: */
			}
			if (nounit) {
			    x[j] *= ap[kk + j - 1];
			}
		    }
		    kk += j;
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__2 = kk + j - 2;
			for (k = kk; k <= i__2; ++k) {
			    x[ix] += temp * ap[k];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= ap[kk + j - 1];
			}
		    }
		    jx += *incx;
		    kk += j;
/* L40: */
		}
	    }
	} else {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			temp = x[j];
			k = kk;
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[i__] += temp * ap[k];
			    --k;
/* L50: */
			}
			if (nounit) {
			    x[j] *= ap[kk - *n + j];
			}
		    }
		    kk -= *n - j + 1;
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__1 = kk - (*n - (j + 1));
			for (k = kk; k >= i__1; --k) {
			    x[ix] += temp * ap[k];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= ap[kk - *n + j];
			}
		    }
		    jx -= *incx;
		    kk -= *n - j + 1;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame(uplo, "U")) {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= ap[kk];
		    }
		    k = kk - 1;
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += ap[k] * x[i__];
			--k;
/* L90: */
		    }
		    x[j] = temp;
		    kk -= j;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= ap[kk];
		    }
		    i__1 = kk - j + 1;
		    for (k = kk - 1; k >= i__1; --k) {
			ix -= *incx;
			temp += ap[k] * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
		    kk -= j;
/* L120: */
		}
	    }
	} else {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= ap[kk];
		    }
		    k = kk + 1;
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += ap[k] * x[i__];
			++k;
/* L130: */
		    }
		    x[j] = temp;
		    kk += *n - j + 1;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= ap[kk];
		    }
		    i__2 = kk + *n - j;
		    for (k = kk + 1; k <= i__2; ++k) {
			ix += *incx;
			temp += ap[k] * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
		    kk += *n - j + 1;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTPMV . */

} /* dtpmv_ */

/* Subroutine */ int dtpsv_(uplo, trans, diag, n, ap, x, incx)
char *uplo, *trans, *diag;
integer *n;
doublereal *ap, *x;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
    static integer kk, ix, jx, kx;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTPSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix, supplied in packed form. 
*/

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  AP     - DOUBLE PRECISION array of DIMENSION at least */
/*           ( ( n*( n + 1 ) )/2 ). */
/*           Before entry with  UPLO = 'U' or 'u', the array AP must */
/*           contain the upper triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
/*           respectively, and so on. */
/*           Before entry with UPLO = 'L' or 'l', the array AP must */
/*           contain the lower triangular matrix packed sequentially, */
/*           column by column, so that AP( 1 ) contains a( 1, 1 ), */
/*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
/*           respectively, and so on. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of 
*/
/*           A are not referenced, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten 
*/
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --x;
    --ap;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*incx == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla("DTPSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame(diag, "N");

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of AP are */
/*     accessed sequentially with one pass through AP. */

    if (lsame(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame(uplo, "U")) {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			if (nounit) {
			    x[j] /= ap[kk];
			}
			temp = x[j];
			k = kk - 1;
			for (i__ = j - 1; i__ >= 1; --i__) {
			    x[i__] -= temp * ap[k];
			    --k;
/* L10: */
			}
		    }
		    kk -= j;
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			if (nounit) {
			    x[jx] /= ap[kk];
			}
			temp = x[jx];
			ix = jx;
			i__1 = kk - j + 1;
			for (k = kk - 1; k >= i__1; --k) {
			    ix -= *incx;
			    x[ix] -= temp * ap[k];
/* L30: */
			}
		    }
		    jx -= *incx;
		    kk -= j;
/* L40: */
		}
	    }
	} else {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			if (nounit) {
			    x[j] /= ap[kk];
			}
			temp = x[j];
			k = kk + 1;
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * ap[k];
			    ++k;
/* L50: */
			}
		    }
		    kk += *n - j + 1;
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			if (nounit) {
			    x[jx] /= ap[kk];
			}
			temp = x[jx];
			ix = jx;
			i__2 = kk + *n - j;
			for (k = kk + 1; k <= i__2; ++k) {
			    ix += *incx;
			    x[ix] -= temp * ap[k];
/* L70: */
			}
		    }
		    jx += *incx;
		    kk += *n - j + 1;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame(uplo, "U")) {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    k = kk;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= ap[k] * x[i__];
			++k;
/* L90: */
		    }
		    if (nounit) {
			temp /= ap[kk + j - 1];
		    }
		    x[j] = temp;
		    kk += j;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    i__2 = kk + j - 2;
		    for (k = kk; k <= i__2; ++k) {
			temp -= ap[k] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= ap[kk + j - 1];
		    }
		    x[jx] = temp;
		    jx += *incx;
		    kk += j;
/* L120: */
		}
	    }
	} else {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    k = kk;
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= ap[k] * x[i__];
			--k;
/* L130: */
		    }
		    if (nounit) {
			temp /= ap[kk - *n + j];
		    }
		    x[j] = temp;
		    kk -= *n - j + 1;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    i__1 = kk - (*n - (j + 1));
		    for (k = kk; k >= i__1; --k) {
			temp -= ap[k] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= ap[kk - *n + j];
		    }
		    x[jx] = temp;
		    jx -= *incx;
		    kk -= *n - j + 1;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTPSV . */

} /* dtpsv_ */

/* Subroutine */ int dtrmv_(uplo, trans, diag, n, a, lda, x, incx)
char *uplo, *trans, *diag;
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
integer *incx;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j;
    static integer ix, jx, kx;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRMV  performs one of the matrix-vector operations */

/*     x := A*x,   or   x := A'*x, */

/*  where x is an n element vector and  A is an n by n unit, or non-unit, 
*/
/*  upper or lower triangular matrix. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   x := A*x. */

/*              TRANS = 'T' or 't'   x := A'*x. */

/*              TRANS = 'C' or 'c'   x := A'*x. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper 
*/
/*           triangular matrix and the strictly lower triangular part of 
*/
/*           A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower 
*/
/*           triangular matrix and the strictly upper triangular part of 
*/
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of 
*/
/*           A are not referenced either, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element vector x. On exit, X is overwritten with the */
/*           tranformed vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla("DTRMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame(diag, "N");

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame(trans, "N")) {

/*        Form  x := A*x. */

	if (lsame(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
/* L10: */
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[i__] += temp * a[i__ + j * a_dim1];
/* L50: */
			}
			if (nounit) {
			    x[j] *= a[j + j * a_dim1];
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[ix] += temp * a[i__ + j * a_dim1];
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a[j + j * a_dim1];
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= *incx;
			temp += a[i__ + j * a_dim1] * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a[i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += *incx;
			temp += a[i__ + j * a_dim1] * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRMV . */

} /* dtrmv_ */

/* Subroutine */ int dtrsv_(uplo, trans, diag, n, a, lda, x, incx)
char *uplo, *trans, *diag;
integer *n;
doublereal *a;
integer *lda;
doublereal *x;
integer *incx;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j;
    static integer ix, jx, kx;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRSV  solves one of the systems of equations */

/*     A*x = b,   or   A'*x = b, */

/*  where b and x are n element vectors and A is an n by n unit, or */
/*  non-unit, upper or lower triangular matrix. */

/*  No test for singularity or near-singularity is included in this */
/*  routine. Such tests must be performed before calling this routine. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix is an upper or */
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the equations to be solved as */
/*           follows: */

/*              TRANS = 'N' or 'n'   A*x = b. */

/*              TRANS = 'T' or 't'   A'*x = b. */

/*              TRANS = 'C' or 'c'   A'*x = b. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit */
/*           triangular as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the order of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/*           upper triangular part of the array A must contain the upper 
*/
/*           triangular matrix and the strictly lower triangular part of 
*/
/*           A is not referenced. */
/*           Before entry with UPLO = 'L' or 'l', the leading n by n */
/*           lower triangular part of the array A must contain the lower 
*/
/*           triangular matrix and the strictly upper triangular part of 
*/
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u', the diagonal elements of 
*/
/*           A are not referenced either, but are assumed to be unity. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, n ). */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the n */
/*           element right-hand side vector b. On exit, X is overwritten 
*/
/*           with the solution vector x. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */


/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;

    /* Function Body */
    info = 0;
    if (! lsame(uplo, "U") && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla("DTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame(diag, "N");

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			if (nounit) {
			    x[j] /= a[j + j * a_dim1];
			}
			temp = x[j];
			for (i__ = j - 1; i__ >= 1; --i__) {
			    x[i__] -= temp * a[i__ + j * a_dim1];
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			if (nounit) {
			    x[jx] /= a[j + j * a_dim1];
			}
			temp = x[jx];
			ix = jx;
			for (i__ = j - 1; i__ >= 1; --i__) {
			    ix -= *incx;
			    x[ix] -= temp * a[i__ + j * a_dim1];
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			if (nounit) {
			    x[j] /= a[j + j * a_dim1];
			}
			temp = x[j];
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    x[i__] -= temp * a[i__ + j * a_dim1];
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			if (nounit) {
			    x[jx] /= a[j + j * a_dim1];
			}
			temp = x[jx];
			ix = jx;
			i__2 = *n;
			for (i__ = j + 1; i__ <= i__2; ++i__) {
			    ix += *incx;
			    x[ix] -= temp * a[i__ + j * a_dim1];
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= a[i__ + j * a_dim1] * x[i__];
/* L90: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = kx;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp -= a[i__ + j * a_dim1] * x[ix];
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= a[i__ + j * a_dim1] * x[i__];
/* L130: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = kx;
		    i__1 = j + 1;
		    for (i__ = *n; i__ >= i__1; --i__) {
			temp -= a[i__ + j * a_dim1] * x[ix];
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= a[j + j * a_dim1];
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSV . */

} /* dtrsv_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
