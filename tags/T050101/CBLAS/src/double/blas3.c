/* PAK (04/27/98) */

/* blas3.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "blas123.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

/* Subroutine */ int dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, 
	beta, c__, ldc)
char *transa, *transb;
integer *m, *n, *k;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
doublereal *beta, *c__;
integer *ldc;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer info;
    static logical nota, notb;
    static doublereal temp;
    static integer i__, j, l, ncola;
    static integer nrowa, nrowb;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DGEMM  performs one of the matrix-matrix operations */

/*     C := alpha*op( A )*op( B ) + beta*C, */

/*  where  op( X ) is one of */

/*     op( X ) = X   or   op( X ) = X', */

/*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
*/
/*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. 
*/

/*  Parameters */
/*  ========== */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n',  op( A ) = A. */

/*              TRANSA = 'T' or 't',  op( A ) = A'. */

/*              TRANSA = 'C' or 'c',  op( A ) = A'. */

/*           Unchanged on exit. */

/*  TRANSB - CHARACTER*1. */
/*           On entry, TRANSB specifies the form of op( B ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSB = 'N' or 'n',  op( B ) = B. */

/*              TRANSB = 'T' or 't',  op( B ) = B'. */

/*              TRANSB = 'C' or 'c',  op( B ) = B'. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry,  M  specifies  the number  of rows  of the  matrix 
*/
/*           op( A )  and of the  matrix  C.  M  must  be at least  zero. 
*/
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N  specifies the number  of columns of the matrix 
*/
/*           op( B ) and the number of columns of the matrix C. N must be 
*/
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry,  K  specifies  the number of columns of the matrix 
*/
/*           op( A ) and the number of rows of the matrix op( B ). K must 
*/
/*           be at least  zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
*/
/*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k 
*/
/*           part of the array  A  must contain the matrix  A,  otherwise 
*/
/*           the leading  k by m  part of the array  A  must contain  the 
*/
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then 
*/
/*           LDA must be at least  max( 1, m ), otherwise  LDA must be at 
*/
/*           least  max( 1, k ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is 
*/
/*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n 
*/
/*           part of the array  B  must contain the matrix  B,  otherwise 
*/
/*           the leading  n by k  part of the array  B  must contain  the 
*/
/*           matrix B. */
/*           Unchanged on exit. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared 
*/
/*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then 
*/
/*           LDB must be at least  max( 1, k ), otherwise  LDB must be at 
*/
/*           least  max( 1, n ). */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is 
*/
/*           supplied as zero then C need not be set on input. */
/*           Unchanged on exit. */

/*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
/*           Before entry, the leading  m by n  part of the array  C must 
*/
/*           contain the matrix  C,  except when  beta  is zero, in which 
*/
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n  matrix 
*/
/*           ( alpha*op( A )*op( B ) + beta*C ). */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared 
*/
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least 
*/
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not 
*/
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows 
*/
/*     and  columns of  A  and the  number of  rows  of  B  respectively. 
*/

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;

    /* Function Body */
    nota = lsame(transa, "N");
    notb = lsame(transb, "N");
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! lsame(transa, "C") && ! lsame(transa, "T")) {
	info = 1;
    } else if (! notb && ! lsame(transb, "C") && ! lsame(transb, 
	    "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < max(1,nrowa)) {
	info = 8;
    } else if (*ldb < max(1,nrowb)) {
	info = 10;
    } else if (*ldc < max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla("DGEMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

/*     And if  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L50: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[l + j * b_dim1] != 0.) {
			temp = *alpha * b[l + j * b_dim1];
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
/* L100: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {

/*           Form  C := alpha*A*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L130: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b[j + l * b_dim1] != 0.) {
			temp = *alpha * b[j + l * b_dim1];
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
/* L180: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }

    return 0;

/*     End of DGEMM . */

} /* dgemm_ */

/* Subroutine */ int dsyr2k_(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, 
	c__, ldc)
char *uplo, *trans;
integer *n, *k;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
doublereal *beta, *c__;
integer *ldc;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j, l;
    static integer nrowa;
    static logical upper;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */
/*           lower triangular part of the array C must contain the lower 
*/
/*           triangular part  of the  symmetric matrix  and the strictly 
*/
/*           upper triangular part of C is not referenced.  On exit, the 
*/
/*           lower triangular part of the array  C is overwritten by the 
*/
/*           lower triangular part of the updated matrix. */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared 
*/
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least 
*/
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */


/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame(uplo, "U");

    info = 0;
    if (! upper && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,nrowa)) {
	info = 9;
    } else if (*ldc < max(1,*n)) {
	info = 12;
    }
    if (info != 0) {
	xerbla("DSYR2K", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame(trans, "N")) {

/*        Form  C := alpha*A*B' + alpha*B*A' + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
			temp1 = *alpha * b[j + l * b_dim1];
			temp2 = *alpha * a[j + l * a_dim1];
			i__3 = j;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
				    i__ + l * a_dim1] * temp1 + b[i__ + l * 
				    b_dim1] * temp2;
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
			temp1 = *alpha * b[j + l * b_dim1];
			temp2 = *alpha * a[j + l * a_dim1];
			i__3 = *n;
			for (i__ = j; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
				    i__ + l * a_dim1] * temp1 + b[i__ + l * 
				    b_dim1] * temp2;
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*B + alpha*B'*A + C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
			temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
/* L190: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
				temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ *alpha * temp1 + *alpha * temp2;
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    temp1 = 0.;
		    temp2 = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
			temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
/* L220: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
				temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ *alpha * temp1 + *alpha * temp2;
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return 0;

/*     End of DSYR2K. */

} /* dsyr2k_ */

/* Subroutine */ int dsyrk_(uplo, trans, n, k, alpha, a, lda, beta, c__, ldc)
char *uplo, *trans;
integer *n, *k;
doublereal *alpha, *a;
integer *lda;
doublereal *beta, *c__;
integer *ldc;
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, l;
    static integer nrowa;
    static logical upper;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DSYRK  performs one of the symmetric rank k operations */

/*     C := alpha*A*A' + beta*C, */

/*  or */

/*     C := alpha*A'*A + beta*C, */

/*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix 
*/
/*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix 
*/
/*  in the second case. */

/*  Parameters */
/*  ========== */

/*  UPLO   - CHARACTER*1. */
/*           On  entry,   UPLO  specifies  whether  the  upper  or  lower 
*/
/*           triangular  part  of the  array  C  is to be  referenced  as 
*/
/*           follows: */

/*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C 
*/
/*                                  is to be referenced. */

/*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C 
*/
/*                                  is to be referenced. */

/*           Unchanged on exit. */

/*  TRANS  - CHARACTER*1. */
/*           On entry,  TRANS  specifies the operation to be performed as 
*/
/*           follows: */

/*              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C. */

/*              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C. */

/*              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C. */

/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry,  N specifies the order of the matrix C.  N must be 
*/
/*           at least zero. */
/*           Unchanged on exit. */

/*  K      - INTEGER. */
/*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number 
*/
/*           of  columns   of  the   matrix   A,   and  on   entry   with 
*/
/*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number 
*/
/*           of rows of the matrix  A.  K must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
*/
/*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
*/
/*           part of the array  A  must contain the matrix  A,  otherwise 
*/
/*           the leading  k by n  part of the array  A  must contain  the 
*/
/*           matrix A. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
*/
/*           then  LDA must be at least  max( 1, n ), otherwise  LDA must 
*/
/*           be at least  max( 1, k ). */
/*           Unchanged on exit. */

/*  BETA   - DOUBLE PRECISION. */
/*           On entry, BETA specifies the scalar beta. */
/*           Unchanged on exit. */

/*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
/*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n 
*/
/*           upper triangular part of the array C must contain the upper 
*/
/*           triangular part  of the  symmetric matrix  and the strictly 
*/
/*           lower triangular part of C is not referenced.  On exit, the 
*/
/*           upper triangular part of the array  C is overwritten by the 
*/
/*           upper triangular part of the updated matrix. */
/*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n 
*/
/*           lower triangular part of the array C must contain the lower 
*/
/*           triangular part  of the  symmetric matrix  and the strictly 
*/
/*           upper triangular part of C is not referenced.  On exit, the 
*/
/*           lower triangular part of the array  C is overwritten by the 
*/
/*           lower triangular part of the updated matrix. */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared 
*/
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least 
*/
/*           max( 1, n ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = lsame(uplo, "U");

    info = 0;
    if (! upper && ! lsame(uplo, "L")) {
	info = 1;
    } else if (! lsame(trans, "N") && ! lsame(trans, "T") &&
	     ! lsame(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldc < max(1,*n)) {
	info = 10;
    }
    if (info != 0) {
	xerbla("DSYRK ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame(trans, "N")) {

/*        Form  C := alpha*A*A' + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L100: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != 0.) {
			temp = *alpha * a[j + l * a_dim1];
			i__3 = j;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L150: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (a[j + l * a_dim1] != 0.) {
			temp = *alpha * a[j + l * a_dim1];
			i__3 = *n;
			for (i__ = j; i__ <= i__3; ++i__) {
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*A + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
/* L190: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
/* L220: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = *alpha * temp;
		    } else {
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return 0;

/*     End of DSYRK . */

} /* dsyrk_ */

/* Subroutine */ int dsymm_(side, uplo, m, n, alpha, a, lda, b, ldb, beta, 
	c__, ldc)
char *side, *uplo;
integer *m, *n;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
doublereal *beta, *c__;
integer *ldc;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer info;
    static doublereal temp1, temp2;
    static integer i__, j, k;
    static integer nrowa;
    static logical upper;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */
/*           Before entry, the leading  m by n  part of the array  C must 
*/
/*           contain the matrix  C,  except when  beta  is zero, in which 
*/
/*           case C need not be set on entry. */
/*           On exit, the array  C  is overwritten by the  m by n updated 
*/
/*           matrix. */

/*  LDC    - INTEGER. */
/*           On entry, LDC specifies the first dimension of C as declared 
*/
/*           in  the  calling  (sub)  program.   LDC  must  be  at  least 
*/
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Set NROWA as the number of rows of A. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = c_dim1 + 1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame(side, "L")) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    upper = lsame(uplo, "U");

/*     Test the input parameters. */

    info = 0;
    if (! lsame(side, "L") && ! lsame(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame(uplo, "L")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,nrowa)) {
	info = 7;
    } else if (*ldb < max(1,*m)) {
	info = 9;
    } else if (*ldc < max(1,*m)) {
	info = 12;
    }
    if (info != 0) {
	xerbla("DSYMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (lsame(side, "L")) {

/*        Form  C := alpha*A*B + beta*C. */

	if (upper) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp1 = *alpha * b[i__ + j * b_dim1];
		    temp2 = 0.;
		    i__3 = i__ - 1;
		    for (k = 1; k <= i__3; ++k) {
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
/* L50: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
		    }
/* L60: */
		}
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		for (i__ = *m; i__ >= 1; --i__) {
		    temp1 = *alpha * b[i__ + j * b_dim1];
		    temp2 = 0.;
		    i__2 = *m;
		    for (k = i__ + 1; k <= i__2; ++k) {
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
/* L80: */
		    }
		    if (*beta == 0.) {
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
		    } else {
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
    } else {

/*        Form  C := alpha*B*A + beta*C. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    temp1 = *alpha * a[j + j * a_dim1];
	    if (*beta == 0.) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = temp1 * b[i__ + j * b_dim1];
/* L110: */
		}
	    } else {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] + 
			    temp1 * b[i__ + j * b_dim1];
/* L120: */
		}
	    }
	    i__2 = j - 1;
	    for (k = 1; k <= i__2; ++k) {
		if (upper) {
		    temp1 = *alpha * a[k + j * a_dim1];
		} else {
		    temp1 = *alpha * a[j + k * a_dim1];
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
/* L130: */
		}
/* L140: */
	    }
	    i__2 = *n;
	    for (k = j + 1; k <= i__2; ++k) {
		if (upper) {
		    temp1 = *alpha * a[j + k * a_dim1];
		} else {
		    temp1 = *alpha * a[k + j * a_dim1];
		}
		i__3 = *m;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
/* L150: */
		}
/* L160: */
	    }
/* L170: */
	}
    }

    return 0;

/*     End of DSYMM . */

} /* dsymm_ */

/* Subroutine */ int dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, 
	ldb)
char *side, *uplo, *transa, *diag;
integer *m, *n;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
    static logical lside;
    static integer nrowa;
    static logical upper;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRSM  solves one of the matrix equations */

/*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, */

/*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or 
*/
/*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of 
*/

/*     op( A ) = A   or   op( A ) = A'. */

/*  The matrix X is overwritten on B. */

/*  Parameters */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry, SIDE specifies whether op( A ) appears on the left 
*/
/*           or right of X as follows: */

/*              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */

/*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix A is an upper or 
*/
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n'   op( A ) = A. */

/*              TRANSA = 'T' or 't'   op( A ) = A'. */

/*              TRANSA = 'C' or 'c'   op( A ) = A'. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit triangular 
*/
/*           as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of B. M must be at 
*/
/*           least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of B.  N must be 
*/
/*           at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is 
*/
/*           zero then  A is not referenced and  B need not be set before 
*/
/*           entry. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m 
*/
/*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. 
*/
/*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k 
*/
/*           upper triangular part of the array  A must contain the upper 
*/
/*           triangular matrix  and the strictly lower triangular part of 
*/
/*           A is not referenced. */
/*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k 
*/
/*           lower triangular part of the array  A must contain the lower 
*/
/*           triangular matrix  and the strictly upper triangular part of 
*/
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of 
*/
/*           A  are not referenced either,  but are assumed to be  unity. 
*/
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
*/
/*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' 
*/
/*           then LDA must be at least max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
/*           Before entry,  the leading  m by n part of the array  B must 
*/
/*           contain  the  right-hand  side  matrix  B,  and  on exit  is 
*/
/*           overwritten by the solution matrix  X. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared 
*/
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least 
*/
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */


/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    lside = lsame(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame(diag, "N");
    upper = lsame(uplo, "U");

    info = 0;
    if (! lside && ! lsame(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame(uplo, "L")) {
	info = 2;
    } else if (! lsame(transa, "N") && ! lsame(transa, "T") 
	    && ! lsame(transa, "C")) {
	info = 3;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla("DTRSM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (lsame(transa, "N")) {

/*           Form  B := alpha*inv( A )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L30: */
			}
		    }
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != 0.) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__2 = k - 1;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L40: */
			    }
			}
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L70: */
			}
		    }
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != 0.) {
			    if (nounit) {
				b[k + j * b_dim1] /= a[k + k * a_dim1];
			    }
			    i__3 = *m;
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
/* L80: */
			    }
			}
/* L90: */
		    }
/* L100: */
		}
	    }
	} else {

/*           Form  B := alpha*inv( A' )*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__3 = i__ - 1;
			for (k = 1; k <= i__3; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L110: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L120: */
		    }
/* L130: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = *alpha * b[i__ + j * b_dim1];
			i__2 = *m;
			for (k = i__ + 1; k <= i__2; ++k) {
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L140: */
			}
			if (nounit) {
			    temp /= a[i__ + i__ * a_dim1];
			}
			b[i__ + j * b_dim1] = temp;
/* L150: */
		    }
/* L160: */
		}
	    }
	}
    } else {
	if (lsame(transa, "N")) {

/*           Form  B := alpha*B*inv( A ). */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L170: */
			}
		    }
		    i__2 = j - 1;
		    for (k = 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L180: */
			    }
			}
/* L190: */
		    }
		    if (nounit) {
			temp = 1. / a[j + j * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L200: */
			}
		    }
/* L210: */
		}
	    } else {
		for (j = *n; j >= 1; --j) {
		    if (*alpha != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
/* L220: */
			}
		    }
		    i__1 = *n;
		    for (k = j + 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    if (nounit) {
			temp = 1. / a[j + j * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    }
	} else {

/*           Form  B := alpha*B*inv( A' ). */

	    if (upper) {
		for (k = *n; k >= 1; --k) {
		    if (nounit) {
			temp = 1. / a[k + k * a_dim1];
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L270: */
			}
		    }
		    i__1 = k - 1;
		    for (j = 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = a[j + k * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L280: */
			    }
			}
/* L290: */
		    }
		    if (*alpha != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L300: */
			}
		    }
/* L310: */
		}
	    } else {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    if (nounit) {
			temp = 1. / a[k + k * a_dim1];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L320: */
			}
		    }
		    i__2 = *n;
		    for (j = k + 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = a[j + k * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
/* L330: */
			    }
			}
/* L340: */
		    }
		    if (*alpha != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
/* L350: */
			}
		    }
/* L360: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSM . */

} /* dtrsm_ */

/* Subroutine */ int dtrmm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, 
	ldb)
char *side, *uplo, *transa, *diag;
integer *m, *n;
doublereal *alpha, *a;
integer *lda;
doublereal *b;
integer *ldb;
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
    static logical lside;
    static integer nrowa;
    static logical upper;
    static logical nounit;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRMM  performs one of the matrix-matrix operations */

/*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ), */

/*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or 
*/
/*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of 
*/

/*     op( A ) = A   or   op( A ) = A'. */

/*  Parameters */
/*  ========== */

/*  SIDE   - CHARACTER*1. */
/*           On entry,  SIDE specifies whether  op( A ) multiplies B from 
*/
/*           the left or right as follows: */

/*              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */

/*              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */

/*           Unchanged on exit. */

/*  UPLO   - CHARACTER*1. */
/*           On entry, UPLO specifies whether the matrix A is an upper or 
*/
/*           lower triangular matrix as follows: */

/*              UPLO = 'U' or 'u'   A is an upper triangular matrix. */

/*              UPLO = 'L' or 'l'   A is a lower triangular matrix. */

/*           Unchanged on exit. */

/*  TRANSA - CHARACTER*1. */
/*           On entry, TRANSA specifies the form of op( A ) to be used in 
*/
/*           the matrix multiplication as follows: */

/*              TRANSA = 'N' or 'n'   op( A ) = A. */

/*              TRANSA = 'T' or 't'   op( A ) = A'. */

/*              TRANSA = 'C' or 'c'   op( A ) = A'. */

/*           Unchanged on exit. */

/*  DIAG   - CHARACTER*1. */
/*           On entry, DIAG specifies whether or not A is unit triangular 
*/
/*           as follows: */

/*              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */

/*              DIAG = 'N' or 'n'   A is not assumed to be unit */
/*                                  triangular. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of B. M must be at 
*/
/*           least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of B.  N must be 
*/
/*           at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is 
*/
/*           zero then  A is not referenced and  B need not be set before 
*/
/*           entry. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m 
*/
/*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. 
*/
/*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k 
*/
/*           upper triangular part of the array  A must contain the upper 
*/
/*           triangular matrix  and the strictly lower triangular part of 
*/
/*           A is not referenced. */
/*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k 
*/
/*           lower triangular part of the array  A must contain the lower 
*/
/*           triangular matrix  and the strictly upper triangular part of 
*/
/*           A is not referenced. */
/*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of 
*/
/*           A  are not referenced either,  but are assumed to be  unity. 
*/
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
*/
/*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' 
*/
/*           then LDA must be at least max( 1, n ). */
/*           Unchanged on exit. */

/*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
/*           Before entry,  the leading  m by n part of the array  B must 
*/
/*           contain the matrix  B,  and  on exit  is overwritten  by the 
*/
/*           transformed matrix. */

/*  LDB    - INTEGER. */
/*           On entry, LDB specifies the first dimension of B as declared 
*/
/*           in  the  calling  (sub)  program.   LDB  must  be  at  least 
*/
/*           max( 1, m ). */
/*           Unchanged on exit. */


/*  Level 3 Blas routine. */

/*  -- Written on 8-February-1989. */
/*     Jack Dongarra, Argonne National Laboratory. */
/*     Iain Duff, AERE Harwell. */
/*     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/*     Sven Hammarling, Numerical Algorithms Group Ltd. */


/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Local Scalars .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = b_dim1 + 1;
    b -= b_offset;

    /* Function Body */
    lside = lsame(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame(diag, "N");
    upper = lsame(uplo, "U");

    info = 0;
    if (! lside && ! lsame(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame(uplo, "L")) {
	info = 2;
    } else if (! lsame(transa, "N") && ! lsame(transa, "T") 
	    && ! lsame(transa, "C")) {
	info = 3;
    } else if (! lsame(diag, "U") && ! lsame(diag, "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla("DTRMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (lsame(transa, "N")) {

/*           Form  B := alpha*A*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b[k + j * b_dim1] != 0.) {
			    temp = *alpha * b[k + j * b_dim1];
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
/* L30: */
			    }
			    if (nounit) {
				temp *= a[k + k * a_dim1];
			    }
			    b[k + j * b_dim1] = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (b[k + j * b_dim1] != 0.) {
			    temp = *alpha * b[k + j * b_dim1];
			    b[k + j * b_dim1] = temp;
			    if (nounit) {
				b[k + j * b_dim1] *= a[k + k * a_dim1];
			    }
			    i__2 = *m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {

/*           Form  B := alpha*A'*B. */

	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L90: */
			}
			b[i__ + j * b_dim1] = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b[i__ + j * b_dim1];
			if (nounit) {
			    temp *= a[i__ + i__ * a_dim1];
			}
			i__3 = *m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
/* L120: */
			}
			b[i__ + j * b_dim1] = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (lsame(transa, "N")) {

/*           Form  B := alpha*B*A. */

	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L150: */
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    temp = *alpha * a[k + j * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a[j + j * a_dim1];
		    }
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
/* L190: */
		    }
		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (a[k + j * a_dim1] != 0.) {
			    temp = *alpha * a[k + j * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {

/*           Form  B := alpha*B*A'. */

	    if (upper) {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = *alpha * a[j + k * a_dim1];
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (temp != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (a[j + k * a_dim1] != 0.) {
			    temp = *alpha * a[j + k * a_dim1];
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a[k + k * a_dim1];
		    }
		    if (temp != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRMM . */

} /* dtrmm_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
