/* PAK (07/24/98) */

/* la_dpotf2.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>

#include "linpack.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

#include "utilsx.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b10 = -1.;
static doublereal c_b12 = 1.;

/* Subroutine */ int dpotf2(uplo, n, a, lda, info)
char *uplo;
integer *n;
doublereal *a;
integer *lda, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static logical upper;
    static doublereal ajj;


/*  -- LAPACK routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DPOTF2 computes the Cholesky factorization of a real symmetric */
/*  positive definite matrix A. */

/*  The factorization has the form */
/*     A = U' * U ,  if UPLO = 'U', or */
/*     A = L  * L',  if UPLO = 'L', */
/*  where U is an upper triangular matrix and L is lower triangular. */

/*  This is the unblocked version of the algorithm, calling Level 2 BLAS. 
*/

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          Specifies whether the upper or lower triangular part of the */
/*          symmetric matrix A is stored. */
/*          = 'U':  Upper triangular */
/*          = 'L':  Lower triangular */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
*/
/*          n by n upper triangular part of A contains the upper */
/*          triangular part of the matrix A, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading n by n lower triangular part of A contains the lower 
*/
/*          triangular part of the matrix A, and the strictly upper */
/*          triangular part of A is not referenced. */

/*          On exit, if INFO = 0, the factor U or L from the Cholesky */
/*          factorization A = U'*U  or A = L*L'. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -k, the k-th argument had an illegal value */
/*          > 0: if INFO = k, the leading minor of order k is not */
/*               positive definite, and the factorization could not be */
/*               completed. */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = lsame(uplo, "U");
    if (! upper && ! lsame(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("DPOTF2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;
	    ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j * a_dim1 + 1], &c__1, 
		    &a[j * a_dim1 + 1], &c__1);
	    if (ajj <= 0.) {
		a[j + j * a_dim1] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of row J. */

	    if (j < *n) {
		i__2 = j - 1;
		i__3 = *n - j;
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[(j + 1) * a_dim1 
			+ 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b12, &a[j + (
			j + 1) * a_dim1], lda);
		i__2 = *n - j;
		d__1 = 1. / ajj;
		dscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
	    }
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;
	    ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j + a_dim1], lda, &a[j 
		    + a_dim1], lda);
	    if (ajj <= 0.) {
		a[j + j * a_dim1] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of column J. */

	    if (j < *n) {
		i__2 = *n - j;
		i__3 = j - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[j + 1 + 
			a_dim1], lda, &a[j + a_dim1], lda, &c_b12, &a[j + 1 + 
			j * a_dim1], &c__1);
		i__2 = *n - j;
		d__1 = 1. / ajj;
		dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return 0;

/*     End of DPOTF2 */

} /* dpotf2_ */


#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
