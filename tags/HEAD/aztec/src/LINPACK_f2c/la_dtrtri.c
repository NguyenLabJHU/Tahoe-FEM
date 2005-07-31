/* PAK (07/24/98) */

/* la_dtrtri.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "linpack.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

#include "utilsx.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b18 = 1.;
static doublereal c_b22 = -1.;

/* Subroutine */ int dtrtri(uplo, diag, n, a, lda, info)
char *uplo, *diag;
integer *n;
doublereal *a;
integer *lda, *info;
{
    /* System generated locals */
    char a__1[3];
    integer a_dim1, a_offset, i__1, i__3, i__4, i__5;

    /* Local variables */
    static integer j;
    static logical upper;
    static integer jb, nb, nn;
    static logical nounit;


/*  -- LAPACK routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  DTRTRI computes the inverse of a real upper or lower triangular */
/*  matrix A. */

/*  This is the Level 3 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  A is upper triangular; */
/*          = 'L':  A is lower triangular. */

/*  DIAG    (input) CHARACTER*1 */
/*          = 'N':  A is non-unit triangular; */
/*          = 'U':  A is unit triangular. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the triangular matrix A.  If UPLO = 'U', the */
/*          leading N-by-N upper triangular part of the array A contains 
*/
/*          the upper triangular matrix, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading N-by-N lower triangular part of the array A contains 
*/
/*          the lower triangular matrix, and the strictly upper */
/*          triangular part of A is not referenced.  If DIAG = 'U', the */
/*          diagonal elements of A are also not referenced and are */
/*          assumed to be 1. */
/*          On exit, the (triangular) inverse of the original matrix, in 
*/
/*          the same storage format. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */
/*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular */
/*               matrix is singular and its inverse can not be computed. 
*/

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
    nounit = lsame(diag, "N");
    if (! upper && ! lsame(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! lsame(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("DTRTRI", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (a[*info + *info * a_dim1] == 0.) {
		return 0;
	    }
/* L10: */
	}
	*info = 0;
    }

/*     Determine the block size for this environment. */

/* Writing concatenation */
	a__1[0] = *uplo;
	a__1[1] = *diag;
	a__1[2] = '\0';
    nb = ilaenv(&c__1, "DTRTRI", a__1, n, &c_n1, &c_n1, &c_n1);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

	dtrti2(uplo, diag, n, &a[a_offset], lda, info);
    } else {

/*        Use blocked code */

	if (upper) {

/*           Compute inverse of upper triangular matrix */

	    i__1 = *n;
	    i__3 = nb;
	    for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* Computing MIN */
		i__4 = nb, i__5 = *n - j + 1;
		jb = min(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

		i__4 = j - 1;
		dtrmm_("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &a[a_offset], lda, &a[j * a_dim1 + 1], lda);
		i__4 = j - 1;
		dtrsm_("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], lda);

/*              Compute inverse of current diagonal block */

		dtrti2("Upper", diag, &jb, &a[j + j * a_dim1], lda, info);
/* L20: */
	    }
	} else {

/*           Compute inverse of lower triangular matrix */

	    nn = (*n - 1) / nb * nb + 1;
	    i__3 = -nb;
	    for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
/* Computing MIN */
		i__1 = nb, i__4 = *n - j + 1;
		jb = min(i__1,i__4);
		if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block co
lumn */

		    i__1 = *n - j - jb + 1;
		    dtrmm_("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &a[j + jb + (j + jb) * a_dim1], lda, &a[j 
			    + jb + j * a_dim1], lda);
		    i__1 = *n - j - jb + 1;
		    dtrsm_("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &a[j + j * a_dim1], lda, &a[j + jb + j * 
			    a_dim1], lda);
		}

/*              Compute inverse of current diagonal block */

		dtrti2("Lower", diag, &jb, &a[j + j * a_dim1], lda, info);
/* L30: */
	    }
	}
    }

    return 0;

/*     End of DTRTRI */

} /* dtrtri_ */


#ifdef __MWERKS__
#pragma require_prototypes reset
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
