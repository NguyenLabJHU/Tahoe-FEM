/* $Id: dpotrf.c,v 1.1 2005-01-03 00:46:56 paklein Exp $ */
/* dpotrf.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "pspases_f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b13 = -1.;
static doublereal c_b14 = 1.;

/*<       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO ) >*/
/* Subroutine */ int dpotrf_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer j, jb, nb;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dpotf2_(char *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     March 31, 1993 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER          UPLO >*/
/*<       INTEGER            INFO, LDA, N >*/
/*     .. */
/*     .. Array Arguments .. */
/*<       DOUBLE PRECISION   A( LDA, * ) >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DPOTRF computes the Cholesky factorization of a real symmetric */
/*  positive definite matrix A. */

/*  The factorization has the form */
/*     A = U**T * U,  if UPLO = 'U', or */
/*     A = L  * L**T,  if UPLO = 'L', */
/*  where U is an upper triangular matrix and L is lower triangular. */

/*  This is the block version of the algorithm, calling Level 3 BLAS. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          = 'U':  Upper triangle of A is stored; */
/*          = 'L':  Lower triangle of A is stored. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/*          N-by-N upper triangular part of A contains the upper */
/*          triangular part of the matrix A, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading N-by-N lower triangular part of A contains the lower */
/*          triangular part of the matrix A, and the strictly upper */
/*          triangular part of A is not referenced. */

/*          On exit, if INFO = 0, the factor U or L from the Cholesky */
/*          factorization A = U**T*U or A = L*L**T. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
/*          > 0:  if INFO = i, the leading minor of order i is not */
/*                positive definite, and the factorization could not be */
/*                completed. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*<       DOUBLE PRECISION   ONE >*/
/*<       PARAMETER          ( ONE = 1.0D+0 ) >*/
/*     .. */
/*     .. Local Scalars .. */
/*<       LOGICAL            UPPER >*/
/*<       INTEGER            J, JB, NB >*/
/*     .. */
/*     .. External Functions .. */
/*<       LOGICAL            LSAME >*/
/*<       INTEGER            ILAENV >*/
/*<       EXTERNAL           LSAME, ILAENV >*/
/*     .. */
/*     .. External Subroutines .. */
/*<       EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA >*/
/*     .. */
/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          MAX, MIN >*/
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

/*<       INFO = 0 >*/
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
/*<       UPPER = LSAME( UPLO, 'U' ) >*/
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
/*<       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN >*/
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
/*<          INFO = -1 >*/
	*info = -1;
/*<       ELSE IF( N.LT.0 ) THEN >*/
    } else if (*n < 0) {
/*<          INFO = -2 >*/
	*info = -2;
/*<       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN >*/
    } else if (*lda < max(1,*n)) {
/*<          INFO = -4 >*/
	*info = -4;
/*<       END IF >*/
    }
/*<       IF( INFO.NE.0 ) THEN >*/
    if (*info != 0) {
/*<          CALL XERBLA( 'DPOTRF', -INFO ) >*/
	i__1 = -(*info);
	xerbla_("DPOTRF", &i__1, (ftnlen)6);
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }

/*     Quick return if possible */

/*<    >*/
    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

/*<       NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 ) >*/
    nb = ilaenv_(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
/*<       IF( NB.LE.1 .OR. NB.GE.N ) THEN >*/
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

/*<          CALL DPOTF2( UPLO, N, A, LDA, INFO ) >*/
	dpotf2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
/*<       ELSE >*/
    } else {

/*        Use blocked code. */

/*<          IF( UPPER ) THEN >*/
	if (upper) {

/*           Compute the Cholesky factorization A = U'*U. */

/*<             DO 10 J = 1, N, NB >*/
	    i__1 = *n;
	    i__2 = nb;
	    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/*<                JB = MIN( NB, N-J+1 ) >*/
/* Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = min(i__3,i__4);
/*<    >*/
		i__3 = j - 1;
		dsyrk_("Upper", "Transpose", &jb, &i__3, &c_b13, &a[j * 
			a_dim1 + 1], lda, &c_b14, &a[j + j * a_dim1], lda, (
			ftnlen)5, (ftnlen)9);
/*<                CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO ) >*/
		dpotf2_("Upper", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
			5);
/*<    >*/
		if (*info != 0) {
		    goto L30;
		}
/*<                IF( J+JB.LE.N ) THEN >*/
		if (j + jb <= *n) {

/*                 Compute the current block row. */

/*<    >*/
		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    dgemm_("Transpose", "No transpose", &jb, &i__3, &i__4, &
			    c_b13, &a[j * a_dim1 + 1], lda, &a[(j + jb) * 
			    a_dim1 + 1], lda, &c_b14, &a[j + (j + jb) * 
			    a_dim1], lda, (ftnlen)9, (ftnlen)12);
/*<    >*/
		    i__3 = *n - j - jb + 1;
		    dtrsm_("Left", "Upper", "Transpose", "Non-unit", &jb, &
			    i__3, &c_b14, &a[j + j * a_dim1], lda, &a[j + (j 
			    + jb) * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
			    ftnlen)9, (ftnlen)8);
/*<                END IF >*/
		}
/*<    10       CONTINUE >*/
/* L10: */
	    }

/*<          ELSE >*/
	} else {

/*           Compute the Cholesky factorization A = L*L'. */

/*<             DO 20 J = 1, N, NB >*/
	    i__2 = *n;
	    i__1 = nb;
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/*<                JB = MIN( NB, N-J+1 ) >*/
/* Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = min(i__3,i__4);
/*<    >*/
		i__3 = j - 1;
		dsyrk_("Lower", "No transpose", &jb, &i__3, &c_b13, &a[j + 
			a_dim1], lda, &c_b14, &a[j + j * a_dim1], lda, (
			ftnlen)5, (ftnlen)12);
/*<                CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO ) >*/
		dpotf2_("Lower", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
			5);
/*<    >*/
		if (*info != 0) {
		    goto L30;
		}
/*<                IF( J+JB.LE.N ) THEN >*/
		if (j + jb <= *n) {

/*                 Compute the current block column. */

/*<    >*/
		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    dgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &
			    c_b13, &a[j + jb + a_dim1], lda, &a[j + a_dim1], 
			    lda, &c_b14, &a[j + jb + j * a_dim1], lda, (
			    ftnlen)12, (ftnlen)9);
/*<    >*/
		    i__3 = *n - j - jb + 1;
		    dtrsm_("Right", "Lower", "Transpose", "Non-unit", &i__3, &
			    jb, &c_b14, &a[j + j * a_dim1], lda, &a[j + jb + 
			    j * a_dim1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)9,
			     (ftnlen)8);
/*<                END IF >*/
		}
/*<    20       CONTINUE >*/
/* L20: */
	    }
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       GO TO 40 >*/
    goto L40;

/*<    30 CONTINUE >*/
L30:
/*<       INFO = INFO + J - 1 >*/
    *info = *info + j - 1;

/*<    40 CONTINUE >*/
L40:
/*<       RETURN >*/
    return 0;

/*     End of DPOTRF */

/*<       END >*/
} /* dpotrf_ */

/*<       SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO ) >*/
/* Subroutine */ int dpotf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer j;
    doublereal ajj;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     February 29, 1992 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER          UPLO >*/
/*<       INTEGER            INFO, LDA, N >*/
/*     .. */
/*     .. Array Arguments .. */
/*<       DOUBLE PRECISION   A( LDA, * ) >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  DPOTF2 computes the Cholesky factorization of a real symmetric */
/*  positive definite matrix A. */

/*  The factorization has the form */
/*     A = U' * U ,  if UPLO = 'U', or */
/*     A = L  * L',  if UPLO = 'L', */
/*  where U is an upper triangular matrix and L is lower triangular. */

/*  This is the unblocked version of the algorithm, calling Level 2 BLAS. */

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
/*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/*          n by n upper triangular part of A contains the upper */
/*          triangular part of the matrix A, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading n by n lower triangular part of A contains the lower */
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

/*  ===================================================================== */

/*     .. Parameters .. */
/*<       DOUBLE PRECISION   ONE, ZERO >*/
/*<       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 ) >*/
/*     .. */
/*     .. Local Scalars .. */
/*<       LOGICAL            UPPER >*/
/*<       INTEGER            J >*/
/*<       DOUBLE PRECISION   AJJ >*/
/*     .. */
/*     .. External Functions .. */
/*<       LOGICAL            LSAME >*/
/*<       DOUBLE PRECISION   DDOT >*/
/*<       EXTERNAL           LSAME, DDOT >*/
/*     .. */
/*     .. External Subroutines .. */
/*<       EXTERNAL           DGEMV, DSCAL, XERBLA >*/
/*     .. */
/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          MAX, SQRT >*/
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

/*<       INFO = 0 >*/
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
/*<       UPPER = LSAME( UPLO, 'U' ) >*/
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
/*<       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN >*/
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
/*<          INFO = -1 >*/
	*info = -1;
/*<       ELSE IF( N.LT.0 ) THEN >*/
    } else if (*n < 0) {
/*<          INFO = -2 >*/
	*info = -2;
/*<       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN >*/
    } else if (*lda < max(1,*n)) {
/*<          INFO = -4 >*/
	*info = -4;
/*<       END IF >*/
    }
/*<       IF( INFO.NE.0 ) THEN >*/
    if (*info != 0) {
/*<          CALL XERBLA( 'DPOTF2', -INFO ) >*/
	i__1 = -(*info);
	xerbla_("DPOTF2", &i__1, (ftnlen)6);
/*<          RETURN >*/
	return 0;
/*<       END IF >*/
    }

/*     Quick return if possible */

/*<    >*/
    if (*n == 0) {
	return 0;
    }

/*<       IF( UPPER ) THEN >*/
    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

/*<          DO 10 J = 1, N >*/
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness. */

/*<             AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 ) >*/
	    i__2 = j - 1;
	    ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j * a_dim1 + 1], &c__1, 
		    &a[j * a_dim1 + 1], &c__1);
/*<             IF( AJJ.LE.ZERO ) THEN >*/
	    if (ajj <= 0.) {
/*<                A( J, J ) = AJJ >*/
		a[j + j * a_dim1] = ajj;
/*<                GO TO 30 >*/
		goto L30;
/*<             END IF >*/
	    }
/*<             AJJ = SQRT( AJJ ) >*/
	    ajj = sqrt(ajj);
/*<             A( J, J ) = AJJ >*/
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of row J. */

/*<             IF( J.LT.N ) THEN >*/
	    if (j < *n) {
/*<    >*/
		i__2 = j - 1;
		i__3 = *n - j;
		dgemv_("Transpose", &i__2, &i__3, &c_b13, &a[(j + 1) * a_dim1 
			+ 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b14, &a[j + (
			j + 1) * a_dim1], lda, (ftnlen)9);
/*<                CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA ) >*/
		i__2 = *n - j;
		d__1 = 1. / ajj;
		dscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
/*<             END IF >*/
	    }
/*<    10    CONTINUE >*/
/* L10: */
	}
/*<       ELSE >*/
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

/*<          DO 20 J = 1, N >*/
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness. */

/*<    >*/
	    i__2 = j - 1;
	    ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j + a_dim1], lda, &a[j 
		    + a_dim1], lda);
/*<             IF( AJJ.LE.ZERO ) THEN >*/
	    if (ajj <= 0.) {
/*<                A( J, J ) = AJJ >*/
		a[j + j * a_dim1] = ajj;
/*<                GO TO 30 >*/
		goto L30;
/*<             END IF >*/
	    }
/*<             AJJ = SQRT( AJJ ) >*/
	    ajj = sqrt(ajj);
/*<             A( J, J ) = AJJ >*/
	    a[j + j * a_dim1] = ajj;

/*           Compute elements J+1:N of column J. */

/*<             IF( J.LT.N ) THEN >*/
	    if (j < *n) {
/*<    >*/
		i__2 = *n - j;
		i__3 = j - 1;
		dgemv_("No transpose", &i__2, &i__3, &c_b13, &a[j + 1 + 
			a_dim1], lda, &a[j + a_dim1], lda, &c_b14, &a[j + 1 + 
			j * a_dim1], &c__1, (ftnlen)12);
/*<                CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 ) >*/
		i__2 = *n - j;
		d__1 = 1. / ajj;
		dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
/*<             END IF >*/
	    }
/*<    20    CONTINUE >*/
/* L20: */
	}
/*<       END IF >*/
    }
/*<       GO TO 40 >*/
    goto L40;

/*<    30 CONTINUE >*/
L30:
/*<       INFO = J >*/
    *info = j;

/*<    40 CONTINUE >*/
L40:
/*<       RETURN >*/
    return 0;

/*     End of DPOTF2 */

/*<       END >*/
} /* dpotf2_ */

/*<    >*/
integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__;
    char c1[1], c2[2], c3[3], c4[2];
    integer ic, nb, iz, nx;
    logical cname, sname;
    integer nbmin;
    char subnam[6];


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER*( * )    NAME, OPTS >*/
/*<       INTEGER            ISPEC, N1, N2, N3, N4 >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV is called from the LAPACK routines to choose problem-dependent */
/*  parameters for the local environment.  See ISPEC for a description of */
/*  the parameters. */

/*  This version provides a set of parameters which should give good, */
/*  but not optimal, performance on many of the currently available */
/*  computers.  Users are encouraged to modify this subroutine to set */
/*  the tuning parameters for their particular machine using the option */
/*  and problem size information in the arguments. */

/*  This routine will not function correctly if it is converted to all */
/*  lower case.  Converting it to all upper case is allowed. */

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies the parameter to be returned as the value of */
/*          ILAENV. */
/*          = 1: the optimal blocksize; if this value is 1, an unblocked */
/*               algorithm will give the best performance. */
/*          = 2: the minimum block size for which the block routine */
/*               should be used; if the usable block size is less than */
/*               this value, an unblocked routine should be used. */
/*          = 3: the crossover point (in a block routine, for N less */
/*               than this value, an unblocked routine should be used) */
/*          = 4: the number of shifts, used in the nonsymmetric */
/*               eigenvalue routines */
/*          = 5: the minimum column dimension for blocking to be used; */
/*               rectangular blocks must have dimension at least k by m, */
/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
/*          = 6: the crossover point for the SVD (when reducing an m by n */
/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
/*               this value, a QR factorization is used first to reduce */
/*               the matrix to a triangular form.) */
/*          = 7: the number of processors */
/*          = 8: the crossover point for the multishift QR and QZ methods */
/*               for nonsymmetric eigenvalue problems. */

/*  NAME    (input) CHARACTER*(*) */
/*          The name of the calling subroutine, in either upper case or */
/*          lower case. */

/*  OPTS    (input) CHARACTER*(*) */
/*          The character options to the subroutine NAME, concatenated */
/*          into a single character string.  For example, UPLO = 'U', */
/*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/*          be specified as OPTS = 'UTN'. */

/*  N1      (input) INTEGER */
/*  N2      (input) INTEGER */
/*  N3      (input) INTEGER */
/*  N4      (input) INTEGER */
/*          Problem dimensions for the subroutine NAME; these may not all */
/*          be required. */

/* (ILAENV) (output) INTEGER */
/*          >= 0: the value of the parameter specified by ISPEC */
/*          < 0:  if ILAENV = -k, the k-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  The following conventions have been used when calling ILAENV from the */
/*  LAPACK routines: */
/*  1)  OPTS is a concatenation of all of the character options to */
/*      subroutine NAME, in the same order that they appear in the */
/*      argument list for NAME, even if they are not used in determining */
/*      the value of the parameter specified by ISPEC. */
/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
/*      that they appear in the argument list for NAME.  N1 is used */
/*      first, N2 second, and so on, and unused problem dimensions are */
/*      passed a value of -1. */
/*  3)  The parameter value returned by ILAENV is checked for validity in */
/*      the calling subroutine.  For example, ILAENV is used to retrieve */
/*      the optimal blocksize for STRTRI as follows: */

/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*<       LOGICAL            CNAME, SNAME >*/
/*<       CHARACTER*1        C1 >*/
/*<       CHARACTER*2        C2, C4 >*/
/*<       CHARACTER*3        C3 >*/
/*<       CHARACTER*6        SUBNAM >*/
/*<       INTEGER            I, IC, IZ, NB, NBMIN, NX >*/
/*     .. */
/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          CHAR, ICHAR, INT, MIN, REAL >*/
/*     .. */
/*     .. Executable Statements .. */

/*<       GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC >*/
    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
    }

/*     Invalid value for ISPEC */

/*<       ILAENV = -1 >*/
    ret_val = -1;
/*<       RETURN >*/
    return ret_val;

/*<   100 CONTINUE >*/
L100:

/*     Convert NAME to upper case if the first character is lower case. */

/*<       ILAENV = 1 >*/
    ret_val = 1;
/*<       SUBNAM = NAME >*/
    s_copy(subnam, name__, (ftnlen)6, name_len);
/*<       IC = ICHAR( SUBNAM( 1:1 ) ) >*/
    ic = *(unsigned char *)subnam;
/*<       IZ = ICHAR( 'Z' ) >*/
    iz = 'Z';
/*<       IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN >*/
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

/*<          IF( IC.GE.97 .AND. IC.LE.122 ) THEN >*/
	if (ic >= 97 && ic <= 122) {
/*<             SUBNAM( 1:1 ) = CHAR( IC-32 ) >*/
	    *(unsigned char *)subnam = (char) (ic - 32);
/*<             DO 10 I = 2, 6 >*/
	    for (i__ = 2; i__ <= 6; ++i__) {
/*<                IC = ICHAR( SUBNAM( I:I ) ) >*/
		ic = *(unsigned char *)&subnam[i__ - 1];
/*<    >*/
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/*<    10       CONTINUE >*/
/* L10: */
	    }
/*<          END IF >*/
	}

/*<       ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN >*/
    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

/*<    >*/
	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
/*<             SUBNAM( 1:1 ) = CHAR( IC+64 ) >*/
	    *(unsigned char *)subnam = (char) (ic + 64);
/*<             DO 20 I = 2, 6 >*/
	    for (i__ = 2; i__ <= 6; ++i__) {
/*<                IC = ICHAR( SUBNAM( I:I ) ) >*/
		ic = *(unsigned char *)&subnam[i__ - 1];
/*<    >*/
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		}
/*<    20       CONTINUE >*/
/* L20: */
	    }
/*<          END IF >*/
	}

/*<       ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN >*/
    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

/*<          IF( IC.GE.225 .AND. IC.LE.250 ) THEN >*/
	if (ic >= 225 && ic <= 250) {
/*<             SUBNAM( 1:1 ) = CHAR( IC-32 ) >*/
	    *(unsigned char *)subnam = (char) (ic - 32);
/*<             DO 30 I = 2, 6 >*/
	    for (i__ = 2; i__ <= 6; ++i__) {
/*<                IC = ICHAR( SUBNAM( I:I ) ) >*/
		ic = *(unsigned char *)&subnam[i__ - 1];
/*<    >*/
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/*<    30       CONTINUE >*/
/* L30: */
	    }
/*<          END IF >*/
	}
/*<       END IF >*/
    }

/*<       C1 = SUBNAM( 1:1 ) >*/
    *(unsigned char *)c1 = *(unsigned char *)subnam;
/*<       SNAME = C1.EQ.'S' .OR. C1.EQ.'D' >*/
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
/*<       CNAME = C1.EQ.'C' .OR. C1.EQ.'Z' >*/
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
/*<    >*/
    if (! (cname || sname)) {
	return ret_val;
    }
/*<       C2 = SUBNAM( 2:3 ) >*/
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
/*<       C3 = SUBNAM( 4:6 ) >*/
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
/*<       C4 = C3( 2:3 ) >*/
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

/*<       GO TO ( 110, 200, 300 ) ISPEC >*/
    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

/*<   110 CONTINUE >*/
L110:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

/*<       NB = 1 >*/
    nb = 1;

/*<       IF( C2.EQ.'GE' ) THEN >*/
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 64 >*/
		nb = 64;
/*<             ELSE >*/
	    } else {
/*<                NB = 64 >*/
		nb = 64;
/*<             END IF >*/
	    }
/*<    >*/
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 32 >*/
		nb = 32;
/*<             ELSE >*/
	    } else {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'HRD' ) THEN >*/
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 32 >*/
		nb = 32;
/*<             ELSE >*/
	    } else {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'BRD' ) THEN >*/
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 32 >*/
		nb = 32;
/*<             ELSE >*/
	    } else {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'TRI' ) THEN >*/
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 64 >*/
		nb = 64;
/*<             ELSE >*/
	    } else {
/*<                NB = 64 >*/
		nb = 64;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'PO' ) THEN >*/
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 64 >*/
		nb = 64;
/*<             ELSE >*/
	    } else {
/*<                NB = 64 >*/
		nb = 64;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'SY' ) THEN >*/
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 64 >*/
		nb = 64;
/*<             ELSE >*/
	    } else {
/*<                NB = 64 >*/
		nb = 64;
/*<             END IF >*/
	    }
/*<          ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN >*/
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NB = 1 >*/
	    nb = 1;
/*<          ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN >*/
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NB = 64 >*/
	    nb = 64;
/*<          END IF >*/
	}
/*<       ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN >*/
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NB = 64 >*/
	    nb = 64;
/*<          ELSE IF( C3.EQ.'TRD' ) THEN >*/
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NB = 1 >*/
	    nb = 1;
/*<          ELSE IF( C3.EQ.'GST' ) THEN >*/
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NB = 64 >*/
	    nb = 64;
/*<          END IF >*/
	}
/*<       ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN >*/
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3( 1:1 ).EQ.'G' ) THEN >*/
	if (*(unsigned char *)c3 == 'G') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3( 1:1 ).EQ.'M' ) THEN >*/
	} else if (*(unsigned char *)c3 == 'M') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN >*/
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3( 1:1 ).EQ.'G' ) THEN >*/
	if (*(unsigned char *)c3 == 'G') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3( 1:1 ).EQ.'M' ) THEN >*/
	} else if (*(unsigned char *)c3 == 'M') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NB = 32 >*/
		nb = 32;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'GB' ) THEN >*/
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                IF( N4.LE.64 ) THEN >*/
		if (*n4 <= 64) {
/*<                   NB = 1 >*/
		    nb = 1;
/*<                ELSE >*/
		} else {
/*<                   NB = 32 >*/
		    nb = 32;
/*<                END IF >*/
		}
/*<             ELSE >*/
	    } else {
/*<                IF( N4.LE.64 ) THEN >*/
		if (*n4 <= 64) {
/*<                   NB = 1 >*/
		    nb = 1;
/*<                ELSE >*/
		} else {
/*<                   NB = 32 >*/
		    nb = 32;
/*<                END IF >*/
		}
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'PB' ) THEN >*/
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                IF( N2.LE.64 ) THEN >*/
		if (*n2 <= 64) {
/*<                   NB = 1 >*/
		    nb = 1;
/*<                ELSE >*/
		} else {
/*<                   NB = 32 >*/
		    nb = 32;
/*<                END IF >*/
		}
/*<             ELSE >*/
	    } else {
/*<                IF( N2.LE.64 ) THEN >*/
		if (*n2 <= 64) {
/*<                   NB = 1 >*/
		    nb = 1;
/*<                ELSE >*/
		} else {
/*<                   NB = 32 >*/
		    nb = 32;
/*<                END IF >*/
		}
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'TR' ) THEN >*/
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRI' ) THEN >*/
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 64 >*/
		nb = 64;
/*<             ELSE >*/
	    } else {
/*<                NB = 64 >*/
		nb = 64;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'LA' ) THEN >*/
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'UUM' ) THEN >*/
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NB = 64 >*/
		nb = 64;
/*<             ELSE >*/
	    } else {
/*<                NB = 64 >*/
		nb = 64;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN >*/
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'EBZ' ) THEN >*/
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NB = 1 >*/
	    nb = 1;
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       ILAENV = NB >*/
    ret_val = nb;
/*<       RETURN >*/
    return ret_val;

/*<   200 CONTINUE >*/
L200:

/*     ISPEC = 2:  minimum block size */

/*<       NBMIN = 2 >*/
    nbmin = 2;
/*<       IF( C2.EQ.'GE' ) THEN >*/
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
/*<    >*/
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             ELSE >*/
	    } else {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'HRD' ) THEN >*/
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             ELSE >*/
	    } else {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'BRD' ) THEN >*/
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             ELSE >*/
	    } else {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'TRI' ) THEN >*/
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             ELSE >*/
	    } else {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'SY' ) THEN >*/
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRF' ) THEN >*/
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NBMIN = 8 >*/
		nbmin = 8;
/*<             ELSE >*/
	    } else {
/*<                NBMIN = 8 >*/
		nbmin = 8;
/*<             END IF >*/
	    }
/*<          ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN >*/
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NBMIN = 2 >*/
	    nbmin = 2;
/*<          END IF >*/
	}
/*<       ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN >*/
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRD' ) THEN >*/
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NBMIN = 2 >*/
	    nbmin = 2;
/*<          END IF >*/
	}
/*<       ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN >*/
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3( 1:1 ).EQ.'G' ) THEN >*/
	if (*(unsigned char *)c3 == 'G') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3( 1:1 ).EQ.'M' ) THEN >*/
	} else if (*(unsigned char *)c3 == 'M') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN >*/
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3( 1:1 ).EQ.'G' ) THEN >*/
	if (*(unsigned char *)c3 == 'G') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3( 1:1 ).EQ.'M' ) THEN >*/
	} else if (*(unsigned char *)c3 == 'M') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NBMIN = 2 >*/
		nbmin = 2;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       ILAENV = NBMIN >*/
    ret_val = nbmin;
/*<       RETURN >*/
    return ret_val;

/*<   300 CONTINUE >*/
L300:

/*     ISPEC = 3:  crossover point */

/*<       NX = 0 >*/
    nx = 0;
/*<       IF( C2.EQ.'GE' ) THEN >*/
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
/*<    >*/
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NX = 128 >*/
		nx = 128;
/*<             ELSE >*/
	    } else {
/*<                NX = 128 >*/
		nx = 128;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'HRD' ) THEN >*/
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NX = 128 >*/
		nx = 128;
/*<             ELSE >*/
	    } else {
/*<                NX = 128 >*/
		nx = 128;
/*<             END IF >*/
	    }
/*<          ELSE IF( C3.EQ.'BRD' ) THEN >*/
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             IF( SNAME ) THEN >*/
	    if (sname) {
/*<                NX = 128 >*/
		nx = 128;
/*<             ELSE >*/
	    } else {
/*<                NX = 128 >*/
		nx = 128;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( C2.EQ.'SY' ) THEN >*/
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( SNAME .AND. C3.EQ.'TRD' ) THEN >*/
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NX = 1 >*/
	    nx = 1;
/*<          END IF >*/
	}
/*<       ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN >*/
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3.EQ.'TRD' ) THEN >*/
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
/*<             NX = 1 >*/
	    nx = 1;
/*<          END IF >*/
	}
/*<       ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN >*/
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3( 1:1 ).EQ.'G' ) THEN >*/
	if (*(unsigned char *)c3 == 'G') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NX = 128 >*/
		nx = 128;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN >*/
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
/*<          IF( C3( 1:1 ).EQ.'G' ) THEN >*/
	if (*(unsigned char *)c3 == 'G') {
/*<    >*/
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
/*<                NX = 128 >*/
		nx = 128;
/*<             END IF >*/
	    }
/*<          END IF >*/
	}
/*<       END IF >*/
    }
/*<       ILAENV = NX >*/
    ret_val = nx;
/*<       RETURN >*/
    return ret_val;

/*<   400 CONTINUE >*/
L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

/*<       ILAENV = 6 >*/
    ret_val = 6;
/*<       RETURN >*/
    return ret_val;

/*<   500 CONTINUE >*/
L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

/*<       ILAENV = 2 >*/
    ret_val = 2;
/*<       RETURN >*/
    return ret_val;

/*<   600 CONTINUE  >*/
L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

/*<       ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 ) >*/
    ret_val = (integer) ((real) min(*n1,*n2) * 1.6f);
/*<       RETURN >*/
    return ret_val;

/*<   700 CONTINUE >*/
L700:

/*     ISPEC = 7:  number of processors (not used) */

/*<       ILAENV = 1 >*/
    ret_val = 1;
/*<       RETURN >*/
    return ret_val;

/*<   800 CONTINUE >*/
L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

/*<       ILAENV = 50 >*/
    ret_val = 50;
/*<       RETURN >*/
    return ret_val;

/*     End of ILAENV */

/*<       END >*/
} /* ilaenv_ */

/*<       LOGICAL          FUNCTION LSAME( CA, CB ) >*/
logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    integer inta, intb, zcode;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER          CA, CB >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  LSAME returns .TRUE. if CA is the same letter as CB regardless of */
/*  case. */

/*  Arguments */
/*  ========= */

/*  CA      (input) CHARACTER*1 */
/*  CB      (input) CHARACTER*1 */
/*          CA and CB specify the single characters to be compared. */

/* ===================================================================== */

/*     .. Intrinsic Functions .. */
/*<       INTRINSIC          ICHAR >*/
/*     .. */
/*     .. Local Scalars .. */
/*<       INTEGER            INTA, INTB, ZCODE >*/
/*     .. */
/*     .. Executable Statements .. */

/*     Test if the characters are equal */

/*<       LSAME = CA.EQ.CB >*/
    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
/*<    >*/
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

/*<       ZCODE = ICHAR( 'Z' ) >*/
    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
/*     machines, on which ICHAR returns a value with bit 8 set. */
/*     ICHAR('A') on Prime machines returns 193 which is the same as */
/*     ICHAR('A') on an EBCDIC machine. */

/*<       INTA = ICHAR( CA ) >*/
    inta = *(unsigned char *)ca;
/*<       INTB = ICHAR( CB ) >*/
    intb = *(unsigned char *)cb;

/*<       IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN >*/
    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
/*        upper case 'Z'. */

/*<          IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32 >*/
	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
/*<          IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32 >*/
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

/*<       ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN >*/
    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
/*        upper case 'Z'. */

/*<    >*/
	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
/*<    >*/
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

/*<       ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN >*/
    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
/*        plus 128 of either lower or upper case 'Z'. */

/*<          IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32 >*/
	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
/*<          IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32 >*/
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
/*<       END IF >*/
    }
/*<       LSAME = INTA.EQ.INTB >*/
    ret_val = inta == intb;

/*     RETURN */

/*     End of LSAME */

/*<       END >*/
    return ret_val;
} /* lsame_ */

/*<       SUBROUTINE XERBLA( SRNAME, INFO ) >*/
/* Subroutine */ int xerbla_(char *srname, integer *info, ftnlen srname_len)
{
	int i;
	printf("function \"");
	for (i = 0; i < srname_len; i++)
		printf("%c", srname[i]);
	printf("\" hit error code %d", *info);
	abort();
	return 0;

# if 0
    /* Format strings */
    static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu"
	    "mber \002,i2,\002 had \002,\002an illegal value\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*<       CHARACTER*6        SRNAME >*/
/*<       INTEGER            INFO >*/
/*     .. */

/*  Purpose */
/*  ======= */

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */

/*  Installers may consider modifying the STOP statement in order to */
/*  call system-specific exception-handling facilities. */

/*  Arguments */
/*  ========= */

/*  SRNAME  (input) CHARACTER*6 */
/*          The name of the routine which called XERBLA. */

/*  INFO    (input) INTEGER */
/*          The position of the invalid parameter in the parameter list */
/*          of the calling routine. */

/* ===================================================================== */

/*     .. Executable Statements .. */

/*<       WRITE( *, FMT = 9999 )SRNAME, INFO >*/
    s_wsfe(&io___24);
    do_fio(&c__1, srname, (ftnlen)6);
    do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
    e_wsfe();

/*<       STOP >*/
    s_stop("", (ftnlen)0);

/*<  9 >*/

/*     End of XERBLA */

/*<       END >*/
    return 0;
#endif
} /* xerbla_ */

