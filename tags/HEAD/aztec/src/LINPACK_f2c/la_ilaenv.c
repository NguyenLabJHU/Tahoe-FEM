/* PAK (07/24/98) */

/* la_ilaenv.f -- translated by f2c (version 19960611).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <ctype.h>
#include <string.h>

#include "linpack.h"

#ifdef __MWERKS__
#pragma warn_unusedarg off
#endif /* __MWERKS__ */

#include "utilsx.h"

integer ilaenv(ispec, name__, opts, n1, n2, n3, n4)
integer *ispec;
char *name__, *opts;
integer *n1, *n2, *n3, *n4;
{
    /* System generated locals */
    integer ret_val;
    int sdex;

    /* Local variables */
    static integer i__;
    static logical cname, sname;
    static integer nbmin;
    static integer ic, nb, iz, nx;
    static char subnam[7];
    static char c1[2], c2[3], c3[4], c4[3];

	/* insert end of string markers */
    subnam[6] = '\0';
    c1[1] = '\0';
    c2[2] = '\0';
    c3[3] = '\0';
    c4[2] = '\0';


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV is called from the LAPACK routines to choose problem-dependent 
*/
/*  parameters for the local environment.  See ISPEC for a description of 
*/
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
/*          = 1: the optimal blocksize; if this value is 1, an unblocked 
*/
/*               algorithm will give the best performance. */
/*          = 2: the minimum block size for which the block routine */
/*               should be used; if the usable block size is less than */
/*               this value, an unblocked routine should be used. */
/*          = 3: the crossover point (in a block routine, for N less */
/*               than this value, an unblocked routine should be used) */
/*          = 4: the number of shifts, used in the nonsymmetric */
/*               eigenvalue routines */
/*          = 5: the minimum column dimension for blocking to be used; */
/*               rectangular blocks must have dimension at least k by m, 
*/
/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) 
*/
/*          = 6: the crossover point for the SVD (when reducing an m by n 
*/
/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds 
*/
/*               this value, a QR factorization is used first to reduce */
/*               the matrix to a triangular form.) */
/*          = 7: the number of processors */
/*          = 8: the crossover point for the multishift QR and QZ methods 
*/
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
/*          Problem dimensions for the subroutine NAME; these may not all 
*/
/*          be required. */

/* (ILAENV) (output) INTEGER */
/*          >= 0: the value of the parameter specified by ISPEC */
/*          < 0:  if ILAENV = -k, the k-th argument had an illegal value. 
*/

/*  Further Details */
/*  =============== */

/*  The following conventions have been used when calling ILAENV from the 
*/
/*  LAPACK routines: */
/*  1)  OPTS is a concatenation of all of the character options to */
/*      subroutine NAME, in the same order that they appear in the */
/*      argument list for NAME, even if they are not used in determining 
*/
/*      the value of the parameter specified by ISPEC. */
/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order 
*/
/*      that they appear in the argument list for NAME.  N1 is used */
/*      first, N2 second, and so on, and unused problem dimensions are */
/*      passed a value of -1. */
/*  3)  The parameter value returned by ILAENV is checked for validity in 
*/
/*      the calling subroutine.  For example, ILAENV is used to retrieve 
*/
/*      the optimal blocksize for STRTRI as follows: */

/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

/*  ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    switch ((int)*ispec) {
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

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    
    /* copy into subnam as CAPS */
    sdex = 0;
    while (name__[sdex] != '\0')
    {
    	subnam[sdex] = toupper(name__[sdex]);
    	sdex++;
    }

    c1[0] = subnam[0];
    sname = (c1[0] == 'S' || c1[0] == 'D');
    cname = (c1[0] == 'C' || c1[0] == 'Z');
        
    if (! (cname || sname)) {
	return ret_val;
    }
    
    /* get substrings */
    c2[0] = subnam[1]; c2[1] = subnam[2];
    c3[0] = subnam[3]; c3[0] = subnam[4]; c3[0] = subnam[5];
    c4[0] = c3[1]; c4[0] = c3[2]; 

    switch ((int)*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size */

/*     In these examples, separate code is provided for setting NB for */
/*     real and complex.  We assume that NB will take the same value in */
/*     single or double precision. */

    nb = 1;

    if (strcmp(c2, "GE") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (strcmp(c3, "QRF") == 0 || strcmp(c3, "RQF") 
		== 0 || strcmp(c3, "LQF") == 0 || strcmp(c3, "QLF") == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (strcmp(c3, "HRD") == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (strcmp(c3, "BRD") == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (strcmp(c3, "TRI") == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (strcmp(c2, "PO") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (strcmp(c2, "SY") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && strcmp(c3, "TRD") == 0) {
	    nb = 1;
	} else if (sname && strcmp(c3, "GST") == 0) {
	    nb = 64;
	}
    } else if (cname && strcmp(c2, "HE") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    nb = 64;
	} else if (strcmp(c3, "TRD") == 0) {
	    nb = 1;
	} else if (strcmp(c3, "GST") == 0) {
	    nb = 64;
	}
    } else if (sname && strcmp(c2, "OR") == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nb = 32;
	    }
	}
    } else if (cname && strcmp(c2, "UN") == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nb = 32;
	    }
	}
    } else if (strcmp(c2, "GB") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (strcmp(c2, "PB") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (strcmp(c2, "TR") == 0) {
	if (strcmp(c3, "TRI") == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (strcmp(c2, "LA") == 0) {
	if (strcmp(c3, "UUM") == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && strcmp(c2, "ST") == 0) {
	if (strcmp(c3, "EBZ") == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (strcmp(c2, "GE") == 0) {
	if (strcmp(c3, "QRF") == 0 || strcmp(c3, "RQF") == 0 || 
		strcmp(c3, "LQF") == 0 || strcmp(c3, "QLF") == 
		0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (strcmp(c3, "HRD") == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (strcmp(c3, "BRD") == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (strcmp(c3, "TRI") == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (strcmp(c2, "SY") == 0) {
	if (strcmp(c3, "TRF") == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && strcmp(c3, "TRD") == 0) {
	    nbmin = 2;
	}
    } else if (cname && strcmp(c2, "HE") == 0) {
	if (strcmp(c3, "TRD") == 0) {
	    nbmin = 2;
	}
    } else if (sname && strcmp(c2, "OR") == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && strcmp(c2, "UN") == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (strcmp(c2, "GE") == 0) {
	if (strcmp(c3, "QRF") == 0 || strcmp(c3, "RQF") == 0 || 
		strcmp(c3, "LQF") == 0 || strcmp(c3, "QLF") == 
		0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (strcmp(c3, "HRD") == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (strcmp(c3, "BRD") == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (strcmp(c2, "SY") == 0) {
	if (sname && strcmp(c3, "TRD") == 0) {
	    nx = 1;
	}
    } else if (cname && strcmp(c2, "HE") == 0) {
	if (strcmp(c3, "TRD") == 0) {
	    nx = 1;
	}
    } else if (sname && strcmp(c2, "OR") == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nx = 128;
	    }
	}
    } else if (cname && strcmp(c2, "UN") == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (strcmp(c4, "QR") == 0 || strcmp(c4, "RQ") == 0 
		    || strcmp(c4, "LQ") == 0 || strcmp(c4, "QL")
		     == 0 || strcmp(c4, "HR") == 0 || strcmp(c4, "TR") == 0 || 
		     strcmp(c4, "BR") == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) min(*n1,*n2) * (float)1.6);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */

#ifdef __MWERKS__
#pragma warn_unusedarg reset
#endif /* __MWERKS__ */
