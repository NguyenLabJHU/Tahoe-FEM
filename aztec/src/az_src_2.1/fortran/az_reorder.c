/* az_reorder.f -- translated by f2c (version 19971204).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"


/*     This stuff comes from George & Liu. It basically corresponds */
/*     to computing a Reverse Cuthill-McKee Ordering corresponding */
/*     to the graph of a matrix. */

/* Subroutine */ int az_rcm__(integer *root, integer *xadj, integer *adjncy, 
	integer *mask, integer *perm, integer *ccsize, integer *deg)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer node, fnbr, lnbr, i__, j, k, l;
    extern /* Subroutine */ int az_degree__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer lperm, jstop, jstrt, lbegin, lvlend, nbr;

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */

    /* Parameter adjustments */
    --deg;
    --perm;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    az_degree__(root, &xadj[1], &adjncy[1], &mask[1], &deg[1], ccsize, &perm[
	    1]);

    mask[*root] = 0;
    if (*ccsize <= 1) {
	return 0;
    }
    lvlend = 0;
    lnbr = 1;

L100:
    lbegin = lvlend + 1;

    lvlend = lnbr;

    i__1 = lvlend;
    for (i__ = lbegin; i__ <= i__1; ++i__) {
	node = perm[i__];
	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	fnbr = lnbr + 1;

	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    nbr = adjncy[j];
	    if (mask[nbr] == 0) {
		goto L200;
	    }
	    ++lnbr;
	    mask[nbr] = 0;
	    perm[lnbr] = nbr;
L200:
	    ;
	}

	if (fnbr >= lnbr) {
	    goto L600;
	}
	k = fnbr;

L300:
	l = k;

	++k;
	nbr = perm[k];

L400:
	if (l < fnbr) {
	    goto L500;
	}

	lperm = perm[l];
	if (deg[lperm] <= deg[nbr]) {
	    goto L500;
	}
	perm[l + 1] = lperm;
	--l;
	goto L400;

L500:
	perm[l + 1] = nbr;

	if (k < lnbr) {
	    goto L300;
	}
L600:
	;
    }

    if (lnbr > lvlend) {
	goto L100;
    }
    k = *ccsize / 2;
    l = *ccsize;

    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lperm = perm[l];
	perm[l] = perm[i__];
	perm[i__] = lperm;
	--l;
/* L700: */
    }


    return 0;
} /* az_rcm__ */

/* Subroutine */ int az_degree__(integer *root, integer *xadj, integer *
	adjncy, integer *mask, integer *deg, integer *ccsize, integer *ls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ideg, node, i__, j, jstop, jstrt, lbegin, lvlend, lvsize, 
	    nbr;

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

    /* Parameter adjustments */
    --ls;
    --deg;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    ls[1] = *root;
    xadj[*root] = -xadj[*root];
    lvlend = 0;
    *ccsize = 1;

L100:
    lbegin = lvlend + 1;

    lvlend = *ccsize;
    i__1 = lvlend;
    for (i__ = lbegin; i__ <= i__1; ++i__) {
	node = ls[i__];
	jstrt = -xadj[node];
	jstop = (i__2 = xadj[node + 1], abs(i__2)) - 1;
	ideg = 0;
	if (jstop < jstrt) {
	    goto L300;
	}

	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    nbr = adjncy[j];
	    if (mask[nbr] == 0) {
		goto L200;
	    }
	    ++ideg;
	    if (xadj[nbr] < 0) {
		goto L200;
	    }
	    xadj[nbr] = -xadj[nbr];
	    ++(*ccsize);
	    ls[*ccsize] = nbr;
L200:
	    ;
	}

L300:
	deg[node] = ideg;

/* L400: */
    }

    lvsize = *ccsize - lvlend;
    if (lvsize > 0) {
	goto L100;
    }

    i__1 = *ccsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	node = ls[i__];
	xadj[node] = -xadj[node];
/* L500: */
    }


    return 0;
} /* az_degree__ */

/* Subroutine */ int az_fnroot__(integer *root, integer *xadj, integer *
	adjncy, integer *mask, integer *nlvl, integer *xls, integer *ls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ndeg, node, j, k, nabor, kstop, jstrt, kstrt;
    extern /* Subroutine */ int az_rootls__(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer mindeg, ccsize, nunlvl;

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ls;
    --xls;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    az_rootls__(root, &xadj[1], &adjncy[1], &mask[1], nlvl, &xls[1], &ls[1]);

    ccsize = xls[*nlvl + 1] - 1;
    if (*nlvl == 1 || *nlvl == ccsize) {
	return 0;
    }

L100:
    jstrt = xls[*nlvl];

    mindeg = ccsize;
    *root = ls[jstrt];
    if (ccsize == jstrt) {
	goto L400;
    }

    i__1 = ccsize;
    for (j = jstrt; j <= i__1; ++j) {
	node = ls[j];
	ndeg = 0;
	kstrt = xadj[node];
	kstop = xadj[node + 1] - 1;

	i__2 = kstop;
	for (k = kstrt; k <= i__2; ++k) {
	    nabor = adjncy[k];
	    if (mask[nabor] > 0) {
		++ndeg;
	    }
/* L200: */
	}

	if (ndeg >= mindeg) {
	    goto L300;
	}
	*root = node;
	mindeg = ndeg;
L300:
	;
    }

L400:
    az_rootls__(root, &xadj[1], &adjncy[1], &mask[1], &nunlvl, &xls[1], &ls[1]
	    );

    if (nunlvl <= *nlvl) {
	return 0;
    }
    *nlvl = nunlvl;
    if (*nlvl < ccsize) {
	goto L100;
    }


    return 0;
} /* az_fnroot__ */

/* Subroutine */ int az_rootls__(integer *root, integer *xadj, integer *
	adjncy, integer *mask, integer *nlvl, integer *xls, integer *ls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer node, i__, j, jstop, jstrt, lbegin, ccsize, lvlend, lvsize,
	     nbr;

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ls;
    --xls;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    mask[*root] = 0;
    ls[1] = *root;
    *nlvl = 0;
    lvlend = 0;
    ccsize = 1;

L200:
    lbegin = lvlend + 1;

    lvlend = ccsize;
    ++(*nlvl);
    xls[*nlvl] = lbegin;

    i__1 = lvlend;
    for (i__ = lbegin; i__ <= i__1; ++i__) {
	node = ls[i__];
	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	if (jstop < jstrt) {
	    goto L400;
	}

	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    nbr = adjncy[j];
	    if (mask[nbr] == 0) {
		goto L300;
	    }
	    ++ccsize;
	    ls[ccsize] = nbr;
	    mask[nbr] = 0;
L300:
	    ;
	}

L400:
	;
    }
    lvsize = ccsize - lvlend;
    if (lvsize > 0) {
	goto L200;
    }
    xls[*nlvl + 1] = lvlend + 1;

    i__1 = ccsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	node = ls[i__];
	mask[node] = 1;
/* L500: */
    }


    return 0;
} /* az_rootls__ */
