/* $Id: genleaftree.c,v 1.1 2005-01-03 00:59:50 paklein Exp $ */
/* genleaftree.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "pspases_f2c.h"

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   genleaftree.f                                                           +/ */
/* /+                                                                           +/ */
/* /+   Written by Mahesh Joshi, U of MN.                                       +/ */
/* /+                                                                           +/ */
/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+ This code is meant to be used solely for educational, research, and       +/ */
/* /+ benchmarking purposes by non-profit institutions and US government        +/ */
/* /+ agencies only.  Use by any other organization requires prior written      +/ */
/* /+ permission from both IBM Corporation and the University of Minnesota.     +/ */
/* /+ The software may not be sold or redistributed.  One may make copies       +/ */
/* /+ of the software or modify it for their use provided that the copies,      +/ */
/* /+ modified or otherwise, are not sold or distributed, are used under the    +/ */
/* /+ same terms and conditions, and this notice and any part of the source     +/ */
/* /+ code that follows this notice are not separated.                          +/ */
/* /+                                                                           +/ */
/* /+ As unestablished research software, this code is provided on an           +/ */
/* /+ ``as is'' basis without warranty of any kind, either expressed or         +/ */
/* /+ implied, including but not limited to implied warranties of               +/ */
/* /+ merchantability and fitness for a particular purpose.  IBM does not       +/ */
/* /+ warrant that the functions contained in this software will meet the       +/ */
/* /+ user's requirements or that the operation of its routines will be         +/ */
/* /+ uninterrupted or error-free.  Acceptance and use of this program          +/ */
/* /+ constitutes the user's understanding that he/she will have no recourse    +/ */
/* /+ to IBM for any actual or consequential damages, including, but not        +/ */
/* /+ limited to, lost profits or savings, arising out of the use or inability  +/ */
/* /+ to use these libraries.  Even if the user informs IBM of the possibility  +/ */
/* /+ of such damages, IBM expects the user to accept the risk of any such      +/ */
/* /+ harm, or the user shall not attempt to use these libraries for any        +/ */
/* /+ purpose.                                                                  +/ */
/* /+                                                                           +/ */
/* /+ The downloading, compiling, or executing any part of this software        +/ */
/* /+ constitutes an implicit agreement to these terms.  These terms and        +/ */
/* /+ conditions are subject to change at any time without prior notice.        +/ */
/* /+                                                                           +/ */
/* /+***************************************************************************+/ */
/* /+ $Id: genleaftree.c,v 1.1 2005-01-03 00:59:50 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<       subroutine GENLEAFTREE(beginnode,size,parent,ptrs,inds,setid) >*/
/* Subroutine */ int genleaftree_(integer *beginnode, integer *size, integer *
	parent, integer *ptrs, integer *inds, integer *setid)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, k, l, node, tpar, lset;

/*<       integer beginnode, size >*/
/*<       integer parent(0:*), ptrs(2,0:*), inds(*), setid(0:*) >*/
/*<       integer i,j,k,l,lset,node,endnode,tpar >*/
/*<       tpar = -1 >*/
    /* Parameter adjustments */
    --inds;
    --ptrs;

    /* Function Body */
    tpar = -1;
/*<       do 30 i = beginnode, beginnode+size-1 >*/
    i__1 = *beginnode + *size - 1;
    for (i__ = *beginnode; i__ <= i__1; ++i__) {
/*<         parent(i) = tpar >*/
	parent[i__] = tpar;
/*<         setid(i) = tpar >*/
	setid[i__] = tpar;
/*<         k = ptrs(1,i) >*/
	k = ptrs[(i__ << 1) + 1];
/*<         do 40 j = k, k + ptrs(2,i) - 1 >*/
	i__2 = k + ptrs[(i__ << 1) + 2] - 1;
	for (j = k; j <= i__2; ++j) {
/*<           l = inds(j) >*/
	    l = inds[j];
/*<           if (l .ge. i) go to 30 >*/
	    if (l >= i__) {
		goto L30;
	    }
/*<           lset = l >*/
	    lset = l;
/*<           node = setid(l) >*/
	    node = setid[l];
/*<           do while (node .ne. tpar) >*/
	    while(node != tpar) {
/*<             setid(lset) = i >*/
		setid[lset] = i__;
/*<             lset = node >*/
		lset = node;
/*<             node = setid(node) >*/
		node = setid[node];
/*<           end do >*/
	    }
/*<           if(lset.ne.i) then >*/
	    if (lset != i__) {
/*<             parent(lset) = i >*/
		parent[lset] = i__;
/*<             setid(lset) = i >*/
		setid[lset] = i__;
/*<           end if >*/
	    }
/*< 40      continue >*/
/* L40: */
	}
/*< 30    continue >*/
L30:
	;
    }
/*<       end >*/
    return 0;
} /* genleaftree_ */

