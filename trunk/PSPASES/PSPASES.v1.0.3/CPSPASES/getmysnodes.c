/* $Id: getmysnodes.c,v 1.1 2005-01-03 00:33:58 paklein Exp $ */
/* getmysnodes.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "pspases_f2c.h"

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   getmysnodes.f                                                           +/ */
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
/* /+ $Id: getmysnodes.c,v 1.1 2005-01-03 00:33:58 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int getmysnodes_(integer *root, integer *sup, integer *tinds,
	 integer *tptrs, integer *n, integer *supsize, integer *mysnodes, 
	integer *nsupnode, integer *lrud, integer *dd, integer *maxhsize, 
	integer *maxvsize, integer *ns, integer *myid)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, jd, jl, jr, ju, node, nbrp, myup, level, level2, myleft, 
	    rowcol, mydown, supptr, myright;

/*<       implicit none >*/
/*<       integer root,N,supsize,nsupnode,dd,maxhsize,maxvsize,ns,myid >*/
/*<       integer sup(*),tinds(*),tptrs(3,0:*),mysnodes(*) >*/
/*<       integer lrud(*) >*/
/*<       integer mydown,myup,myright,myleft,level,level2,supptr,rowcol >*/
/*<       integer i,j,node,nbrp,jl,jr,ju,jd >*/
/*<       maxhsize = 0 >*/
    /* Parameter adjustments */
    --lrud;
    --mysnodes;
    --tptrs;
    --tinds;
    --sup;

    /* Function Body */
    *maxhsize = 0;
/*<       maxvsize = 0 >*/
    *maxvsize = 0;
/*<       j = 1 >*/
    j = 1;
/*<       node = root >*/
    node = *root;
/*<       mysnodes(j) = node >*/
    mysnodes[j] = node;
/*<       j = j+1 >*/
    ++j;
/*<       supptr = tptrs(3,node) >*/
    supptr = tptrs[node * 3 + 3];
/*<       maxhsize = max(sup(supptr+1),maxhsize) >*/
/* Computing MAX */
    i__1 = sup[supptr + 1];
    *maxhsize = max(i__1,*maxhsize);
/*<       maxvsize = max(sup(supptr+3),maxvsize) >*/
/* Computing MAX */
    i__1 = sup[supptr + 3];
    *maxvsize = max(i__1,*maxvsize);
/*<       node = sup(supptr)  >*/
    node = sup[supptr];
/*<       do i=1,dd >*/
    i__1 = *dd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         do while (tptrs(2,node).eq.1)  >*/
	while(tptrs[node * 3 + 2] == 1) {
/*<           node = tinds(tptrs(1,node))  >*/
	    node = tinds[tptrs[node * 3 + 1]];
/*<           mysnodes(j) = node >*/
	    mysnodes[j] = node;
/*<           j = j+1 >*/
	    ++j;
/*<           supptr = tptrs(3,node) >*/
	    supptr = tptrs[node * 3 + 3];
/*<           maxhsize = max(sup(supptr+1),maxhsize) >*/
/* Computing MAX */
	    i__2 = sup[supptr + 1];
	    *maxhsize = max(i__2,*maxhsize);
/*<           maxvsize = max(sup(supptr+3),maxvsize) >*/
/* Computing MAX */
	    i__2 = sup[supptr + 3];
	    *maxvsize = max(i__2,*maxvsize);
/*<           node = sup(supptr)  >*/
	    node = sup[supptr];
/*<         end do >*/
	}
/*<         if(iand(myid,ishft(1,dd-i)).eq.0) then        >*/
	if ((*myid & lbit_shift((ftnlen)1, *dd - i__)) == 0) {
/*<           node = tinds(tptrs(1,node)) >*/
	    node = tinds[tptrs[node * 3 + 1]];
/*<         else                                     >*/
	} else {
/*<           node = tinds(tptrs(1,node)+1) >*/
	    node = tinds[tptrs[node * 3 + 1] + 1];
/*<         end if >*/
	}
/*<         mysnodes(j) = node >*/
	mysnodes[j] = node;
/*<         j = j+1 >*/
	++j;
/*<         supptr = tptrs(3,node) >*/
	supptr = tptrs[node * 3 + 3];
/*<         maxhsize = max(sup(supptr+1),maxhsize) >*/
/* Computing MAX */
	i__2 = sup[supptr + 1];
	*maxhsize = max(i__2,*maxhsize);
/*<         maxvsize = max(sup(supptr+3),maxvsize) >*/
/* Computing MAX */
	i__2 = sup[supptr + 3];
	*maxvsize = max(i__2,*maxvsize);
/*<         node = sup(supptr) >*/
	node = sup[supptr];
/*<       end do >*/
    }
/*<       nsupnode = j-1  >*/
    *nsupnode = j - 1;
/*<       rowcol = 0 >*/
    rowcol = 0;
/*<       myleft = myid >*/
    myleft = *myid;
/*<       myright = myid >*/
    myright = *myid;
/*<       myup = myid >*/
    myup = *myid;
/*<       mydown = myid >*/
    mydown = *myid;
/*<       jl = 0 >*/
    jl = 0;
/*<       jr = 0 >*/
    jr = 0;
/*<       ju = 1 >*/
    ju = 1;
/*<       jd = 1 >*/
    jd = 1;
/*<       nbrp = 1 >*/
    nbrp = 1;
/*<       do level=dd-1,0,-1 >*/
    for (level = *dd - 1; level >= 0; --level) {
/*<         rowcol = 1-rowcol >*/
	rowcol = 1 - rowcol;
/*<         level2 = ishft(dd-level-1,-1) >*/
	level2 = lbit_shift(*dd - level - 1, (ftnlen)-1);
/*<         if(rowcol.eq.1) then       >*/
	if (rowcol == 1) {
/*<           if(myleft.ge.myid) then >*/
	    if (myleft >= *myid) {
/*<             myleft = ieor(myleft,ishft(1,jl)) >*/
		myleft ^= lbit_shift((ftnlen)1, jl);
/*<             jl = jl+2 >*/
		jl += 2;
/*<           endif >*/
	    }
/*<           if(myright.le.myid) then >*/
	    if (myright <= *myid) {
/*<             myright = ieor(myright,ishft(1,jr)) >*/
		myright ^= lbit_shift((ftnlen)1, jr);
/*<             jr = jr+2 >*/
		jr += 2;
/*<           end if >*/
	    }
/*<         else                   >*/
	} else {
/*<           if(myup.ge.myid) then >*/
	    if (myup >= *myid) {
/*<             myup = ieor(myup,ishft(1,ju)) >*/
		myup ^= lbit_shift((ftnlen)1, ju);
/*<             ju = ju+2 >*/
		ju += 2;
/*<           end if >*/
	    }
/*<           if(mydown.le.myid) then >*/
	    if (mydown <= *myid) {
/*<             mydown = ieor(mydown,ishft(1,jd)) >*/
		mydown ^= lbit_shift((ftnlen)1, jd);
/*<             jd = jd+2 >*/
		jd += 2;
/*<           end if >*/
	    }
/*<         endif >*/
	}
/*<         lrud(nbrp)   = myleft >*/
	lrud[nbrp] = myleft;
/*<         lrud(nbrp+1) = myright >*/
	lrud[nbrp + 1] = myright;
/*<         lrud(nbrp+2) = myup >*/
	lrud[nbrp + 2] = myup;
/*<         lrud(nbrp+3) = mydown >*/
	lrud[nbrp + 3] = mydown;
/*<         nbrp = nbrp + 4 >*/
	nbrp += 4;
/*<       end do >*/
    }
/*<       end >*/
    return 0;
} /* getmysnodes_ */

