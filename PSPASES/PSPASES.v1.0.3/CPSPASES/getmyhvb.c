/* $Id: getmyhvb.c,v 1.1 2005-01-02 23:52:31 paklein Exp $ */
/* getmyhvb.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "pspases_f2c.h"

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   getmyhvb.f                                                              +/ */
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
/* /+ $Id: getmyhvb.c,v 1.1 2005-01-02 23:52:31 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int getmyhvb_(integer *mysnodes, integer *nsupnode, integer *
	sup, integer *ssize, integer *supinds, integer *supindsize, integer *
	tptrs, integer *tinds, integer *n, integer *dd, integer *lgblk, 
	integer *hvbndry, integer *hvbsize, integer *myidh, integer *myidv, 
	integer *myid)
{

    /* Local variables */
    integer poshdiff, posvdiff, i__, j, nb, ph, is, pv, bhp, pch, bip, bhs, 
	    bis, nvb, bvp, pcv, bvs, pos, node, nvbg, hvbp, nvbt, sptr, bnode,
	     level, hsize, vsize, level2, pclean, bmaskh, bmaskv, rowcol, 
	    supbot, vsizer, suptop, supptr, supsiz, pcleanh, pactual, istflag,
	     mymaskh, lbotsiz, mymaskv, currpos;

/*<       implicit none >*/
/*<       integer lendp >*/
/*<       parameter(lendp=8) >*/
/*<       integer nsupnode,ssize,supindsize,dd,lgblk >*/
/*<       integer hvbsize,myidv,myidh,N >*/
/*<       integer sup(*), supinds(*) >*/
/*<       integer tptrs(3,0:*),tinds(*) >*/
/*<       integer mysnodes(*),hvbndry(*) >*/
/*<       integer node,i,j,k,rowcol,level,level2,nb,nb_per_proc >*/
/*<       integer myleft,myright,myup,mydown,bmaskh,hsize,bmaskv,vsize >*/
/*<       integer hvbp,supptr,supbot,supsiz,lbotsiz,sptr,is,bnode >*/
/*<       integer bhs,bvs,bhp,bvp,pclean,currpos,pos,nvbg,pactual >*/
/*<       integer ph,pch,pv,pcv,nvb,nvbt,vsizer,bip,bis,pcleanh >*/
/*<       integer mymaskh,mymaskv,suptop,istflag >*/
/*<       integer myid,nbleft,nbleft1,nbleft2,poshdiff,posvdiff >*/
/*<       level = dd >*/
    /* Parameter adjustments */
    --hvbndry;
    --tinds;
    --tptrs;
    --supinds;
    --sup;
    --mysnodes;

    /* Function Body */
    level = *dd;
/*<       rowcol = 0 >*/
    rowcol = 0;
/*<       bmaskh = 0 >*/
    bmaskh = 0;
/*<       bmaskv = 0 >*/
    bmaskv = 0;
/*<       hsize = 1 >*/
    hsize = 1;
/*<       vsize = 1 >*/
    vsize = 1;
/*<       mymaskh = 0 >*/
    mymaskh = 0;
/*<       mymaskv = 0 >*/
    mymaskv = 0;
/*<       hvbp = 1 >*/
    hvbp = 1;
/*<       do is = nsupnode-1,1,-1 >*/
    for (is = *nsupnode - 1; is >= 1; --is) {
/*<         suptop = mysnodes(is) >*/
	suptop = mysnodes[is];
/*<         supptr = tptrs(3,suptop) >*/
	supptr = tptrs[suptop * 3 + 3];
/*<         if (tptrs(2,sup(supptr)).ne.1) then >*/
	if (tptrs[sup[supptr] * 3 + 2] != 1) {
/*<           level = level-1 >*/
	    --level;
/*<           rowcol = 1-rowcol >*/
	    rowcol = 1 - rowcol;
/*<           level2 = ishft(dd-level-1,-1) >*/
	    level2 = lbit_shift(*dd - level - 1, (ftnlen)-1);
/*<           if(rowcol.eq.1) then       >*/
	    if (rowcol == 1) {
/*<             bmaskh = ior(bmaskh,ishft(1,level2)) >*/
		bmaskh |= lbit_shift((ftnlen)1, level2);
/*<             hsize = ishft(hsize,1) >*/
		hsize <<= 1;
/*<             mymaskh = iand(myidh,bmaskh) >*/
		mymaskh = *myidh & bmaskh;
/*<           else                   >*/
	    } else {
/*<             bmaskv = ior(bmaskv,ishft(1,level2)) >*/
		bmaskv |= lbit_shift((ftnlen)1, level2);
/*<             vsize = ishft(vsize,1) >*/
		vsize <<= 1;
/*<             mymaskv = iand(myidv,bmaskv) >*/
		mymaskv = *myidv & bmaskv;
/*<           end if >*/
	    }
/*<         end if >*/
	}
/*<         supbot = sup(supptr) >*/
	supbot = sup[supptr];
/*<         supsiz = sup(supptr+1) >*/
	supsiz = sup[supptr + 1];
/*<         sptr   = sup(supptr+2) >*/
	sptr = sup[supptr + 2];
/*<         lbotsiz = sup(supptr+3) >*/
	lbotsiz = sup[supptr + 3];
/*<         istflag = 0 >*/
	istflag = 0;
/*<         bhs = hvbp+6 >*/
	bhs = hvbp + 6;
/*<         bhp = bhs >*/
	bhp = bhs;
/*<         i = 0 >*/
	i__ = 0;
/*<         bnode = ishft(supinds(sptr),-lgblk) >*/
	bnode = lbit_shift(supinds[sptr], -(*lgblk));
/*<         pch = iand(bnode,bmaskh) >*/
	pch = bnode & bmaskh;
/*<         pclean  = pch >*/
	pclean = pch;
/*<         pactual = pch >*/
	pactual = pch;
/*<         currpos = mod(mymaskh+hsize-pclean,hsize) >*/
	currpos = (mymaskh + hsize - pclean) % hsize;
/*<         pos = 0 >*/
	pos = 0;
/*<         do while (i .lt. supsiz) >*/
	while(i__ < supsiz) {
/*<           j = i >*/
	    j = i__;
/*<           i = i+1 >*/
	    ++i__;
/*<           do while (i.lt.supsiz) >*/
	    while(i__ < supsiz) {
/*<             bnode = ishft(supinds(sptr+i),-lgblk) >*/
		bnode = lbit_shift(supinds[sptr + i__], -(*lgblk));
/*<             ph = iand(bnode,bmaskh) >*/
		ph = bnode & bmaskh;
/*<             if(pch.ne.ph) goto 10 >*/
		if (pch != ph) {
		    goto L10;
		}
/*<             i = i + 1 >*/
		++i__;
/*<           end do >*/
	    }
/*<  10          if(pos.ge.currpos) then   >*/
L10:
	    if (pos >= currpos) {
/*<             hvbndry(bhp) = supinds(sptr+j) >*/
		hvbndry[bhp] = supinds[sptr + j];
/*<             hvbndry(bhp+2) = supinds(sptr+i-1) >*/
		hvbndry[bhp + 2] = supinds[sptr + i__ - 1];
/*<             if(pos.eq.currpos) then >*/
		if (pos == currpos) {
/*<               hvbndry(bhp+1) = i-j >*/
		    hvbndry[bhp + 1] = i__ - j;
/*<             else >*/
		} else {
/*<               hvbndry(bhp+1) = 0 >*/
		    hvbndry[bhp + 1] = 0;
/*<             end if >*/
		}
/*<             bhp = bhp+3 >*/
		bhp += 3;
/*<             currpos = currpos + hsize >*/
		currpos += hsize;
/*<           end if >*/
	    }
/*<           pclean = mod(pclean+1,hsize) >*/
	    pclean = (pclean + 1) % hsize;
/*<           pos = pos + 1 + mod(ph+hsize-pclean,hsize) >*/
	    pos = pos + 1 + (ph + hsize - pclean) % hsize;
/*<           pclean = ph >*/
	    pclean = ph;
/*<           pch = ph >*/
	    pch = ph;
/*<         end do >*/
	}
/*<         hvbndry(bhs-1) = (bhp-bhs)/3  >*/
	hvbndry[bhs - 1] = (bhp - bhs) / 3;
/*<         bis = bhp+1 >*/
	bis = bhp + 1;
/*<         bvs = bhp+4+lbotsiz >*/
	bvs = bhp + 4 + lbotsiz;
/*<         bip = bis >*/
	bip = bis;
/*<         bvp = bvs >*/
	bvp = bvs;
/*<         node  = supinds(sptr) >*/
	node = supinds[sptr];
/*<         i = 0 >*/
	i__ = 0;
/*<         pos = 0 >*/
	pos = 0;
/*<         bnode = ishft(node,-lgblk) >*/
	bnode = lbit_shift(node, -(*lgblk));
/*<         nvb = 0 >*/
	nvb = 0;
/*<         nvbt = 0 >*/
	nvbt = 0;
/*<         nvbg = 0 >*/
	nvbg = 0;
/*<         pcv = iand(bnode,bmaskv) >*/
	pcv = bnode & bmaskv;
/*<         pch = iand(bnode,bmaskh) >*/
	pch = bnode & bmaskh;
/*<         currpos = mod(mymaskv+vsize-pcv,vsize) >*/
	currpos = (mymaskv + vsize - pcv) % vsize;
/*<         pclean = pcv >*/
	pclean = pcv;
/*<         pcleanh = pch >*/
	pcleanh = pch;
/*<         do while(i.lt.supsiz) >*/
	while(i__ < supsiz) {
/*<           j = i >*/
	    j = i__;
/*<           if(pos.eq.currpos) then >*/
	    if (pos == currpos) {
/*<             hvbndry(bip) = node >*/
		hvbndry[bip] = node;
/*<             bip = bip+1 >*/
		++bip;
/*<           end if >*/
	    }
/*<           i = i+1 >*/
	    ++i__;
/*<           do while (i.lt.supsiz) >*/
	    while(i__ < supsiz) {
/*<             node  = supinds(sptr+i) >*/
		node = supinds[sptr + i__];
/*<             bnode = ishft(node,-lgblk) >*/
		bnode = lbit_shift(node, -(*lgblk));
/*<             ph = iand(bnode,bmaskh) >*/
		ph = bnode & bmaskh;
/*<             pv = iand(bnode,bmaskv) >*/
		pv = bnode & bmaskv;
/*<             if((ph.ne.pch).or.(pv.ne.pcv)) goto 20 >*/
		if (ph != pch || pv != pcv) {
		    goto L20;
		}
/*<             if(pos.eq.currpos) then >*/
		if (pos == currpos) {
/*<               hvbndry(bip) = node >*/
		    hvbndry[bip] = node;
/*<               bip = bip+1 >*/
		    ++bip;
/*<             end if >*/
		}
/*<             i = i + 1 >*/
		++i__;
/*<           end do >*/
	    }
/*<  20          if(pos.ge.currpos) then   >*/
L20:
	    if (pos >= currpos) {
/*<             hvbndry(bvp) = supinds(sptr+j) >*/
		hvbndry[bvp] = supinds[sptr + j];
/*<             hvbndry(bvp+2) = supinds(sptr+i-1) >*/
		hvbndry[bvp + 2] = supinds[sptr + i__ - 1];
/*<             if(pos.eq.currpos) then >*/
		if (pos == currpos) {
/*<               hvbndry(bvp+1) = i-j >*/
		    hvbndry[bvp + 1] = i__ - j;
/*<             else >*/
		} else {
/*<               hvbndry(bvp+1) = 0 >*/
		    hvbndry[bvp + 1] = 0;
/*<             end if >*/
		}
/*<             bvp = bvp+3 >*/
		bvp += 3;
/*<             nvbt = nvbt+1 >*/
		++nvbt;
/*<             currpos = currpos + vsize >*/
		currpos += vsize;
/*<           end if >*/
	    }
/*<           pcleanh = mod(pcleanh+1,hsize) >*/
	    pcleanh = (pcleanh + 1) % hsize;
/*<           if(pv.ne.pcv .and. i.lt.supsiz) then >*/
	    if (pv != pcv && i__ < supsiz) {
/*<             pclean = mod(pclean+1,vsize) >*/
		pclean = (pclean + 1) % vsize;
/*<             posvdiff = mod(pv+vsize-pclean,vsize) >*/
		posvdiff = (pv + vsize - pclean) % vsize;
/*<             poshdiff = mod(ph+hsize-pcleanh,hsize) >*/
		poshdiff = (ph + hsize - pcleanh) % hsize;
/*<             pos = pos + 1 + posvdiff >*/
		pos = pos + 1 + posvdiff;
/*<             if(poshdiff.ne.posvdiff) then  >*/
		if (poshdiff != posvdiff) {
/*<               hvbndry(bvp)   = node >*/
		    hvbndry[bvp] = node;
/*<               hvbndry(bvp+1) = 0 >*/
		    hvbndry[bvp + 1] = 0;
/*<               hvbndry(bvp+2) = node >*/
		    hvbndry[bvp + 2] = node;
/*<               bvp = bvp+3 >*/
		    bvp += 3;
/*<               nvbt = nvbt+1 >*/
		    ++nvbt;
/*<               nvbg = nvbg+vsize  >*/
		    nvbg += vsize;
/*<             end if >*/
		}
/*<             pclean = pv >*/
		pclean = pv;
/*<           else >*/
	    } else {
/*<             pos = pos + vsize >*/
		pos += vsize;
/*<           end if >*/
	    }
/*<           pcleanh = ph >*/
	    pcleanh = ph;
/*<           nvbg = nvbg+1 >*/
	    ++nvbg;
/*<           pch = ph >*/
	    pch = ph;
/*<           pcv = pv >*/
	    pcv = pv;
/*<         end do >*/
	}
/*<         nvb = nvbt >*/
	nvb = nvbt;
/*<         vsizer = 0 >*/
	vsizer = 0;
/*<         i = supsiz >*/
	i__ = supsiz;
/*<         node  = supinds(sptr+i) >*/
	node = supinds[sptr + i__];
/*<         bnode = ishft(node,-lgblk) >*/
	bnode = lbit_shift(node, -(*lgblk));
/*<         pclean = iand(bnode,bmaskv) >*/
	pclean = bnode & bmaskv;
/*<         pactual = pclean >*/
	pactual = pclean;
/*<         currpos = mod(mymaskv+vsize-pclean,vsize) >*/
	currpos = (mymaskv + vsize - pclean) % vsize;
/*<         pos = 0 >*/
	pos = 0;
/*<         do while(i.lt.lbotsiz) >*/
	while(i__ < lbotsiz) {
/*<           j = i >*/
	    j = i__;
/*<           if(pos.eq.currpos) then >*/
	    if (pos == currpos) {
/*<             hvbndry(bip) = node >*/
		hvbndry[bip] = node;
/*<             bip = bip+1 >*/
		++bip;
/*<           end if >*/
	    }
/*<           i = i+1  >*/
	    ++i__;
/*<           do while(i.lt.lbotsiz) >*/
	    while(i__ < lbotsiz) {
/*<             node  = supinds(sptr+i) >*/
		node = supinds[sptr + i__];
/*<             bnode = ishft(node,-lgblk) >*/
		bnode = lbit_shift(node, -(*lgblk));
/*<             pv = iand(bnode,bmaskv) >*/
		pv = bnode & bmaskv;
/*<             if(pv.ne.pactual) goto 30 >*/
		if (pv != pactual) {
		    goto L30;
		}
/*<             if(pos.eq.currpos) then >*/
		if (pos == currpos) {
/*<               hvbndry(bip) = node >*/
		    hvbndry[bip] = node;
/*<               bip = bip+1 >*/
		    ++bip;
/*<             end if >*/
		}
/*<             i = i+1 >*/
		++i__;
/*<           end do >*/
	    }
/*<  30          if(pos.ge.currpos) then >*/
L30:
	    if (pos >= currpos) {
/*<             hvbndry(bvp) = supinds(sptr+j) >*/
		hvbndry[bvp] = supinds[sptr + j];
/*<             hvbndry(bvp+2) = supinds(sptr+i-1) >*/
		hvbndry[bvp + 2] = supinds[sptr + i__ - 1];
/*<             if(pos.eq.currpos) then >*/
		if (pos == currpos) {
/*<               hvbndry(bvp+1) = i-j >*/
		    hvbndry[bvp + 1] = i__ - j;
/*<             else >*/
		} else {
/*<               hvbndry(bvp+1) = 0 >*/
		    hvbndry[bvp + 1] = 0;
/*<             end if >*/
		}
/*<             vsizer = vsizer+hvbndry(bvp+1) >*/
		vsizer += hvbndry[bvp + 1];
/*<             nvb = nvb+1 >*/
		++nvb;
/*<             bvp = bvp+3 >*/
		bvp += 3;
/*<             currpos = currpos + vsize >*/
		currpos += vsize;
/*<           end if >*/
	    }
/*<           pclean = mod(pclean+1,vsize) >*/
	    pclean = (pclean + 1) % vsize;
/*<           pactual = pv >*/
	    pactual = pv;
/*<           pos = pos + 1 + mod(pv+vsize-pclean,vsize) >*/
	    pos = pos + 1 + (pv + vsize - pclean) % vsize;
/*<           pclean = pv >*/
	    pclean = pv;
/*<         end do >*/
	}
/*<         nb=1 >*/
	nb = 1;
/*<         bnode = ishft(suptop,-lgblk) >*/
	bnode = lbit_shift(suptop, -(*lgblk));
/*<         pv = iand(bnode,bmaskv) >*/
	pv = bnode & bmaskv;
/*<         ph = iand(bnode,bmaskh) >*/
	ph = bnode & bmaskh;
/*<         if(pv.eq.mymaskv) then >*/
	if (pv == mymaskv) {
/*<           istflag = 2 >*/
	    istflag = 2;
/*<           if(ph.eq.mymaskh) istflag = 1 >*/
	    if (ph == mymaskh) {
		istflag = 1;
	    }
/*<         end if >*/
	}
/*<         hvbndry(bvs-1)  = nvb       >*/
	hvbndry[bvs - 1] = nvb;
/*<         hvbndry(bvs-2)  = nvb-nvbt  >*/
	hvbndry[bvs - 2] = nvb - nvbt;
/*<         hvbndry(bvs-3)  = vsizer    >*/
	hvbndry[bvs - 3] = vsizer;
/*<         hvbndry(hvbp)   = bvp+1 >*/
	hvbndry[hvbp] = bvp + 1;
/*<         hvbndry(hvbp+1) = bvs >*/
	hvbndry[hvbp + 1] = bvs;
/*<         hvbndry(hvbp+2) = nb >*/
	hvbndry[hvbp + 2] = nb;
/*<         hvbndry(hvbp+3) = istflag >*/
	hvbndry[hvbp + 3] = istflag;
/*<         hvbndry(hvbp+4) = lbotsiz >*/
	hvbndry[hvbp + 4] = lbotsiz;
/*<         hvbndry(bis-1)  = bip-bis   >*/
	hvbndry[bis - 1] = bip - bis;
/*<         hvbndry(bvp) = hvbp >*/
	hvbndry[bvp] = hvbp;
/*<         hvbp = bvp+1 >*/
	hvbp = bvp + 1;
/*<       end do >*/
    }
/*<       hvbsize = hvbp-1 >*/
    *hvbsize = hvbp - 1;
/*<       end  >*/
    return 0;
} /* getmyhvb_ */

