/* $Id: parsymb.c,v 1.5 2005-01-15 08:18:28 paklein Exp $ */
/* parsymb.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_int.h"
#include <stdio.h>

/* Table of constant values */
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c__11 = 11;
static integer c__21 = 21;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   parsymb.f                                                               +/ */
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
/* /+ $Id: parsymb.c,v 1.5 2005-01-15 08:18:28 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

/*<    >*/
/* Subroutine */ int parsymb_(integer *root, integer *aptrs, integer *ainds, 
	integer *tptrs, integer *tinds, integer *lptrs, integer *linds, 
	integer *lindsptr, integer *lvalsptr, integer *nnz, doublereal *lbal, 
	integer *sizes, integer *sup, integer *supinds, integer *supptr, 
	integer *supindsptr, integer *iu, integer *lsize, integer *dd, 
	integer *n, integer *myid, integer *lgblk, integer *myidr, integer *
	myidc, integer *mystak, integer *resdcol, integer *info, integer *
	cinfo, MPI_Comm *comm)
	
/* dummy arguments added for f2c
	, integer *temp1, integer *temp2, integer *sinds, 
	integer *sptrs, doublereal *temp3, doublereal *temp4, doublereal *
	svals)
*/
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    logical ihadlast;
    integer leafroot, ipartner;
    doublereal opcountl, ssthresh;
    integer i__, j, k, l;
    extern /* Subroutine */ int symbolic6_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, integer *), 
	    mergelists_(integer *, integer *, integer *
	    , integer *, integer *, integer *);
    integer firstindex, resdcolsiz, recvcolsiz;
    extern /* Subroutine */ int all_to_all_union_hc_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    MPI_Comm *);
    integer node, ierr, nnzl;
    integer level, scptr, level2;
    integer bmaskc, bmaskr;
    doublereal osnode;
    integer rowcol, stakst, supbot;
    doublereal opxtra;
    integer suptop, supsiz, kidnode, mymaskc, /* mpistat[4], */ mymaskr, lbotsiz, 
	    stakptr, nuptopl, nnzxtra;
	MPI_Status mpistat;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer itag >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       parameter(itag=1) >*/
/*<       double precision SSTDEF >*/
/*<       parameter(SSTDEF=2.0d2) >*/
/*<       integer lsize,info >*/
/*<       integer node,level,rowcol,dd,N,myid,lgblk,root >*/
/*<       integer supptr,supindsptr,iu,lindsptr,lvalsptr >*/
/*<       integer myidr,myidc,bmaskr,bmaskc,mymaskr,mymaskc >*/
/*<       integer stakptr,kidnode,scptr,comm >*/
/*<       integer ns,is1,is4,nnz,nnzxtra,nnzl >*/
/*<       double precision opcount,opcountl,opxtra,ssthresh,lbal(0:*) >*/
/*<       double precision osnode >*/
/*<       integer i,j,k,l,ipartner,nbrecv,level2,leafroot >*/
/*<       integer lbotsiz , resdcolsiz, recvcolsiz >*/
/*<       integer stakst,supbot,suptop,supsiz >*/
/*<       integer aptrs(2,0:*),ainds(*) >*/
/*<       integer tptrs(3,0:*),tinds(*) >*/
/*<       integer mystak(*) , resdcol(*) >*/
/*<       integer lptrs(3,0:*),linds(*) >*/
/*<       integer sup(*),supinds(*) >*/
/*<       integer cinfo(0:*) >*/
/*<       logical ihadlast >*/
/*<       integer firstindex >*/
/*<       integer nuptopl >*/
/*     integer, allocatable :: temp1(:),temp2(:),sinds(:),sptrs(:) */
/*<       integer temp1(*),temp2(*),sinds(*),sptrs(*) >*/
/*     double precision, allocatable :: temp3(:),temp4(:),svals(:) */
/*<       double precision temp3(*),temp4(*),svals(*) >*/
/*<       integer sizes(0:*),nrecv,maxisizel,maxisizeg,sisize,isize >*/
/*<       integer tsize,siptr,n1,n2,m,ti >*/
/*<       integer mpistat(MPI_STATUS_SIZE),ierr >*/
/*<       node = root >*/
    /* Parameter adjustments */
/*    --svals; */
/*    --temp4; */
/*    --temp3; */
/*    --sptrs; */
/*    --sinds; */
/*    --temp2; */
/*    --temp1; */
    --resdcol;
    --mystak;
    --supinds;
    --sup;
    --linds;
    --lptrs;
    --tinds;
    --tptrs;
    --ainds;
    --aptrs;

    /* Function Body */
    node = *root;
/*<       j=1 >*/
    j = 1;
/*<       do i=1,dd >*/
    i__1 = *dd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         k = j >*/
	k = j;
/*<         do while(tptrs(2,node).eq.1) >*/
	while(tptrs[node * 3 + 2] == 1) {
/*<           mystak(j) = node >*/
	    mystak[j] = node;
/*<           node = tinds(tptrs(1,node)) >*/
	    node = tinds[tptrs[node * 3 + 1]];
/*<           j = j+1 >*/
	    ++j;
/*<         end do >*/
	}
/*<         mystak(j) = node >*/
	mystak[j] = node;
/*<         j = j+1 >*/
	++j;
/*<         mystak(j) = j-k >*/
	mystak[j] = j - k;
/*<         j = j+1 >*/
	++j;
/*<         if(tptrs(2,node).gt.2) then >*/
	if (tptrs[node * 3 + 2] > 2) {
		printf("%d: tree not binarized at level %d for node %d", *myid, i__, node);
		MPI_Abort(*comm, 0);
	}
/*<         if(iand(myid,ishft(1,dd-i)).eq.0) then  >*/
	if ((*myid & lbit_shift((ftnlen)1, *dd - i__)) == 0) {
/*<           node = tinds(tptrs(1,node))  >*/
	    node = tinds[tptrs[node * 3 + 1]];
/*<         else  >*/
	} else {
/*<           node = tinds(tptrs(1,node)+1)  >*/
	    node = tinds[tptrs[node * 3 + 1] + 1];
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       stakptr = j-1  >*/
    stakptr = j - 1;
/*<       do i=0,N-1 >*/
    i__1 = *n - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         lptrs(2,i) = 0 >*/
	lptrs[i__ * 3 + 2] = 0;
/*<       end do >*/
    }
/*<       lvalsptr = 1 >*/
    *lvalsptr = 1;
/*<       nnzxtra = 0 >*/
    nnzxtra = 0;
/*<       nnzl = 1 >*/
    nnzl = 1;
/*<       opxtra = 1.d0 >*/
    opxtra = 1.;
/*<       opcountl = 0.d0 >*/
    opcountl = 0.;
/*<       supptr = 1 >*/
    *supptr = 1;
/*<       lindsptr = 1 >*/
    *lindsptr = 1;
/*<       ssthresh = SSTDEF >*/
    ssthresh = 200.;
/*<       info = 0 >*/
    *info = 0;
/*<       leafroot = node >*/
    leafroot = node;
/*<    >*/
    symbolic6_(&aptrs[1], &ainds[1], &lptrs[1], &linds[1], &sup[1], &resdcol[
	    1], &tptrs[1], &tinds[1], &nnzl, &node, supptr, lvalsptr, &
	    opcountl, &mystak[j], lindsptr, &nnzxtra, &opxtra, &ssthresh, 
	    lsize, info);
/*<       if(info.eq.1) return >*/
    if (*info == 1) {
	return 0;
    }
/*<       j = sup(supptr-4)   >*/
    j = sup[*supptr - 4];
/*<       iu = ishft(lptrs(2,j),1) >*/
    *iu = lptrs[j * 3 + 2] << 1;
/*<       supindsptr = 1 >*/
    *supindsptr = 1;
/*<       lbal(myid) = opcountl >*/
    lbal[*myid] = opcountl;
/*<       if (dd.ne.0) then  >*/
    if (*dd != 0) {
/*<       rowcol = 0  >*/
	rowcol = 0;
/*<       bmaskc = 0 >*/
	bmaskc = 0;
/*<       bmaskr = 0 >*/
	bmaskr = 0;
/*<       mymaskc = 0 >*/
	mymaskc = 0;
/*<       mymaskr = 0 >*/
	mymaskr = 0;
/*<       ipartner = ieor(myid,1) >*/
	ipartner = *myid ^ 1;
/*<       scptr = lptrs(3,node)+1 >*/
	scptr = lptrs[node * 3 + 3] + 1;
/*<       lbotsiz = lptrs(2,node)-1 >*/
	lbotsiz = lptrs[node * 3 + 2] - 1;
/*<    >*/


/*      call mpi_sendrecv(linds(scptr),lbotsiz,MPI_INTEGER,ipartner,
     +                 itag,supinds(supindsptr),N,MPI_INTEGER, 
     +                 ipartner,itag,comm, mpistat,ierr) */
/*	mpi_sendrecv__(&linds[scptr], &lbotsiz, &c__11, &ipartner, &c__1, &
		supinds[*supindsptr], n, &c__11, &ipartner, &c__1, comm, 
		mpistat, &ierr); */
	MPI_Sendrecv(&linds[scptr], lbotsiz, MPI_INT, ipartner, 1, 
		      &supinds[*supindsptr], *n, MPI_INT, ipartner, 1, *comm, &mpistat);

/*<       call mpi_get_count(mpistat,MPI_INTEGER,recvcolsiz,ierr) >*/
	myMPI_Get_count(&mpistat, MPI_INT, &recvcolsiz);

/*<    >*/
	mergelists_(&linds[scptr], &lbotsiz, &supinds[*supindsptr], &
		recvcolsiz, &resdcol[1], &resdcolsiz);
/*<       nuptopl = ishft(1,dd) >*/
	nuptopl = lbit_shift((ftnlen)1, *dd);
/*<       do level=dd-1,0,-1 >*/
	for (level = *dd - 1; level >= 0; --level) {
/*<         stakst = stakptr >*/
	    stakst = stakptr;
/*<         supsiz = mystak(stakptr) >*/
	    supsiz = mystak[stakptr];
/*<         if(mod(supsiz,2).eq.0) then  >*/
	    if (supsiz % 2 == 0) {
/*<           do stakptr = stakst-1,stakst-supsiz,-2 >*/
		i__1 = stakst - supsiz;
		for (stakptr = stakst - 1; stakptr >= i__1; stakptr += -2) {
/*<             kidnode = mystak(stakptr) >*/
		    kidnode = mystak[stakptr];
/*<    >*/
		    mergelists_(&resdcol[1], &resdcolsiz, &ainds[aptrs[(
			    kidnode << 1) + 1]], &aptrs[(kidnode << 1) + 2], &
			    linds[*lindsptr], &lbotsiz);
/*<             kidnode = mystak(stakptr-1) >*/
		    kidnode = mystak[stakptr - 1];
/*<    >*/
		    mergelists_(&linds[*lindsptr], &lbotsiz, &ainds[aptrs[(
			    kidnode << 1) + 1]], &aptrs[(kidnode << 1) + 2], &
			    resdcol[1], &resdcolsiz);
/*<           end do >*/
		}
/*<    >*/
		i__1 = *dd - level;
		all_to_all_union_hc_(&resdcol[1], &supinds[*supindsptr], &
			linds[*lindsptr], n, &resdcolsiz, &lbotsiz, myid, &
			i__1, comm);
/*<         else  >*/
	    } else {
/*<           do stakptr = stakst-1,stakst-supsiz+1,-2 >*/
		i__1 = stakst - supsiz + 1;
		for (stakptr = stakst - 1; stakptr >= i__1; stakptr += -2) {
/*<             kidnode = mystak(stakptr) >*/
		    kidnode = mystak[stakptr];
/*<    >*/
		    mergelists_(&resdcol[1], &resdcolsiz, &ainds[aptrs[(
			    kidnode << 1) + 1]], &aptrs[(kidnode << 1) + 2], &
			    linds[*lindsptr], &lbotsiz);
/*<             kidnode = mystak(stakptr-1) >*/
		    kidnode = mystak[stakptr - 1];
/*<    >*/
		    mergelists_(&linds[*lindsptr], &lbotsiz, &ainds[aptrs[(
			    kidnode << 1) + 1]], &aptrs[(kidnode << 1) + 2], &
			    resdcol[1], &resdcolsiz);
/*<           end do >*/
		}
/*<           kidnode = mystak(stakptr) >*/
		kidnode = mystak[stakptr];
/*<           stakptr = stakptr-1 >*/
		--stakptr;
/*<    >*/
		mergelists_(&resdcol[1], &resdcolsiz, &ainds[aptrs[(kidnode <<
			 1) + 1]], &aptrs[(kidnode << 1) + 2], &linds[*
			lindsptr], &lbotsiz);
/*<    >*/
		i__1 = *dd - level;
		all_to_all_union_hc_(&linds[*lindsptr], &supinds[*supindsptr]
			, &resdcol[1], n, &lbotsiz, &lbotsiz, myid, &i__1, 
			comm);
/*<         end if >*/
	    }
/*<         rowcol = 1-rowcol  >*/
	    rowcol = 1 - rowcol;
/*<         level2 = ishft(dd-level-1,-1) >*/
	    level2 = lbit_shift(*dd - level - 1, (ftnlen)-1);
/*<         if(rowcol.eq.1) then  >*/
	    if (rowcol == 1) {
/*<           bmaskc = ior(bmaskc,ishft(1,level2)) >*/
		bmaskc |= lbit_shift((ftnlen)1, level2);
/*<         else >*/
	    } else {
/*<           bmaskr = ior(bmaskr,ishft(1,level2)) >*/
		bmaskr |= lbit_shift((ftnlen)1, level2);
/*<         end if >*/
	    }
/*<         l = supindsptr+lbotsiz >*/
	    l = *supindsptr + lbotsiz;
/*<         k=1 >*/
	    k = 1;
/*<         do i=0,lbotsiz-1 >*/
	    i__1 = lbotsiz - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
/*<           j = supinds(supindsptr+i) >*/
		j = supinds[*supindsptr + i__];
/*<           if(iand(myidr,bmaskr).eq.iand(ishft(j,-lgblk),bmaskr)) then >*/
		if ((*myidr & bmaskr) == (lbit_shift(j, -(*lgblk)) & bmaskr)) 
			{
/*<             supinds(l+k-1) = j >*/
		    supinds[l + k - 1] = j;
/*<             k = k+1 >*/
		    ++k;
/*<           end if >*/
		}
/*<         end do >*/
	    }
/*<         recvcolsiz = k-1 >*/
	    recvcolsiz = k - 1;
/*<         if(iand(myidr,bmaskr).eq.0 .and. iand(myidc,bmaskc).eq.0) then >*/
	    if ((*myidr & bmaskr) == 0 && (*myidc & bmaskc) == 0) {
/*<           osnode = 0.d0 >*/
		osnode = 0.;
/*<           j = lbotsiz >*/
		j = lbotsiz;
/*<           do i=0,supsiz-1 >*/
		i__1 = supsiz - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<             osnode = osnode + dble(j*j) >*/
		    osnode += (doublereal) (j * j);
/*<             j = j-1 >*/
		    --j;
/*<           end do >*/
		}
/*<           j = ishft(myid,level-dd)+nuptopl >*/
		j = lbit_shift(*myid, level - *dd) + nuptopl;
/*<           nuptopl = nuptopl+ishft(1,level) >*/
		nuptopl += lbit_shift((ftnlen)1, level);
/*<           lbal(j) = osnode >*/
		lbal[j] = osnode;
/*<         end if >*/
	    }
/*<         k=1 >*/
	    k = 1;
/*<         j=0 >*/
	    j = 0;
/*<         ihadlast = .false. >*/
	    ihadlast = FALSE_;
/*<         do i = stakst-1,stakst-supsiz,-1 >*/
	    i__1 = stakst - supsiz;
	    for (i__ = stakst - 1; i__ >= i__1; --i__) {
/*<           kidnode = mystak(i) >*/
		kidnode = mystak[i__];
/*<    >*/
		if ((*myidc & bmaskc) == (lbit_shift(kidnode, -(*lgblk)) & 
			bmaskc)) {
/*<    >*/
		    while(supinds[l + k - 1] < kidnode && k <= recvcolsiz) {
/*<               k = k+1 >*/
			++k;
/*<             end do >*/
		    }
/*<             if(j.eq.0) j=k >*/
		    if (j == 0) {
			j = k;
		    }
/*<             lptrs(2,kidnode) = recvcolsiz-k+1 >*/
		    lptrs[kidnode * 3 + 2] = recvcolsiz - k + 1;
/*<             lptrs(3,kidnode) = lindsptr+k-j >*/
		    lptrs[kidnode * 3 + 3] = *lindsptr + k - j;
/*<             cinfo(kidnode) = 1  >*/
		    cinfo[kidnode] = 1;
/*<             if(.not.ihadlast) then >*/
		    if (! ihadlast) {
/*<               lptrs(1,kidnode) = lvalsptr >*/
			lptrs[kidnode * 3 + 1] = *lvalsptr;
/*<               lvalsptr = lvalsptr+lptrs(2,kidnode) >*/
			*lvalsptr += lptrs[kidnode * 3 + 2];
/*<               ihadlast = .true. >*/
			ihadlast = TRUE_;
/*<               firstindex = k >*/
			firstindex = k;
/*<             else >*/
		    } else {
/*<               lvalsptr = lvalsptr+k-firstindex >*/
			*lvalsptr = *lvalsptr + k - firstindex;
/*<               lptrs(1,kidnode) = lvalsptr >*/
			lptrs[kidnode * 3 + 1] = *lvalsptr;
/*<               lvalsptr = lvalsptr+lptrs(2,kidnode) >*/
			*lvalsptr += lptrs[kidnode * 3 + 2];
/*<             end if >*/
		    }
/*<           else >*/
		} else {
/*<             lptrs(2,kidnode) = 0 >*/
		    lptrs[kidnode * 3 + 2] = 0;
/*<             lptrs(3,kidnode) = lindsptr >*/
		    lptrs[kidnode * 3 + 3] = *lindsptr;
/*<             cinfo(kidnode) = 0  >*/
		    cinfo[kidnode] = 0;
/*<             ihadlast = .false. >*/
		    ihadlast = FALSE_;
/*<           end if >*/
		}
/*<           nnzl = nnzl+lptrs(2,kidnode) >*/
		nnzl += lptrs[kidnode * 3 + 2];
/*<         end do >*/
	    }
/*<         if(j.ne.0) then >*/
	    if (j != 0) {
/*<           k = recvcolsiz-j+1 >*/
		k = recvcolsiz - j + 1;
/*<           do i=0,k-1 >*/
		i__1 = k - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<             linds(lindsptr+i) = supinds(l+j-1+i) >*/
		    linds[*lindsptr + i__] = supinds[l + j - 1 + i__];
/*<           end do >*/
		}
/*<           lindsptr = lindsptr+k >*/
		*lindsptr += k;
/*<         end if >*/
	    }
/*<         supbot = mystak(stakst-1) >*/
	    supbot = mystak[stakst - 1];
/*<         suptop = mystak(stakst-supsiz) >*/
	    suptop = mystak[stakst - supsiz];
/*<         tptrs(3,suptop) = supptr >*/
	    tptrs[suptop * 3 + 3] = *supptr;
/*<         tptrs(3,supbot) = supptr >*/
	    tptrs[supbot * 3 + 3] = *supptr;
/*<         sup(supptr) = supbot            >*/
	    sup[*supptr] = supbot;
/*<         sup(supptr+1) = supsiz          >*/
	    sup[*supptr + 1] = supsiz;
/*<         sup(supptr+2) = supindsptr          >*/
	    sup[*supptr + 2] = *supindsptr;
/*<         sup(supptr+3) = lbotsiz          >*/
	    sup[*supptr + 3] = lbotsiz;
/*<         supptr = supptr + 4 >*/
	    *supptr += 4;
/*<         scptr = supindsptr+supsiz >*/
	    scptr = *supindsptr + supsiz;
/*<         supindsptr = l+lbotsiz >*/
	    *supindsptr = l + lbotsiz;
/*<         lbotsiz = lbotsiz-supsiz >*/
	    lbotsiz -= supsiz;
/*<         if(level.ne.0) then >*/
	    if (level != 0) {
/*<           ipartner = ieor(myid,ishft(1,dd-level)) >*/
		ipartner = *myid ^ lbit_shift((ftnlen)1, *dd - level);
/*<    >*/

/*		mpi_sendrecv__(&supinds[scptr], &lbotsiz, &c__11, &ipartner, &
			c__1, &supinds[*supindsptr], n, &c__11, &ipartner, &
			c__1, comm, mpistat, &ierr); */
		MPI_Sendrecv(&supinds[scptr], lbotsiz, MPI_INT, ipartner, 1, 
		             &supinds[*supindsptr], *n, MPI_INT, ipartner, 1, *comm, &mpistat);

/*<           call mpi_get_count(mpistat,MPI_INTEGER,recvcolsiz,ierr) >*/
		myMPI_Get_count(&mpistat, MPI_INT, &recvcolsiz);
		
/*<    >*/
		mergelists_(&supinds[scptr], &lbotsiz, &supinds[*supindsptr], 
			&recvcolsiz, &resdcol[1], &resdcolsiz);
/*<         end if >*/
	    }
/*<       end do >*/
	}
/*<       end if  >*/
    }
/*<       call mpi_reduce(nnzl,nnz,1,MPI_INTEGER,MPI_SUM,0,comm,ierr) >*/
    MPI_Reduce(&nnzl, nnz, 1, MPI_INT, MPI_SUM, 0, *comm);
/*<       end >*/
    return 0;
} /* parsymb_ */
