/* $Id: parelimh1i.c,v 1.1 2005-01-04 17:46:35 paklein Exp $ */
/* parelimh1i.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* debugging */
#undef __DO_DEBUG__
/* #define __DO_DEBUG__ 1 */

#include "mpi.h"
#include "pspases_f2c.h"
#include <stdio.h>

/* Table of constant values */
static integer c__5 = 5;
static doublereal c_b15 = 1.;
static integer c__4 = 4;
static doublereal c_b22 = -1.;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__27 = 27;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   parelimh1i.f                                                            +/ */
/* /+                                                                           +/ */
/* /+   Written by Anshul Gupta, IBM Corp.                                      +/ */
/* /+   Modified by Mahesh Joshi, U of MN.                                      +/ */
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
/* /+ $Id: parelimh1i.c,v 1.1 2005-01-04 17:46:35 paklein Exp $ +/ */
/* /+***************************************************************************+/ */

static integer lbit_shift(integer a, integer b) {
	return b >= 0 ? a << b : (integer)((uinteger)a >> -b);
};

static integer max(integer a, integer b) {
	return (a > b) ? a : b;
}

static integer min(integer a, integer b) {
	return (a < b) ? a : b;
}

/*<    >*/
/* Subroutine */ int parelimh1_(doublereal *wmem, doublereal *vbuf_s__, 
	doublereal *hbuf_s__, doublereal *vbuf_r__, doublereal *hbuf_r__, 
	integer *supinds, integer *tptrs, integer *tinds, integer *cinfo, 
	integer *stak, integer *nstak, doublereal *avals, integer *aptrs, 
	integer *ainds, doublereal *lvals, integer *lptrs, integer *linds, 
	integer *sup, integer *stakptr, integer *nptr, integer *mask1, 
	integer *fptr, integer *ldf, integer *nrows, integer *ncols, integer *
	rsuptr, integer *csuptr, integer *myid, integer *myleft, integer *
	myright, integer *myup, integer *mydown, integer *diagproc, integer *
	bufsize, integer *wmemsize, integer *locinds, integer *prof, 
	doublereal *dfopts, integer *ifopts, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer msgtype1, msgtype2;
    integer npending;
    integer csuptr_u__, rsuptr_u__, j, i__, k, l;
    integer ii, jj, ldb, lda /*, req[4] */;
    MPI_Request req[4];
    
    integer locf, node, rank, info;
    extern /* Subroutine */ int mydc_(integer *, doublereal *, doublereal *);
    integer ireq, ierr, uptr;
    extern /* Subroutine */ int dpack_(doublereal *, doublereal *, integer *, 
	    integer *, integer *), dgemm_(char *, char *, integer *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    integer inode, jnode, nclim;
    extern /* Subroutine */ int myddc_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    integer msgid, index, nrlim, nprof;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    integer nsent1, nsent2, buflen;
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    integer locptr, bottom, nbytes, newsup, newlocf, ncols_u__ /*, statall[16] */
	    /* was [4][4] */, newrank, /* mpistat[4], */ msgtype;
	MPI_Status mpistat[4], statall[16];

/*<       integer supinds(*),tptrs(3,0:*),tinds(*),stak(3,*),nstak(*) >*/
/*<       integer cinfo(0:*) >*/
/*<       double precision lvals(*),avals(*),dfopts(7) >*/
/*<       integer aptrs(2,0:*),ainds(*),lptrs(3,0:*),linds(*),sup(*) >*/
/*<       integer stakptr,nptr,mask1,fptr,ldf,nrows,ncols,comm >*/
/*<       integer rsuptr,csuptr,myid,myleft,myright,myup,mydown >*/
/*<       integer diagproc,bufsize,wmemsize,locinds(*) >*/
/*<       integer prof(4,*),ifopts(5) >*/
/*<       double precision wmem(*),vbuf_s(*),hbuf_s(*)  >*/
/*<       double precision vbuf_r(*),hbuf_r(*) >*/
/*<       integer rank, uptr, bottom, rsuptr_u, csuptr_u, buflen >*/
/*<       include 'mpif.h' >*/
/*<       integer mpistat(MPI_STATUS_SIZE),req(4),ierr >*/
/* -*- fortran -*- */

/* double precision functions */
/*<       integer statall(MPI_STATUS_SIZE,4) >*/
/*<       do ireq=1,4 >*/
    /* Parameter adjustments */
    --ifopts;
    --dfopts;
    prof -= 5;
    --locinds;
    --sup;
    --linds;
    --lptrs;
    --lvals;
    --ainds;
    --aptrs;
    --avals;
    --nstak;
    stak -= 4;
    --tinds;
    --tptrs;
    --supinds;
    --hbuf_r__;
    --vbuf_r__;
    --hbuf_s__;
    --vbuf_s__;
    --wmem;

    /* Function Body */
    for (ireq = 1; ireq <= 4; ++ireq) {
/*<         req(ireq) = MPI_REQUEST_NULL >*/
	req[ireq - 1] = MPI_REQUEST_NULL;
/*<       end do >*/
    }
/*<       buflen = ishft(bufsize,3) >*/
    buflen = *bufsize << 3;
/*<       node = nstak(nptr) >*/
    node = nstak[*nptr];
/*<       nclim = csuptr + ncols >*/
    nclim = *csuptr + *ncols;
/*<       nrlim = rsuptr + nrows >*/
    nrlim = *rsuptr + *nrows;
/*<       newsup = 1 >*/
    newsup = 1;
/*< 1000  rank = 0 >*/
L1000:
    rank = 0;
/*<       bottom = node >*/
    bottom = node;
/*<       j = tptrs(3,node) >*/
    j = tptrs[node * 3 + 3];
/*<       if (j .ne. 0) then >*/
    if (j != 0) {
/*<         if (sup(j) .eq. node) then >*/
	if (sup[j] == node) {
/*<           if (sup(j+1) .eq. 1) then >*/
	    if (sup[j + 1] == 1) {
/*<             nptr = nptr - 1 >*/
		--(*nptr);
/*<             rank = 1 >*/
		rank = 1;
/*<             if (nptr .gt. 0) node = nstak(nptr) >*/
		if (*nptr > 0) {
		    node = nstak[*nptr];
		}
/*<           end if >*/
	    }
/*<         else >*/
	} else {
/*<           nptr = nptr - 1 >*/
	    --(*nptr);
/*<           rank = 1 >*/
	    rank = 1;
/*<           if (nptr .gt. 0) node = nstak(nptr) >*/
	    if (*nptr > 0) {
		node = nstak[*nptr];
	    }
/*<         end if >*/
	}
/*<       end if >*/
    }
/*<       if (rank .eq. 0) then >*/
    if (rank == 0) {
/*<         index = iand(bottom,mask1) >*/
	index = bottom & *mask1;
/*<         nptr = nptr - 1 >*/
	--(*nptr);
/*<         rank = 1 >*/
	rank = 1;
/*<         if (nptr .gt. 0) then >*/
	if (*nptr > 0) {
/*<           node = nstak(nptr) >*/
	    node = nstak[*nptr];
/*<           do while (iand(mask1,node) .eq. index) >*/
	    while((*mask1 & node) == index) {
/*<             if (tptrs(3,node) .ne. 0) then >*/
		if (tptrs[node * 3 + 3] != 0) {
/*<               nptr = nptr - 1 >*/
		    --(*nptr);
/*<               rank = rank + 1 >*/
		    ++rank;
/*<               if (nptr .gt. 0) node = nstak(nptr) >*/
		    if (*nptr > 0) {
			node = nstak[*nptr];
		    }
/*<               goto 10 >*/
		    goto L10;
/*<             end if >*/
		}
/*<             nptr = nptr - 1 >*/
		--(*nptr);
/*<             rank = rank + 1 >*/
		++rank;
/*<             node = nstak(nptr) >*/
		node = nstak[*nptr];
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<       end if >*/
    }
/*< 10    if (newsup .eq. 1) then >*/
L10:
    if (newsup == 1) {
/*<         newsup = 0 >*/
	newsup = 0;
/*<         if (cinfo(bottom) .ne. 1) then >*/
	if (cinfo[bottom] != 1) {
/*<           nprof = 0 >*/
	    nprof = 0;
/*<           if (diagproc .ne. 1) then >*/
	    if (*diagproc != 1) {
/*<             jj = 1 >*/
		jj = 1;
/*<             newlocf = fptr >*/
		newlocf = *fptr;
/*<             i = csuptr >*/
		i__ = *csuptr;
/*<             j = rsuptr >*/
		j = *rsuptr;
/*<             do while (i .lt. nclim) >*/
		while(i__ < nclim) {
/*<               nprof = nprof + 1 >*/
		    ++nprof;
/*<               locf = newlocf >*/
		    locf = newlocf;
/*<               prof(1,nprof) = newlocf >*/
		    prof[(nprof << 2) + 1] = newlocf;
/*<               newrank = 0 >*/
		    newrank = 0;
/*<               inode = supinds(j) >*/
		    inode = supinds[j];
/*<               do while (i+newrank .lt. nclim) >*/
		    while(i__ + newrank < nclim) {
/*<                 if (supinds(i+newrank) .gt. inode) goto 135 >*/
			if (supinds[i__ + newrank] > inode) {
			    goto L135;
			}
/*<                 newrank = newrank + 1 >*/
			++newrank;
/*<               end do >*/
		    }
/*< 135           i = i + newrank >*/
L135:
		    i__ += newrank;
/*<               prof(2,nprof) = newrank >*/
		    prof[(nprof << 2) + 2] = newrank;
/*<               newlocf = newlocf + ldf * newrank >*/
		    newlocf += *ldf * newrank;
/*<               if (i .lt. nclim) then >*/
		    if (i__ < nclim) {
/*<                 jnode = supinds(i) >*/
			jnode = supinds[i__];
/*<                 k = 0 >*/
			k = 0;
/*<                 do while (supinds(j+k) .lt. jnode) >*/
			while(supinds[j + k] < jnode) {
/*<                   k = k + 1 >*/
			    ++k;
/*<                 end do >*/
			}
/*<               else >*/
		    } else {
/*<                 k = nrlim - j >*/
			k = nrlim - j;
/*<               end if >*/
		    }
/*<               prof(3,nprof) = k >*/
		    prof[(nprof << 2) + 3] = k;
/*<               prof(4,nprof) = jj >*/
		    prof[(nprof << 2) + 4] = jj;
/*<               j = j + k >*/
		    j += k;
/*<               newlocf = newlocf + k >*/
		    newlocf += k;
/*<               jj = jj + newrank >*/
		    jj += newrank;
/*<             end do >*/
		}
/*<           end if >*/
	    }
/*<         end if >*/
	}
/*<       end if >*/
    }
/*<       if (cinfo(bottom) .eq. 1) then >*/
    if (cinfo[bottom] == 1) {
/*<         locf = fptr >*/
	locf = *fptr;
/*<         j = rsuptr >*/
	j = *rsuptr;
/*<         do 40 i = csuptr, csuptr + min0(rank,ncols) - 1 >*/
	i__1 = *csuptr + min(rank,*ncols) - 1;
	for (i__ = *csuptr; i__ <= i__1; ++i__) {
/*<           jnode = supinds(i) >*/
	    jnode = supinds[i__];
/*<           jj = locf - j >*/
	    jj = locf - j;
/*<           ii = j >*/
	    ii = j;
/*<           do 30 k = aptrs(1,jnode), aptrs(1,jnode)+aptrs(2,jnode)-1 >*/
	    i__2 = aptrs[(jnode << 1) + 1] + aptrs[(jnode << 1) + 2] - 1;
	    for (k = aptrs[(jnode << 1) + 1]; k <= i__2; ++k) {
/*<             inode = ainds(k) >*/
		inode = ainds[k];
/*<             do while (supinds(ii) .lt. inode) >*/
		while(supinds[ii] < inode) {
/*<               ii = ii + 1 >*/
		    ++ii;
/*<             end do >*/
		}
/*<             wmem(jj+ii) = wmem(jj+ii) + avals(k) >*/
		wmem[jj + ii] += avals[k];
/*< 30        continue >*/
/* L30: */
	    }
/*<           locf = locf + ldf + diagproc >*/
	    locf = locf + *ldf + *diagproc;
/*<           j = j + diagproc >*/
	    j += *diagproc;
/*< 40      continue >*/
/* L40: */
	}
/*< 20      if (diagproc .eq. 1) then >*/
/* L20: */
	if (*diagproc == 1) {
/*<           call dpack(wmem(fptr),vbuf_s,rank,ldf,rank) >*/
	    dpack_(&wmem[*fptr], &vbuf_s__[1], &rank, ldf, &rank);
/*<           call dpotrf('l',rank,vbuf_s,rank,info) >*/
	    dpotrf_("l", &rank, &vbuf_s__[1], &rank, &info, (ftnlen)1);
/*<         if (info.gt.0) goto 1 >*/
	    if (info > 0) {
		goto L1;
	    }
/*<           if (myid .ne. mydown) then >*/
	    if (*myid != *mydown) {
/*<    >*/
		i__1 = rank * rank << 3;
		MPI_Send(&vbuf_s__[1], i__1, MPI_BYTE, *mydown, *myid, *comm);
/*<           end if >*/
	    }
/*<    >*/
	    i__1 = *nrows - rank;
	    dtrsm_("R", "L", "T", "N", &i__1, &rank, &c_b15, &vbuf_s__[1], &
		    rank, &wmem[*fptr + rank], ldf, (ftnlen)1, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
/*<           locf = fptr + rank >*/
	    locf = *fptr + rank;
/*<           uptr = 1 >*/
	    uptr = 1;
/*<           jj = 1 >*/
	    jj = 1;
/*<           do 50 i = 0, rank - 1 >*/
	    i__1 = rank - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
/*<             jnode = supinds(i+csuptr) >*/
		jnode = supinds[i__ + *csuptr];
/*<             call mydc(rank-i,vbuf_s(uptr),lvals(lptrs(1,jnode))) >*/
		i__2 = rank - i__;
		mydc_(&i__2, &vbuf_s__[uptr], &lvals[lptrs[jnode * 3 + 1]]);
/*<    >*/
		i__2 = *nrows - rank;
		myddc_(&i__2, &wmem[locf], &lvals[lptrs[jnode * 3 + 1] + rank 
			- i__], &hbuf_s__[jj]);
/*<             locf = locf + ldf  >*/
		locf += *ldf;
/*<             jj = jj + nrows - rank >*/
		jj = jj + *nrows - rank;
/*<             uptr = uptr + rank + 1 >*/
		uptr = uptr + rank + 1;
/*< 50        continue >*/
/* L50: */
	    }
/*<    >*/
	    i__1 = rank * (*nrows - rank) << 3;
	    MPI_Isend(&hbuf_s__[1], i__1, MPI_BYTE, *myright, *myid, *comm, req);
/*<           i = i + csuptr >*/
	    i__ += *csuptr;
/*<           j = rsuptr + rank >*/
	    j = *rsuptr + rank;
/*<           jnode = supinds(i) >*/
	    jnode = supinds[i__];
/*<           if (i .lt. nclim) then >*/
	    if (i__ < nclim) {
/*<             do while (supinds(j) .lt. jnode) >*/
		while(supinds[j] < jnode) {
/*<               j = j + 1 >*/
		    ++j;
/*<             end do >*/
		}
/*<           else >*/
	    } else {
/*<             j = nrlim >*/
		j = nrlim;
/*<           end if >*/
	    }
/*<           newlocf = locf + j - rsuptr - rank >*/
	    newlocf = locf + j - *rsuptr - rank;
/*<           nrows = nrlim - j  >*/
	    *nrows = nrlim - j;
/*<           uptr = fptr + j - rsuptr >*/
	    uptr = *fptr + j - *rsuptr;
/*<           rsuptr = j >*/
	    *rsuptr = j;
/*<           fptr = newlocf >*/
	    *fptr = newlocf;
/*<           csuptr = csuptr + rank >*/
	    *csuptr += rank;
/*<           ncols = ncols - rank >*/
	    *ncols -= rank;
/*<    >*/
	    i__1 = rank * *nrows << 3;
	    MPI_Isend(&hbuf_s__[1], i__1, MPI_BYTE, *mydown, *myid, *comm, &req[1]);

/*<           call mpi_waitall(4,req,statall,ierr) >*/
	    MPI_Waitall(4, req, statall);

/*<    >*/
	    if (*nrows > 0) {
		dsyrk_("L", "N", nrows, &rank, &c_b22, &wmem[uptr], ldf, &
			c_b15, &wmem[*fptr], ldf, (ftnlen)1, (ftnlen)1);
	    }
/*<         else >*/
	} else {
/*<           msgtype = MPI_ANY_TAG >*/
	    msgtype = MPI_ANY_TAG;
/*<    >*/

	    MPI_Recv(&vbuf_r__[1], buflen, MPI_BYTE, *myup, msgtype, *comm, mpistat);

/*<           call mpi_get_count(mpistat,MPI_BYTE,nbytes,ierr) >*/
	    MPI_Get_count(mpistat, MPI_BYTE, &nbytes);

/*<           msgtype = mpistat(MPI_TAG) >*/
	    msgtype = mpistat->MPI_TAG;
/*<           if (msgtype .ne. mydown) then >*/
	    if (msgtype != *mydown) {
/*<    >*/
		MPI_Send(&vbuf_r__[1], nbytes, MPI_BYTE, *mydown, msgtype, *comm);
/*<           end if >*/
	    }
/*<           if (ncols .gt. 0) then >*/
	    if (*ncols > 0) {
/*<    >*/
		dtrsm_("R", "L", "T", "N", nrows, &rank, &c_b15, &vbuf_r__[1],
			 &rank, &wmem[*fptr], ldf, (ftnlen)1, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
/*<           end if >*/
	    }
/*<           call dpack(wmem(fptr),hbuf_s,nrows,ldf,rank) >*/
	    dpack_(&wmem[*fptr], &hbuf_s__[1], nrows, ldf, &rank);
/*<    >*/
	    i__1 = *nrows * rank << 3;
	    MPI_Isend(&hbuf_s__[1], i__1, MPI_BYTE, *myright, *myid, *comm, req);
/*<           msgtype1 = MPI_ANY_TAG >*/
	    msgtype1 = MPI_ANY_TAG;
/*<    >*/
	    MPI_Irecv(&vbuf_r__[1], buflen, MPI_BYTE, *myup, msgtype1, *comm, &req[1]);
/*<           npending = 2 >*/
	    npending = 2;
/*<           nsent1 = 0 >*/
	    nsent1 = 0;
/*<           do while (npending .gt. 0) >*/
	    while(npending > 0) {
/*<             call mpi_waitany(4,req,msgid,mpistat,ierr) >*/
		MPI_Waitany(4, req, &msgid, mpistat);

/*<             call mpi_get_count(mpistat,MPI_BYTE,nbytes,ierr) >*/
		MPI_Get_count(mpistat, MPI_BYTE, &nbytes);

/*<             if (msgid .eq. 2 .and. nsent1 .eq. 0) then >*/
		if (msgid == 2 && nsent1 == 0) {
/*<               ldb = ishft(nbytes/rank,-3) >*/
		    ldb = lbit_shift(nbytes / rank, (ftnlen)-3);
/*<               msgtype1 = mpistat(MPI_TAG) >*/
		    msgtype1 = mpistat->MPI_TAG;
/*<               if (msgtype1 .ne. mydown) then >*/
		    if (msgtype1 != *mydown) {
/*<    >*/
			MPI_Isend(&vbuf_r__[1], nbytes, MPI_BYTE, *mydown, msgtype1, *comm, &req[2]);
/*<                 npending = npending + 1 >*/
			++npending;
/*<               end if >*/
		    }
/*<             nsent1 = 1 >*/
		    nsent1 = 1;
/*<             end if >*/
		}
/*<             npending = npending - 1 >*/
		--npending;
/*<           end do >*/
	    }
/*<           locf = fptr >*/
	    locf = *fptr;
/*<           k = lptrs(2,supinds(csuptr)) >*/
	    k = lptrs[supinds[*csuptr] * 3 + 2];
/*<           do 60 i = csuptr, csuptr + min0(ncols,rank) - 1 >*/
	    i__1 = *csuptr + min(*ncols,rank) - 1;
	    for (i__ = *csuptr; i__ <= i__1; ++i__) {
/*<             call mydc(k,wmem(locf),lvals(lptrs(1,supinds(i)))) >*/
		mydc_(&k, &wmem[locf], &lvals[lptrs[supinds[i__] * 3 + 1]]);
/*<             locf = locf + ldf >*/
		locf += *ldf;
/*< 60        continue >*/
/* L60: */
	    }
/*<           uptr = fptr >*/
	    uptr = *fptr;
/*<           newlocf = locf >*/
	    newlocf = locf;
/*<           csuptr = i >*/
	    *csuptr = i__;
/*<           if (rank .lt. ncols) then >*/
	    if (rank < *ncols) {
/*<             jj = 1 >*/
		jj = 1;
/*<             jnode = supinds(i) >*/
		jnode = supinds[i__];
/*<             j = rsuptr >*/
		j = *rsuptr;
/*<             do while (supinds(j) .lt. jnode) >*/
		while(supinds[j] < jnode) {
/*<               j = j + 1 >*/
		    ++j;
/*<             end do >*/
		}
/*<             newlocf = newlocf + j - rsuptr >*/
		newlocf = newlocf + j - *rsuptr;
/*<             uptr = uptr + j - rsuptr >*/
		uptr = uptr + j - *rsuptr;
/*<             fptr = newlocf >*/
		*fptr = newlocf;
/*<           end if >*/
	    }
/*<           nrows = nrlim - j >*/
	    *nrows = nrlim - j;
/*<           rsuptr = j >*/
	    *rsuptr = j;
/*<           ncols = ncols - rank >*/
	    *ncols -= rank;
/*<           nprof = 0 >*/
	    nprof = 0;
/*<           do while (i .lt. nclim) >*/
	    while(i__ < nclim) {
/*<             nprof = nprof + 1 >*/
		++nprof;
/*<             locf = newlocf >*/
		locf = newlocf;
/*<             prof(1,nprof) = newlocf >*/
		prof[(nprof << 2) + 1] = newlocf;
/*<             newrank = 0 >*/
		newrank = 0;
/*<             inode = supinds(j) >*/
		inode = supinds[j];
/*<             do while (i+newrank .lt. nclim) >*/
		while(i__ + newrank < nclim) {
/*<               if (supinds(i+newrank) .gt. inode) goto 90 >*/
		    if (supinds[i__ + newrank] > inode) {
			goto L90;
		    }
/*<               newrank = newrank + 1 >*/
		    ++newrank;
/*<             end do >*/
		}
/*< 90          i = i + newrank  >*/
L90:
		i__ += newrank;
/*<             prof(2,nprof) = newrank >*/
		prof[(nprof << 2) + 2] = newrank;
/*<             newlocf = newlocf + ldf * newrank >*/
		newlocf += *ldf * newrank;
/*<             if (i .lt. nclim) then >*/
		if (i__ < nclim) {
/*<               jnode = supinds(i) >*/
		    jnode = supinds[i__];
/*<               k = 0 >*/
		    k = 0;
/*<               do while (supinds(j+k) .lt. jnode) >*/
		    while(supinds[j + k] < jnode) {
/*<                 k = k + 1 >*/
			++k;
/*<               end do >*/
		    }
/*<             else >*/
		} else {
/*<               k = nrlim - j >*/
		    k = nrlim - j;
/*<             end if >*/
		}
/*<             prof(3,nprof) = k >*/
		prof[(nprof << 2) + 3] = k;
/*<             prof(4,nprof) = jj >*/
		prof[(nprof << 2) + 4] = jj;
/*<    >*/
		i__1 = nrlim - j;
		dgemm_("N", "T", &i__1, &newrank, &rank, &c_b22, &wmem[uptr], 
			ldf, &vbuf_r__[jj], &ldb, &c_b15, &wmem[locf], ldf, (
			ftnlen)1, (ftnlen)1);
/*<             j = j + k >*/
		j += k;
/*<             newlocf = newlocf + k >*/
		newlocf += k;
/*<             jj = jj + newrank >*/
		jj += newrank;
/*<             uptr = uptr + k >*/
		uptr += k;
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<         if (ncols .lt. 0) then >*/
	if (*ncols < 0) {
/*<           ncols = 0 >*/
	    *ncols = 0;
/*<           csuptr = nclim >*/
	    *csuptr = nclim;
/*<         end if >*/
	}
/*<         if (nrows .lt. 0 .or. ncols .eq. 0) then >*/
	if (*nrows < 0 || *ncols == 0) {
/*<           nrows = 0 >*/
	    *nrows = 0;
/*<           rsuptr = nrlim >*/
	    *rsuptr = nrlim;
/*<         end if >*/
	}
/*<       else >*/
    } else {
/*<         if (diagproc .eq. 1) then >*/
	if (*diagproc == 1) {
/*<           msgtype = MPI_ANY_TAG  >*/
	    msgtype = MPI_ANY_TAG;
/*<    >*/
	    MPI_Recv(&hbuf_r__[1], buflen, MPI_BYTE, *myleft, msgtype, *comm, mpistat);
	    
/*<           call mpi_get_count(mpistat,MPI_BYTE,nbytes,ierr) >*/
	    MPI_Get_count(mpistat, MPI_BYTE, &nbytes);
	    
/*<           msgtype = mpistat(MPI_TAG) >*/
	    msgtype = mpistat->MPI_TAG;
/*<           ldb = ishft(nbytes/rank,-3) >*/
	    ldb = lbit_shift(nbytes / rank, (ftnlen)-3);
/*<           uptr = 1 + ldb - nrows >*/
	    uptr = ldb + 1 - *nrows;
/*<           if (myid .ne. mydown) then >*/
	    if (*myid != *mydown) {
/*<             if (ldb .eq. nrows) then >*/
		if (ldb == *nrows) {
/*<    >*/
		    MPI_Isend(&hbuf_r__[1], nbytes, MPI_BYTE, *mydown, *myid, *comm, req);
/*<             else >*/
		} else {
/*<               call dpack(hbuf_r(uptr),vbuf_s,nrows,ldb,rank) >*/
		    dpack_(&hbuf_r__[uptr], &vbuf_s__[1], nrows, &ldb, &rank);
/*<    >*/
		    i__1 = *nrows * rank << 3;
		    MPI_Isend(&vbuf_s__[1], i__1, MPI_BYTE, *mydown, *myid, *comm, req);
/*<             end if >*/
		}
/*<           end if >*/
	    }
/*<           if (myright .ne. msgtype) then >*/
	    if (*myright != msgtype) {
/*<    >*/
		MPI_Isend(&hbuf_r__[1], nbytes, MPI_BYTE, *myright, msgtype, *comm, &req[1]);
/*<           end if >*/
	    }
/*<           call mpi_waitall(4,req,statall,ierr) >*/
	    MPI_Waitall(4, req, statall);
/*<    >*/
	    if (*nrows > 0) {
		dsyrk_("L", "N", nrows, &rank, &c_b22, &hbuf_r__[uptr], &ldb, 
			&c_b15, &wmem[*fptr], ldf, (ftnlen)1, (ftnlen)1);
	    }
/*<         else >*/
	} else {
/*<           msgtype1 = MPI_ANY_TAG >*/
	    msgtype1 = MPI_ANY_TAG;
/*<           msgtype2 = MPI_ANY_TAG >*/
	    msgtype2 = MPI_ANY_TAG;
/*<    >*/
	    MPI_Irecv(&hbuf_r__[1], buflen, MPI_BYTE, *myleft, msgtype1, *comm, req);
/*<    >*/
	    MPI_Irecv(&vbuf_r__[1], buflen, MPI_BYTE, *myup, msgtype2, *comm, &req[1]);
/*<           npending = 2 >*/
	    npending = 2;
/*<         nsent1 = 0 >*/
	    nsent1 = 0;
/*<         nsent2 = 0 >*/
	    nsent2 = 0;
/*<           do while (npending .ne. 0) >*/
	    while(npending != 0) {
/*<             call mpi_waitany(4,req,msgid,mpistat,ierr) >*/
		MPI_Waitany(4, req, &msgid, mpistat);
		
/*<             call mpi_get_count(mpistat,MPI_BYTE,nbytes,ierr) >*/
		MPI_Get_count(mpistat, MPI_BYTE, &nbytes);

/*<             npending = npending - 1 >*/
		--npending;
/*<             if (msgid .eq. 1 .and. nsent1 .eq. 0) then >*/
		if (msgid == 1 && nsent1 == 0) {
/*<               lda = ishft(nbytes/rank,-3) >*/
		    lda = lbit_shift(nbytes / rank, (ftnlen)-3);
/*<               msgtype1 = mpistat(MPI_TAG) >*/
		    msgtype1 = mpistat->MPI_TAG;
/*<               if (msgtype1 .ne. myright) then >*/
		    if (msgtype1 != *myright) {
/*<    >*/
			MPI_Isend(&hbuf_r__[1], nbytes, MPI_BYTE, *myright, msgtype1, *comm, &req[2]);

/*<                 npending = npending + 1 >*/
			++npending;
/*<               end if >*/
		    }
/*<             nsent1 = 1 >*/
		    nsent1 = 1;
/*<             end if >*/
		}
/*<             if (msgid .eq. 2 .and. nsent2 .eq. 0) then >*/
		if (msgid == 2 && nsent2 == 0) {
/*<               ldb = ishft(nbytes/rank,-3) >*/
		    ldb = lbit_shift(nbytes / rank, (ftnlen)-3);
/*<               msgtype2 = mpistat(MPI_TAG) >*/
		    msgtype2 = mpistat->MPI_TAG;
/*<               if (msgtype2 .ne. mydown) then >*/
		    if (msgtype2 != *mydown) {
/*<    >*/
			MPI_Isend(&vbuf_r__[1], nbytes, MPI_BYTE, *mydown, msgtype2, *comm, &req[3]);

/*<                 npending = npending + 1 >*/
			++npending;
/*<               end if >*/
		    }
/*<               nsent2 = 1 >*/
		    nsent2 = 1;
/*<             end if >*/
		}
/*<           end do >*/
	    }
/*<           uptr = 1 + lda - nrows >*/
	    uptr = lda + 1 - *nrows;
/*<           j = nrows >*/
	    j = *nrows;
/*<           do 130 i = 1, nprof >*/
	    i__1 = nprof;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/*<    >*/
		dgemm_("N", "T", &j, &prof[(i__ << 2) + 2], &rank, &c_b22, &
			hbuf_r__[uptr], &lda, &vbuf_r__[prof[(i__ << 2) + 4]],
			 &ldb, &c_b15, &wmem[prof[(i__ << 2) + 1]], ldf, (
			ftnlen)1, (ftnlen)1);
/*<             j = j - prof(3,i) >*/
		j -= prof[(i__ << 2) + 3];
/*<             uptr = uptr + prof(3,i) >*/
		uptr += prof[(i__ << 2) + 3];
/*< 130       continue >*/
/* L130: */
	    }
/*<         end if >*/
	}
/*<       end if >*/
    }
/*<       if (tptrs(3,node) .eq. 0) goto 1000  >*/
    if (tptrs[node * 3 + 3] == 0) {
	goto L1000;
    }
/*<       if (nptr .eq. 0) return  >*/
    if (*nptr == 0) {
	return 0;
    }
/*<       if (sup(tptrs(3,node)) .ne. node) goto 1000 >*/
    if (sup[tptrs[node * 3 + 3]] != node) {
	goto L1000;
    }
/*<       if (tptrs(2,node) .eq. 1) then >*/
    if (tptrs[node * 3 + 2] == 1) {
/*<         newsup = 1 >*/
	newsup = 1;
/*<         uptr = fptr >*/
	uptr = *fptr;
/*<         ncols_u = nclim >*/
	ncols_u__ = nclim;
/*<         csuptr_u = csuptr >*/
	csuptr_u__ = *csuptr;
/*<         rsuptr_u = rsuptr >*/
	rsuptr_u__ = *rsuptr;
/*<         stakptr = stakptr - 1 >*/
	--(*stakptr);
/*<         ncols = stak(2,stakptr) >*/
	*ncols = stak[*stakptr * 3 + 2];
/*<         nrows = stak(3,stakptr) >*/
	*nrows = stak[*stakptr * 3 + 3];
/*<         i = tptrs(3,node) + 2 >*/
	i__ = tptrs[node * 3 + 3] + 2;
/*<         csuptr = sup(i) >*/
	*csuptr = sup[i__];
/*<         rsuptr = csuptr + sup(i+1) >*/
	*rsuptr = *csuptr + sup[i__ + 1];
/*<         fptr = wmemsize - ldf * (ncols - 1) - nrows + 1 >*/
	*fptr = *wmemsize - *ldf * (*ncols - 1) - *nrows + 1;
/*<         i = rsuptr_u >*/
	i__ = rsuptr_u__;
/*<         j = nrlim >*/
	j = nrlim;
/*<         k = rsuptr >*/
	k = *rsuptr;
/*<         l = rsuptr + nrows >*/
	l = *rsuptr + *nrows;
/*<         locptr = -1 >*/
	locptr = -1;
/*<         do while (i .lt. j .and. k .lt. l) >*/
	while(i__ < j && k < l) {
/*<           locptr = locptr + 2 >*/
	    locptr += 2;
/*<           do while (i .lt. j .and. supinds(i) .eq. supinds(k)) >*/
	    while(i__ < j && supinds[i__] == supinds[k]) {
/*<             i = i + 1 >*/
		++i__;
/*<             k = k + 1 >*/
		++k;
/*<           end do >*/
	    }
/*<           locinds(locptr) = k - 1 >*/
	    locinds[locptr] = k - 1;
/*<           if (i .lt. j) then >*/
	    if (i__ < j) {
/*<             do while (supinds(i) .ne. supinds(k)) >*/
		while(supinds[i__] != supinds[k]) {
/*<               k = k + 1 >*/
		    ++k;
/*<             end do >*/
		}
/*<             locinds(locptr+1) = k - 1 >*/
		locinds[locptr + 1] = k - 1;
/*<           else >*/
	    } else {
/*<             locinds(locptr+1) = l - 1 >*/
		locinds[locptr + 1] = l - 1;
/*<           end if >*/
	    }
/*<         end do >*/
	}
/*<         locf = fptr >*/
	locf = *fptr;
/*<         nclim = ncols + csuptr >*/
	nclim = *ncols + *csuptr;
/*<         nrlim = nrows + rsuptr >*/
	nrlim = *nrows + *rsuptr;
/*<         jj = csuptr >*/
	jj = *csuptr;
/*<         ii = rsuptr >*/
	ii = *rsuptr;
/*<         do while (jj .lt. nclim .and. csuptr_u .lt. ncols_u) >*/
	while(jj < nclim && csuptr_u__ < ncols_u__) {
/*<           do while (supinds(jj) .lt. supinds(csuptr_u)) >*/
	    while(supinds[jj] < supinds[csuptr_u__]) {
/*<             do while (supinds(ii) .lt. supinds(jj)) >*/
		while(supinds[ii] < supinds[jj]) {
/*<               ii = ii + 1 >*/
		    ++ii;
/*<               locf = locf + 1 >*/
		    ++locf;
/*<             end do >*/
		}
/*<             do 220 i = 0, nrlim - ii - 1 >*/
		i__1 = nrlim - ii - 1;
		for (i__ = 0; i__ <= i__1; ++i__) {
/*<               wmem(locf+i) = 0.d0 >*/
		    wmem[locf + i__] = 0.;
/*< 220         continue >*/
/* L220: */
		}
/*<             jj = jj + 1 >*/
		++jj;
/*<             locf = locf + ldf >*/
		locf += *ldf;
/*<           end do >*/
	    }
/*<           do while (supinds(ii) .lt. supinds(jj)) >*/
	    while(supinds[ii] < supinds[jj]) {
/*<             ii = ii + 1 >*/
		++ii;
/*<             locf = locf + 1 >*/
		++locf;
/*<           end do >*/
	    }
/*<           if (locf .eq. uptr) goto 1000 >*/
	    if (locf == uptr) {
		goto L1000;
	    }
/*<           i = ii >*/
	    i__ = ii;
/*<           locf = locf - ii >*/
	    locf -= ii;
/*<           j = uptr >*/
	    j = uptr;
/*<           do 230 l = 1, locptr, 2 >*/
	    i__1 = locptr;
	    for (l = 1; l <= i__1; l += 2) {
/*<             j = j - i >*/
		j -= i__;
/*<             if (locf .eq. j) goto 310 >*/
		if (locf == j) {
		    goto L310;
		}
/*<             do 240 k = i, locinds(l) >*/
		i__2 = locinds[l];
		for (k = i__; k <= i__2; ++k) {
/*<               wmem(locf+k) = wmem(j+k) >*/
		    wmem[locf + k] = wmem[j + k];
/*< 240         continue >*/
/* L240: */
		}
/*<             j = j + k >*/
		j += k;
/*<             do 250 i = k, locinds(l+1) >*/
		i__2 = locinds[l + 1];
		for (i__ = k; i__ <= i__2; ++i__) {
/*<               wmem(locf+i) = 0.d0 >*/
		    wmem[locf + i__] = 0.;
/*< 250         continue >*/
/* L250: */
		}
/*< 230       continue >*/
/* L230: */
	    }
/*< 310       locf = locf + ldf + ii >*/
L310:
	    locf = locf + *ldf + ii;
/*<           uptr = uptr + ldf >*/
	    uptr += *ldf;
/*<           csuptr_u = csuptr_u + 1 >*/
	    ++csuptr_u__;
/*<           jj = jj + 1 >*/
	    ++jj;
/*<           do while (supinds(rsuptr_u) .lt. supinds(csuptr_u)) >*/
	    while(supinds[rsuptr_u__] < supinds[csuptr_u__]) {
/*<             rsuptr_u = rsuptr_u + 1 >*/
		++rsuptr_u__;
/*<             uptr = uptr + 1 >*/
		++uptr;
/*<           end do >*/
	    }
/*<         end do >*/
	}
/*<         do while (jj .lt. nclim) >*/
	while(jj < nclim) {
/*<           do while (supinds(ii) .lt. supinds(jj)) >*/
	    while(supinds[ii] < supinds[jj]) {
/*<             ii = ii + 1 >*/
		++ii;
/*<             locf = locf + 1 >*/
		++locf;
/*<           end do >*/
	    }
/*<           do 210 i = 0, nrlim - ii - 1 >*/
	    i__1 = nrlim - ii - 1;
	    for (i__ = 0; i__ <= i__1; ++i__) {
/*<             wmem(locf+i) = 0.d0 >*/
		wmem[locf + i__] = 0.;
/*< 210       continue >*/
/* L210: */
	    }
/*<           jj = jj + 1 >*/
	    ++jj;
/*<           locf = locf + ldf >*/
	    locf += *ldf;
/*<         end do >*/
	}
/*<         goto 1000 >*/
	goto L1000;
/*<       end if >*/
    }
/*<       return >*/
    return 0;
/*< 1  >*/
L1:
	printf("%d: non positive definite matrix found at %d", *myid, bottom);
	MPI_Abort(*comm, 27);
/*<       call mpi_abort(comm,27,ierr) >*/
/*<       end  >*/
    return 0;
} /* parelimh1_ */

