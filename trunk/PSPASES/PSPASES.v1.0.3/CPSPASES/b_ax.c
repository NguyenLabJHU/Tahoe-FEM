/* $Id: b_ax.c,v 1.3 2004-12-13 00:27:45 paklein Exp $ */
/* b_ax.f -- translated by f2c (version 20030320).
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
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_b22 = 524288;
static integer c__13 = 13;
static integer c__21 = 21;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   b_ax.f                                                                  +/ */
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
/* /+                                                                           +/em */
/* /+***************************************************************************+/ */
/* /+ $Id: b_ax.c,v 1.3 2004-12-13 00:27:45 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
/* Subroutine */ int db_ax_(integer *n, integer *rowdista, integer *rowdistb,
	 integer *nrhs, integer *aptrs, integer *ainds, doublereal *avals, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, integer *
	myid, integer *pp, doublereal *emax, MPI_Comm *comm)

/* dummy arguments for f2c
	, doublereal *tx, 
	doublereal *tb, doublereal *bmax)
*/
{
	/* debugging */
#ifdef __DO_DEBUG__
	FILE* fp = NULL;
	char file[] = "rankN";
	char ints[] = "0123456789";
	int dummy;
#endif

    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, tx_dim1, tx_offset, tb_dim1, 
	    tb_offset, bmax_dim1, bmax_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j, k;
    integer mynnodesa, mynnodesb, ofs;
    doublereal err;
    integer ierr;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       double precision zero >*/
/* -*- fortran -*- */

/* double precision functions */

/*<       double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       parameter(zero=0.d0) >*/
/*<       integer N,rowdista(0:*),rowdistb(0:*),nrhs,aptrs(2,0:*),pp >*/
/*<       integer ldx,ldb >*/
/*<       integer ainds(*),comm,ierr,i,j,k,mynnodesb,mynnodesa,ofs,myid >*/
/*<       double precision avals(*),b(0:ldb-1,*),x(0:ldx-1,*),err,emax >*/
/*      double precision, allocatable :: tx(:,:),tb(:,:),bmax(:,:) */

	doublereal *tx, *tb, *bmax; 

	/* debugging */
#ifdef __DO_DEBUG__
	file[4] = ints[*myid];
	fp = fopen(file, "w");
	fprintf(fp, "   n = %d\n", *n);
	fprintf(fp, "nrhs = %d\n", *nrhs);
#endif

/*<       double precision tx,tb,bmax >*/
/*<       dimension tx(0:N-1,nrhs) >*/
/*<       dimension tb(0:N-1,nrhs) >*/
/*<       dimension bmax(0:N-1,nrhs) >*/
/*      allocate(tx(0:N-1,nrhs),stat=k) */
	tx = (doublereal*) malloc((*n)*(*nrhs)*sizeof(doublereal));
/*<       if(k.ne.0) then >*/
    /* Function Body */
    if (!tx) {
/*<         print *,'memory allocation failure' >*/
		printf("memory allocation failure\n");
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }

/*      allocate(tb(0:N-1,nrhs),stat=k) */
	tb = (doublereal*) malloc((*n)*(*nrhs)*sizeof(doublereal));
/*<       if(k.ne.0) then >*/
    if (!tb) {
/*<         print *,'memory allocation failure' >*/
		printf("memory allocation failure\n");
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }

/*     allocate(bmax(0:N-1,nrhs),stat=k) */
	bmax = (doublereal*) malloc((*n)*(*nrhs)*sizeof(doublereal));
/*<       if(k.ne.0) then >*/
    if (!bmax) {
/*<         print *,'memory allocation failure' >*/
		printf("memory allocation failure\n");
/*<         call mpi_abort(comm,0,ierr) >*/
		MPI_Abort(*comm, 0);
/*<       end if >*/
    }

	/* zero the arrays */
	psp_clear(  tb, (*n)*(*nrhs));
	psp_clear(  tx, (*n)*(*nrhs));
	psp_clear(bmax, (*n)*(*nrhs));

    /* Parameter adjustments */
    bmax_dim1 = *n - 1 - 0 + 1;
    bmax_offset = 0 + bmax_dim1;
    bmax -= bmax_offset;

    tb_dim1 = *n - 1 - 0 + 1;
    tb_offset = 0 + tb_dim1;
    tb -= tb_offset;

    tx_dim1 = *n - 1 - 0 + 1;
    tx_offset = 0 + tx_dim1;
    tx -= tx_offset;

    --aptrs;
    --ainds;
    --avals;
    b_dim1 = *ldb - 1 - 0 + 1;
    b_offset = 0 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx - 1 - 0 + 1;
    x_offset = 0 + x_dim1;
    x -= x_offset;

/*<       mynnodesb = rowdistb(myid+1)-rowdistb(myid) >*/
    mynnodesb = rowdistb[*myid + 1] - rowdistb[*myid];
/*<       k = rowdistb(myid) >*/
    k = rowdistb[*myid];
/*      gather b. */
/*<       do j=1,nrhs >*/
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
/*<         do i=0,N-1 >*/
	i__2 = *n - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           bmax(i,j) = zero >*/
	    bmax[i__ + j * bmax_dim1] = 0.;
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       do j=1,nrhs >*/
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
/*<         do i=0,mynnodesb-1 >*/
	i__2 = mynnodesb - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           bmax(i+k,j) = b(i,j) >*/
	    bmax[i__ + k + j * bmax_dim1] = b[i__ + j * b_dim1];
/*<         end do >*/
	}
/*<       end do >*/
    }

	/* debugging */
#ifdef __DO_DEBUG__
	fprintf(fp, "b unreduced:\n");
	psp_dump(fp, bmax + bmax_offset, (*n)*(*nrhs));
	fprintf(fp, "\n");
#endif

/*<       call staged_mpirds(bmax,tb,N*nrhs,0,0,comm) >*/
    i__1 = *n * *nrhs;
    staged_mpirds__(&bmax[bmax_offset], &tb[tb_offset], &i__1, &c__0, &c__0, 
	    comm);

/*      gather x. */
/*<       do j=1,nrhs >*/
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
/*<         do i=0,N-1 >*/
	i__2 = *n - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           bmax(i,j) = zero >*/
	    bmax[i__ + j * bmax_dim1] = 0.;
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       do j=1,nrhs >*/
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
/*<         do i=0,mynnodesb-1 >*/
	i__2 = mynnodesb - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           bmax(i+k,j) = x(i,j) >*/
	    bmax[i__ + k + j * bmax_dim1] = x[i__ + j * x_dim1];
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       call staged_mpirds(bmax,tx,N*nrhs,0,0,comm) >*/
    i__1 = *n * *nrhs;
    staged_mpirds__(&bmax[bmax_offset], &tx[tx_offset], &i__1, &c__0, &c__0, 
	    comm);
/*      call b-Ax routine. */
/*<       do j=1,nrhs >*/
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
/*<         do i=0,N-1 >*/
	i__2 = *n - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<           bmax(i,j) = zero >*/
	    bmax[i__ + j * bmax_dim1] = 0.;
/*<         end do >*/
	}
/*<       end do >*/
    }
/*<       mynnodesa = rowdista(myid+1)-rowdista(myid) >*/
    mynnodesa = rowdista[*myid + 1] - rowdista[*myid];
/*<       ofs = rowdista(myid) >*/
    ofs = rowdista[*myid];

	/* debugging */
#ifdef __DO_DEBUG__
	fprintf(fp, "b =\n");	
	psp_dump(fp, tb + tb_offset, (*n)*(*nrhs));
	fprintf(fp, "\n");

	fprintf(fp, "x = \n");
	psp_dump(fp, tx + tx_offset, (*n)*(*nrhs));
	fprintf(fp, "\n");
#endif

/*<       call b_ax(aptrs,ainds,avals,tb,tx,N,mynnodesa,ofs,nrhs,bmax) >*/
    b_ax__(&aptrs[1], &ainds[1], &avals[1], &tb[tb_offset], &tx[tx_offset], n,
	     &mynnodesa, &ofs, nrhs, &bmax[bmax_offset]);

	/* debugging */
#ifdef __DO_DEBUG__
	fprintf(fp, "local b-ax = \n");
	psp_dump(fp, bmax + bmax_offset, (*n)*(*nrhs));
	fprintf(fp, "\n");
#endif

/*<       call staged_mpirds(bmax,tx,N*nrhs,1,0,comm) >*/
    i__1 = *n * *nrhs;
    staged_mpirds__(&bmax[bmax_offset], &tx[tx_offset], &i__1, &c__1, &c__0, 
	    comm);

	/* debugging */
#ifdef __DO_DEBUG__
	fprintf(fp, "b - ax = \n");
	psp_dump(fp, tx + tx_offset, (*n)*(*nrhs));
	fprintf(fp, "\n");
#endif

/*<       if (myid.eq.0) then >*/
    if (*myid == 0) {
/*<         emax = 0.d0 >*/
	*emax = 0.;
/*<         do j=1,nrhs >*/
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
/*<           do i=0,N-1 >*/
	    i__2 = *n - 1;
	    for (i__ = 0; i__ <= i__2; ++i__) {
/*<             err = dabs(tx(i,j)) >*/
		err = tx[i__ + j * tx_dim1];
		if (err < 0.0) err = -err;
/*<             if(err.gt.emax) emax = err >*/
		if (err > *emax) {
		    *emax = err;
		}

/* debugging */
#ifdef __DO_DEBUG__
if (err > *emax)
	dummy = 1;
else
	dummy = 0;
fprintf(fp, "%20.14le %20.14le %d\n", *emax, err, dummy);		
#endif

/*<           end do >*/
	    }
/*<         end do >*/
	}
/*        print *,'max |B - AX| = ',emax */
/*<       end if >*/
    }
/*      deallocate(tx) */
/*      deallocate(tb) */
/*      deallocate(bmax) */

	/* restore offsets */
	tx += tx_offset;
	tb += tb_offset;
	bmax += bmax_offset;

	free(tx);
	free(tb);
	free(bmax);

/*<       end  >*/

	/* debugging */
#ifdef __DO_DEBUG__
	fclose(fp);
#endif

    return 0;
} /* db_ax__ */

/*<       subroutine b_ax(aptrs,ainds,avals,b,x,N,mynnodes,ofs,nrhs,ax) >*/
/* Subroutine */ int b_ax__(integer *aptrs, integer *ainds, doublereal *avals,
	 doublereal *b, doublereal *x, integer *n, integer *mynnodes, integer 
	*ofs, integer *nrhs, doublereal *ax)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, ax_dim1, ax_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    integer i__, j, k, ig, node, ivptr, ivsiz;

/*<       double precision avals(*),b(0:N-1,*),x(0:N-1,*),ax(0:N-1,*) >*/
/*<       integer aptrs(2,0:*),ainds(*),mynnodes,ofs >*/
/*<       do i=0,mynnodes-1 >*/
    /* Parameter adjustments */
    --aptrs;
    --ainds;
    --avals;
    ax_dim1 = *n - 1 - 0 + 1;
    ax_offset = 0 + ax_dim1;
    ax -= ax_offset;
    x_dim1 = *n - 1 - 0 + 1;
    x_offset = 0 + x_dim1;
    x -= x_offset;
    b_dim1 = *n - 1 - 0 + 1;
    b_offset = 0 + b_dim1;
    b -= b_offset;

    /* Function Body */
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         ivptr = aptrs(1,i) >*/
	ivptr = aptrs[(i__ << 1) + 1];
/*<         ivsiz = aptrs(2,i) >*/
	ivsiz = aptrs[(i__ << 1) + 2];
/*<         ig = i+ofs >*/
	ig = i__ + *ofs;
/*<         if(ivsiz.ne.0) then >*/
	if (ivsiz != 0) {
/*<           do while(ainds(ivptr).lt.ig) >*/
	    while(ainds[ivptr] < ig) {
/*<             ivptr = ivptr+1 >*/
		++ivptr;
/*<             ivsiz = ivsiz-1 >*/
		--ivsiz;
/*<           end do >*/
	    }
/*<           if(ainds(ivptr).eq.ig) then >*/
	    if (ainds[ivptr] == ig) {
/*<             do k=1,nrhs >*/
		i__2 = *nrhs;
		for (k = 1; k <= i__2; ++k) {
/*<               ax(ig,k) = ax(ig,k) - avals(ivptr)*x(ig,k) >*/
		    ax[ig + k * ax_dim1] -= avals[ivptr] * x[ig + k * x_dim1];
/*<               ax(ig,k) = ax(ig,k) + b(ig,k) >*/
		    ax[ig + k * ax_dim1] += b[ig + k * b_dim1];
/*<             end do >*/
		}
/*<             ivptr = ivptr+1 >*/
		++ivptr;
/*<             ivsiz = ivsiz-1 >*/
		--ivsiz;
/*<           end if >*/
	    }
/*<           do k=1,nrhs >*/
	    i__2 = *nrhs;
	    for (k = 1; k <= i__2; ++k) {
/*<           do j=0,ivsiz-1 >*/
		i__3 = ivsiz - 1;
		for (j = 0; j <= i__3; ++j) {
/*<             node = ainds(ivptr+j) >*/
		    node = ainds[ivptr + j];
/*<             ax(node,k) = ax(node,k) - avals(ivptr+j)*x(ig,k) >*/
		    ax[node + k * ax_dim1] -= avals[ivptr + j] * x[ig + k * 
			    x_dim1];
/*<             ax(ig,k) = ax(ig,k) - avals(ivptr+j)*x(node,k) >*/
		    ax[ig + k * ax_dim1] -= avals[ivptr + j] * x[node + k * 
			    x_dim1];
/*<           end do >*/
		}
/*<           end do >*/
	    }
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       end >*/
    return 0;
} /* b_ax__ */

/*<       subroutine staged_mpirds(sbuf,rbuf,size,opt,root,comm) >*/
/* Subroutine */ int staged_mpirds__(doublereal *sbuf, doublereal *rbuf, 
	integer *size, integer *opt, integer *root, MPI_Comm *comm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    integer ierr;
    integer remains;

/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer eltlimit >*/
/* -*- fortran -*- */

/* double precision functions */

/*<       double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       parameter(eltlimit=512*1024) ! corresponds to 4 MB  >*/
/*<       integer size,opt,root,comm >*/
/*<       double precision sbuf(*),rbuf(*) >*/
/*<       integer i,remains,ierr >*/
/*<       do i=1,size-eltlimit+1,eltlimit >*/
    /* Parameter adjustments */
    --rbuf;
    --sbuf;

    /* Function Body */
    i__1 = *size - 524287;
    for (i__ = 1; i__ <= i__1; i__ += 524288) {
/*<         if(opt.eq.0) then >*/
	if (*opt == 0) {
/*<    >*/
/*	    mpi_allreduce__(&sbuf[i__], &rbuf[i__], &c_b22, &c__13, &c__21, comm, &ierr); */
		MPI_Allreduce(&sbuf[i__], &rbuf[i__], c_b22, MPI_DOUBLE, MPI_SUM, *comm);
/*<         else >*/
	} else {
/*<    >*/
/*	    mpi_reduce__(&sbuf[i__], &rbuf[i__], &c_b22, &c__13, &c__21, root, comm, &ierr); */
		MPI_Reduce(&sbuf[i__], &rbuf[i__], c_b22, MPI_DOUBLE, MPI_SUM, *root, *comm);
/*<         end if >*/
	}
/*<       end do >*/
    }
/*<       remains = size-i+1 >*/
    remains = *size - i__ + 1;
/*<       if(remains.gt.0) then >*/
    if (remains > 0) {
/*<         if(opt.eq.0) then >*/
	if (*opt == 0) {
/*<    >*/
/*	    mpi_allreduce__(&sbuf[i__], &rbuf[i__], &remains, &c__13, &c__21, comm, &ierr); */
		MPI_Allreduce(&sbuf[i__], &rbuf[i__], remains, MPI_DOUBLE, MPI_SUM, *comm);
/*<         else >*/
	} else {
/*<    >*/
/*	    mpi_reduce__(&sbuf[i__], &rbuf[i__], &remains, &c__13, &c__21, root, comm, &ierr); */
		MPI_Reduce(&sbuf[i__], &rbuf[i__], remains, MPI_DOUBLE, MPI_SUM, *root, *comm);
/*<         end if >*/
	}
/*<       end if >*/
    }
/*<       end >*/
    return 0;
} /* staged_mpirds__ */
