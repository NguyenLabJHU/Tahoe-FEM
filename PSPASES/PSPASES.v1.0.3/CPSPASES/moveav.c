/* moveav.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublecomplex mpi_bottom__;
} mpi_bottom__;

#define mpi_bottom__1 mpi_bottom__

struct {
    doublecomplex mpi_argv_null__;
} mpi_argv_null__;

#define mpi_argv_null__1 mpi_argv_null__

struct {
    doublecomplex mpi_argvs_null__;
} mpi_argvs_null__;

#define mpi_argvs_null__1 mpi_argvs_null__

struct {
    doublecomplex mpi_errcodes_ignore__;
} mpi_errcodes_ignore__;

#define mpi_errcodes_ignore__1 mpi_errcodes_ignore__

struct {
    doublecomplex mpi_status_ignore__;
} mpi_status_ignore__;

#define mpi_status_ignore__1 mpi_status_ignore__

struct {
    doublecomplex mpi_statuses_ignore__;
} mpi_statuses_ignore__;

#define mpi_statuses_ignore__1 mpi_statuses_ignore__

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__11 = 11;
static integer c__13 = 13;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   moveav.f                                                                +/ */
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
/* /+ $Id: moveav.c,v 1.1 2004-12-12 21:50:35 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
/* Subroutine */ int moveav_(integer *n, integer *dd, integer *pp, integer *
	lgblk, integer *myid, integer *rowdista, integer *mynnodes, integer *
	order, integer *aptrs, integer *ainds, doublereal *avals, doublereal *
	pavals, integer *wrkint, integer *maxnzpercol, integer *ranmasks, 
	integer *comm, integer *gorder, integer *whichsnode, integer *tainds, 
	doublereal *sendvals, integer *i2, integer *nsend2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer lbit_shift(integer, integer), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle(void);

    /* Local variables */
    static integer beginrow, i__, j, k, l, m;
    extern /* Subroutine */ int mpi_abort__(integer *, integer *, integer *), 
	    ikeysortf_(integer *, integer *, integer *);
    static integer is1, col, ppc, ppg, ppr, row, ierr, proc, prcv, pscv, psdv,
	     prdv;
    extern /* Subroutine */ int mpi_alltoall__(integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *);
    static integer ptr_sendvals__, ptr_c__, nsend, ptr_r__;
    extern /* Subroutine */ int mpi_alltoallv__(doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *);
    static integer bmaskc, bmaskr, fptr_r__;
    extern /* Subroutine */ int mpi_allgatherv__(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    static integer itainds, pgrsize;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };


/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       integer*8 loc >*/
/* -*- fortran -*- */

/* Copyright (c) 2001-2002 The Trustees of Indiana University. */
/*                         All rights reserved. */
/* Copyright (c) 1998-2001 University of Notre Dame. */
/*                         All rights reserved. */
/* Copyright (c) 1994-1998 The Ohio State University. */
/*                         All rights reserved. */

/* This file is part of the LAM/MPI software package.  For license */
/* information, see the LICENSE file in the top level directory of the */
/* LAM/MPI source distribution. */


/*  $Id: moveav.c,v 1.1 2004-12-12 21:50:35 paklein Exp $ */

/* 	Function:	- LAM/MPI F77 header file */

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* Do ***not*** copy this file to the directory where your Fortran */
/* fortran application is compiled unless it is absolutely necessary!  Most */
/* modern Fortran compilers now support the -I command line flag, which */
/* tells the compiler where to find .h files (specifically, this one).  For */
/* example: */

/*      unix% mpif77 foo.f -o foo -I$LAMHOME/include */

/* will probably do the trick (assuming that you have set LAMHOME */
/* properly). */

/* That being said, LAM's "mpif77" wrapper compiler should */
/* automatically include the -I option for you.  The following command */
/* should be equivalent to the command listed above: */

/*      unix% mpif77 foo.f -o foo */

/* You should not copy this file to your local directory because it is */
/* possible that this file will be changed between versions of LAM/MPI. */
/* Indeed, this mpif.h is incompatible with the mpif.f of other */
/* implementations of MPI.  Using this mpif.h with other implementations */
/* of MPI, or with other versions of LAM/MPI will result in undefined */
/* behavior (to include incorrect results, segmentation faults, */
/* unexplainable "hanging" in your application, etc.).  Always use the */
/* -I command line option instead (or let mpif77 do it for you). */

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/* LAM version */
/* This file is generated from configure; do not edit it manually. */

/*<        integer LAM_RELEASE_VERSION >*/
/*<        integer LAM_ALPHA_VERSION, LAM_BETA_VERSION >*/
/*<        parameter (LAM_MAJOR_VERSION=6) >*/
/*<        parameter (LAM_MINOR_VERSION=5) >*/
/*<        parameter (LAM_RELEASE_VERSION=9) >*/
/*<        parameter (LAM_ALPHA_VERSION=0) >*/
/*<        parameter (LAM_BETA_VERSION=0) >*/

/* MPI version */

/*<        integer MPI_VERSION, MPI_SUBVERSION >*/
/*<        parameter (MPI_VERSION=1) >*/
/*<        parameter (MPI_SUBVERSION=2) >*/

/* misc. constants */

/*<        integer MPI_SUCCESS, MPI_ANY_SOURCE, MPI_ANY_TAG >*/
/*<        integer MPI_PROC_NULL, MPI_MAX_PROCESSOR_NAME >*/
/*<        integer MPI_MAX_ERROR_STRING, MPI_UNDEFINED >*/
/*<        integer MPI_CART, MPI_GRAPH, MPI_KEYVAL_INVALID >*/
/*<        integer MPI_STATUS_SIZE, MPI_SOURCE, MPI_TAG, MPI_ERROR >*/
/*<        integer MPI_TAG_UB, MPI_HOST, MPI_IO, MPI_WTIME_IS_GLOBAL >*/
/*<        integer MPI_UNIVERSE_SIZE, MPI_APPNUM, MPI_WIN_BASE >*/
/*<        integer MPI_WIN_SIZE, MPI_WIN_DISP_UNIT, MPI_BSEND_OVERHEAD >*/
/*<        integer MPI_MAX_INFO_KEY, MPI_MAX_INFO_VAL >*/
/*<        integer MPI_MAX_PORT_NAME, MPI_MAX_OBJECT_NAME >*/
/*<        integer MPI_ORDER_C, MPI_ORDER_FORTRAN >*/
/*<        integer MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC >*/
/*<        integer MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG >*/
/*<        parameter (MPI_SUCCESS=0) >*/
/*<        parameter (MPI_ANY_SOURCE=-1) >*/
/*<        parameter (MPI_ANY_TAG=-1) >*/
/*<        parameter (MPI_PROC_NULL=-2) >*/
/*<        parameter (MPI_MAX_PROCESSOR_NAME=255) >*/
/*<        parameter (MPI_MAX_ERROR_STRING=255) >*/
/*<        parameter (MPI_UNDEFINED=-32766) >*/
/*<        parameter (MPI_CART=1) >*/
/*<        parameter (MPI_GRAPH=2) >*/
/*<        parameter (MPI_KEYVAL_INVALID=-1) >*/
/*<        parameter (MPI_STATUS_SIZE=4) >*/
/*<        parameter (MPI_SOURCE=1) >*/
/*<        parameter (MPI_TAG=2) >*/
/*<        parameter (MPI_ERROR=3) >*/
/*<        parameter (MPI_TAG_UB=0) >*/
/*<        parameter (MPI_HOST=1) >*/
/*<        parameter (MPI_IO=2) >*/
/*<        parameter (MPI_WTIME_IS_GLOBAL=3) >*/
/*<        parameter (MPI_UNIVERSE_SIZE=4) >*/
/*<        parameter (MPI_APPNUM=5) >*/
/*<        parameter (MPI_WIN_BASE=6) >*/
/*<        parameter (MPI_WIN_SIZE=7) >*/
/*<        parameter (MPI_WIN_DISP_UNIT=8) >*/
/*<        parameter (MPI_BSEND_OVERHEAD=40) >*/
/*<        parameter (MPI_MAX_INFO_KEY=35) >*/
/*<        parameter (MPI_MAX_INFO_VAL=255) >*/
/*<        parameter (MPI_MAX_PORT_NAME=35) >*/
/*<        parameter (MPI_MAX_OBJECT_NAME=63) >*/
/*<        parameter (MPI_ORDER_C=0) >*/
/*<        parameter (MPI_ORDER_FORTRAN=1) >*/
/*<        parameter (MPI_DISTRIBUTE_BLOCK=0) >*/
/*<        parameter (MPI_DISTRIBUTE_CYCLIC=1) >*/
/*<        parameter (MPI_DISTRIBUTE_NONE=2) >*/
/*<        parameter (MPI_DISTRIBUTE_DFLT_DARG=-1) >*/

/* global variables */

/*<        double complex MPI_BOTTOM, MPI_ARGV_NULL >*/
/*<        double complex MPI_ARGVS_NULL, MPI_ERRCODES_IGNORE >*/
/*<        double complex MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE >*/
/*<        common/mpi_bottom/MPI_BOTTOM >*/
/*<        common/mpi_argv_null/MPI_ARGV_NULL >*/
/*<        common/mpi_argvs_null/MPI_ARGVS_NULL >*/
/*<        common/mpi_errcodes_ignore/MPI_ERRCODES_IGNORE >*/
/*<        common/mpi_status_ignore/MPI_STATUS_IGNORE >*/
/*<        common/mpi_statuses_ignore/MPI_STATUSES_IGNORE >*/

/* NULL "handles" (indices) */

/*<        integer MPI_GROUP_NULL, MPI_COMM_NULL, MPI_DATATYPE_NULL >*/
/*<        integer MPI_REQUEST_NULL, MPI_OP_NULL, MPI_ERRHANDLER_NULL >*/
/*<        integer MPI_INFO_NULL >*/
/*<        parameter (MPI_GROUP_NULL=-1) >*/
/*<        parameter (MPI_COMM_NULL=-1) >*/
/*<        parameter (MPI_DATATYPE_NULL=-1) >*/
/*<        parameter (MPI_REQUEST_NULL=-1) >*/
/*<        parameter (MPI_OP_NULL=-1) >*/
/*<        parameter (MPI_ERRHANDLER_NULL=-1) >*/
/*<        parameter (MPI_INFO_NULL=-1) >*/

/* MPI_Init_thread constants */

/*<        integer MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED >*/
/*<        integer MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE >*/
/*<        parameter (MPI_THREAD_SINGLE=0) >*/
/*<        parameter (MPI_THREAD_FUNNELED=1) >*/
/*<        parameter (MPI_THREAD_SERIALIZED=2) >*/
/*<        parameter (MPI_THREAD_MULTIPLE=3) >*/

/* error classes */

/*<        integer MPI_ERR_BUFFER, MPI_ERR_COUNT, MPI_ERR_TYPE >*/
/*<        integer MPI_ERR_TAG, MPI_ERR_COMM, MPI_ERR_RANK >*/
/*<        integer MPI_ERR_REQUEST, MPI_ERR_ROOT, MPI_ERR_GROUP >*/
/*<        integer MPI_ERR_OP, MPI_ERR_TOPOLOGY, MPI_ERR_DIMS >*/
/*<        integer MPI_ERR_ARG, MPI_ERR_UNKNOWN, MPI_ERR_TRUNCATE >*/
/*<        integer MPI_ERR_OTHER, MPI_ERR_INTERN, MPI_ERR_IN_STATUS >*/
/*<        integer MPI_ERR_PENDING, MPI_ERR_SYSRESOURCE >*/
/*<        integer MPI_ERR_LOCALDEAD, MPI_ERR_REMOTEDEAD >*/
/*<        integer MPI_ERR_VALUE, MPI_ERR_FLAGS, MPI_ERR_SERVICE >*/
/*<        integer MPI_ERR_NAME, MPI_ERR_SPAWN, MPI_ERR_KEYVAL >*/
/*<        integer MPI_ERR_INFO_NOKEY, MPI_ERR_WIN >*/
/*<        integer MPI_ERR_EPOCH, MPI_ERR_TYPENOTSUP >*/
/*<        integer MPI_ERR_INFO_KEY, MPI_ERR_INFO_VALUE >*/
/*<        integer MPI_ERR_NO_MEM, MPI_ERR_BASE >*/
/*<        integer MPI_ERR_LASTCODE >*/
/*<        parameter (MPI_ERR_BUFFER=1) >*/
/*<        parameter (MPI_ERR_COUNT=2) >*/
/*<        parameter (MPI_ERR_TYPE=3) >*/
/*<        parameter (MPI_ERR_TAG=4) >*/
/*<        parameter (MPI_ERR_COMM=5) >*/
/*<        parameter (MPI_ERR_RANK=6) >*/
/*<        parameter (MPI_ERR_REQUEST=7) >*/
/*<        parameter (MPI_ERR_ROOT=8) >*/
/*<        parameter (MPI_ERR_GROUP=9) >*/
/*<        parameter (MPI_ERR_OP=10) >*/
/*<        parameter (MPI_ERR_TOPOLOGY=11) >*/
/*<        parameter (MPI_ERR_DIMS=12) >*/
/*<        parameter (MPI_ERR_ARG=13) >*/
/*<        parameter (MPI_ERR_UNKNOWN=14) >*/
/*<        parameter (MPI_ERR_TRUNCATE=15) >*/
/*<        parameter (MPI_ERR_OTHER=16) >*/
/*<        parameter (MPI_ERR_INTERN=17) >*/
/*<        parameter (MPI_ERR_IN_STATUS=18) >*/
/*<        parameter (MPI_ERR_PENDING=19) >*/
/*<        parameter (MPI_ERR_SYSRESOURCE=20) >*/
/*<        parameter (MPI_ERR_LOCALDEAD=21) >*/
/*<        parameter (MPI_ERR_REMOTEDEAD=22) >*/
/*<        parameter (MPI_ERR_VALUE=23) >*/
/*<        parameter (MPI_ERR_FLAGS=24) >*/
/*<        parameter (MPI_ERR_SERVICE=25) >*/
/*<        parameter (MPI_ERR_NAME=26) >*/
/*<        parameter (MPI_ERR_SPAWN=27) >*/
/*<        parameter (MPI_ERR_KEYVAL=28) >*/
/*<        parameter (MPI_ERR_INFO_NOKEY=29) >*/
/*<        parameter (MPI_ERR_WIN=30) >*/
/*<        parameter (MPI_ERR_EPOCH=31) >*/
/*<        parameter (MPI_ERR_TYPENOTSUP=32) >*/
/*<        parameter (MPI_ERR_INFO_KEY=33) >*/
/*<        parameter (MPI_ERR_INFO_VALUE=34) >*/
/*<        parameter (MPI_ERR_NO_MEM=35) >*/
/*<        parameter (MPI_ERR_BASE=36) >*/
/*<        parameter (MPI_ERR_LASTCODE=37) >*/

/* comparison results */

/*<        integer MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL >*/
/*<        parameter (MPI_IDENT=1) >*/
/*<        parameter (MPI_CONGRUENT=2) >*/
/*<        parameter (MPI_SIMILAR=3) >*/
/*<        parameter (MPI_UNEQUAL=4) >*/

/* lookup table indices */

/*<        integer MPI_COMM_WORLD, MPI_COMM_SELF >*/
/*<        integer MPI_GROUP_EMPTY >*/
/*<        integer MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN >*/
/*<        parameter (MPI_COMM_WORLD=0) >*/
/*<        parameter (MPI_COMM_SELF=1) >*/
/*<        parameter (MPI_GROUP_EMPTY=2) >*/
/*<        parameter (MPI_ERRORS_ARE_FATAL=3) >*/
/*<        parameter (MPI_ERRORS_RETURN=4) >*/
/*<        integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION >*/
/*<        integer MPI_COMPLEX, MPI_LOGICAL, MPI_CHARACTER >*/
/*<        integer MPI_BYTE, MPI_PACKED, MPI_UB, MPI_LB, MPI_2REAL >*/
/*<        integer MPI_2DOUBLE_PRECISION, MPI_2INTEGER >*/
/*<        integer MPI_DOUBLE_COMPLEX >*/
/*<        parameter (MPI_BYTE=5) >*/
/*<        parameter (MPI_PACKED=6) >*/
/*<        parameter (MPI_UB=7) >*/
/*<        parameter (MPI_LB=8) >*/
/*<        parameter (MPI_CHARACTER=9) >*/
/*<        parameter (MPI_LOGICAL=10) >*/
/*<        parameter (MPI_INTEGER=11) >*/
/*<        parameter (MPI_REAL=12) >*/
/*<        parameter (MPI_DOUBLE_PRECISION=13) >*/
/*<        parameter (MPI_COMPLEX=14) >*/
/*<        parameter (MPI_DOUBLE_COMPLEX=15) >*/
/*<        parameter (MPI_2REAL=16) >*/
/*<        parameter (MPI_2DOUBLE_PRECISION=17) >*/
/*<        parameter (MPI_2INTEGER=18) >*/
/*<        integer MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND >*/
/*<        integer MPI_BAND, MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR >*/
/*<        integer MPI_MAXLOC, MPI_MINLOC, MPI_REPLACE >*/
/*<        parameter (MPI_MAX=19) >*/
/*<        parameter (MPI_MIN=20) >*/
/*<        parameter (MPI_SUM=21) >*/
/*<        parameter (MPI_PROD=22) >*/
/*<        parameter (MPI_LAND=23) >*/
/*<        parameter (MPI_BAND=24) >*/
/*<        parameter (MPI_LOR=25) >*/
/*<        parameter (MPI_BOR=26) >*/
/*<        parameter (MPI_LXOR=27) >*/
/*<        parameter (MPI_BXOR=28) >*/
/*<        parameter (MPI_MAXLOC=29) >*/
/*<        parameter (MPI_MINLOC=30) >*/
/*<        parameter (MPI_REPLACE=31) >*/

/* attribute functions */

/*<        external MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN >*/
/*<        external MPI_COMM_NULL_COPY_FN, MPI_COMM_NULL_DELETE_FN >*/
/*<        external MPI_TYPE_NULL_COPY_FN, MPI_TYPE_NULL_DELETE_FN >*/
/*<        external MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN >*/
/*<        external MPI_DUP_FN, MPI_COMM_DUP_FN >*/
/*<        external MPI_TYPE_DUP_FN, MPI_WIN_DUP_FN >*/

/* double precision functions */

/*<       double precision MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       external MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK >*/
/*<       integer rowdista(0:*),order(0:*),aptrs(2,0:*),ainds(*) >*/
/*<       integer wrkint(0:*),ranmasks(5,0:*) >*/
/*<       integer N,dd,pp,lgblk,myid,mynnodes,maxnzpercol,comm >*/
/*<       double precision avals(*),pavals(*) >*/
/*      integer, allocatable :: gorder(:),whichsnode(:),tainds(:) */
/*<       integer gorder, whichsnode, tainds >*/
/*<       integer i2 >*/
/*<       dimension gorder(0:N-1) >*/
/*<       dimension whichsnode(0:mynnodes-1) >*/
/*<       dimension tainds(2*i2) >*/
/*     double precision, allocatable :: sendvals(:) */
/*<       double precision sendvals >*/
/*<       integer nsend2 >*/
/*<       dimension sendvals(0:nsend2-1) >*/
/*<       integer proc,pgrsize,ierr,bmaskr,bmaskc,row,col >*/
/*<       integer i,j,k,l,m,ptr_r,fptr_r,ptr_c,itainds >*/
/*<       integer is1,nsend,ptr_sendvals >*/
/*<       integer beginrow >*/
/*<       integer pscv,psdv,prcv,prdv,ppr,ppc,ppg >*/
/*<       pscv = 0 >*/
    /* Parameter adjustments */
    --aptrs;
    --ainds;
    --avals;
    --pavals;
    --ranmasks;
    --tainds;

    /* Function Body */
    pscv = 0;
/*<       psdv = pp >*/
    psdv = *pp;
/*<       prcv = 2*pp >*/
    prcv = *pp << 1;
/*<       prdv = 3*pp >*/
    prdv = *pp * 3;
/*<       ppr  = 8*pp >*/
    ppr = *pp << 3;
/*<       ppc  = 9*pp >*/
    ppc = *pp * 9;
/*<       ppg  = 10*pp >*/
    ppg = *pp * 10;
/*<       pgrsize = ishft(1,ishft(dd,-1)) >*/
    pgrsize = lbit_shift((ftnlen)1, lbit_shift(*dd, (ftnlen)-1));
/*     allocate(gorder(0:N-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'Error in allocate' >*/
	s_wsle(&io___10);
	do_lio(&c__9, &c__1, "Error in allocate", (ftnlen)17);
	e_wsle();
/*<         call mpi_abort(comm,1,ierr) >*/
	mpi_abort__(comm, &c__1, &ierr);
/*<       end if >*/
    }
/*<       beginrow = rowdista(myid) >*/
    beginrow = rowdista[*myid];
/*<       do proc=0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(pscv+proc) = rowdista(proc+1)-rowdista(proc) >*/
	wrkint[pscv + proc] = rowdista[proc + 1] - rowdista[proc];
/*<       end do >*/
    }
/*<    >*/
    mpi_allgatherv__(order, mynnodes, &c__11, gorder, &wrkint[pscv], rowdista,
	     &c__11, comm, &ierr);
/*<       do proc=0,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 0; proc <= i__1; ++proc) {
/*<         wrkint(pscv+proc) = 0 >*/
	wrkint[pscv + proc] = 0;
/*<       end do >*/
    }
/*     allocate(whichsnode(0:mynnodes-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'Error in allocate' >*/
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, "Error in allocate", (ftnlen)17);
	e_wsle();
/*<         call mpi_abort(comm,1,ierr) >*/
	mpi_abort__(comm, &c__1, &ierr);
/*<       end if >*/
    }
/*<       i = aptrs(1,mynnodes-1)+aptrs(2,mynnodes-1)-1 >*/
    i__ = aptrs[(*mynnodes - 1 << 1) + 1] + aptrs[(*mynnodes - 1 << 1) + 2] - 
	    1;
/*     allocate(tainds(2*i),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'Error in allocate' >*/
	s_wsle(&io___16);
	do_lio(&c__9, &c__1, "Error in allocate", (ftnlen)17);
	e_wsle();
/*<         call mpi_abort(comm,1,ierr) >*/
	mpi_abort__(comm, &c__1, &ierr);
/*<       end if >*/
    }
/*<       itainds = 1 >*/
    itainds = 1;
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         col = gorder(beginrow+i) >*/
	col = gorder[beginrow + i__];
/*<         j = 0 >*/
	j = 0;
/*<         do while (col.lt.ranmasks(1,j) .or. col.gt.ranmasks(2,j)) >*/
	while(col < ranmasks[j * 5 + 1] || col > ranmasks[j * 5 + 2]) {
/*<           j = j+1 >*/
	    ++j;
/*<         end do >*/
	}
/*<         proc   = ranmasks(3,j) >*/
	proc = ranmasks[j * 5 + 3];
/*<         bmaskr = ranmasks(4,j) >*/
	bmaskr = ranmasks[j * 5 + 4];
/*<         bmaskc = ranmasks(5,j) >*/
	bmaskc = ranmasks[j * 5 + 5];
/*<         fptr_r = wrkint(ppr+proc) >*/
	fptr_r__ = wrkint[ppr + proc];
/*<         ptr_c  = wrkint(ppc+proc)+iand(ishft(col,-lgblk),bmaskc) >*/
	ptr_c__ = wrkint[ppc + proc] + (lbit_shift(col, -(*lgblk)) & bmaskc);
/*<         k = aptrs(1,i) >*/
	k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	l = aptrs[(i__ << 1) + 2];
/*<         do j=0,l-1 >*/
	i__2 = l - 1;
	for (j = 0; j <= i__2; ++j) {
/*<           tainds(itainds+j) = gorder(ainds(k+j)) >*/
	    tainds[itainds + j] = gorder[ainds[k + j]];
/*<         end do >*/
	}
/*<         call ikeysortf(l,tainds(itainds),tainds(itainds+l)) >*/
	ikeysortf_(&l, &tainds[itainds], &tainds[itainds + l]);
/*<         m = 0 >*/
	m = 0;
/*<         do while(tainds(itainds+m).lt.col) >*/
	while(tainds[itainds + m] < col) {
/*<           m = m+1 >*/
	    ++m;
/*<         end do >*/
	}
/*<         whichsnode(i) = m >*/
	whichsnode[i__] = m;
/*<         do j=m,l-1 >*/
	i__2 = l - 1;
	for (j = m; j <= i__2; ++j) {
/*<           row = tainds(itainds+j) >*/
	    row = tainds[itainds + j];
/*<           ptr_r = fptr_r+iand(ishft(row,-lgblk),bmaskr) >*/
	    ptr_r__ = fptr_r__ + (lbit_shift(row, -(*lgblk)) & bmaskr);
/*<           proc = wrkint(ppg+ptr_c*pgrsize+ptr_r) >*/
	    proc = wrkint[ppg + ptr_c__ * pgrsize + ptr_r__];
/*<           wrkint(pscv+proc) = wrkint(pscv+proc)+1 >*/
	    ++wrkint[pscv + proc];
/*<           tainds(itainds+j) = proc >*/
	    tainds[itainds + j] = proc;
/*<         end do >*/
	}
/*<         itainds = itainds+2*l >*/
	itainds += l << 1;
/*<       end do  >*/
    }
/*     deallocate(gorder) */
/*<       wrkint(psdv) = 0 >*/
    wrkint[psdv] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         wrkint(psdv+proc) = wrkint(psdv+proc-1)+wrkint(pscv+proc-1) >*/
	wrkint[psdv + proc] = wrkint[psdv + proc - 1] + wrkint[pscv + proc - 
		1];
/*<         wrkint(pscv+proc-1) = 0 >*/
	wrkint[pscv + proc - 1] = 0;
/*<       end do >*/
    }
/*<       nsend = wrkint(psdv+pp-1)+wrkint(pscv+pp-1) >*/
    nsend = wrkint[psdv + *pp - 1] + wrkint[pscv + *pp - 1];
/*<       wrkint(pscv+pp-1) = 0 >*/
    wrkint[pscv + *pp - 1] = 0;
/*     allocate(sendvals(0:nsend-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'Error in allocate' >*/
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, "Error in allocate", (ftnlen)17);
	e_wsle();
/*<         call mpi_abort(comm,1,ierr) >*/
	mpi_abort__(comm, &c__1, &ierr);
/*<       end if >*/
    }
/*<       itainds = 1 >*/
    itainds = 1;
/*<       do i=0,mynnodes-1 >*/
    i__1 = *mynnodes - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         k = aptrs(1,i) >*/
	k = aptrs[(i__ << 1) + 1];
/*<         l = aptrs(2,i) >*/
	l = aptrs[(i__ << 1) + 2];
/*<         m = whichsnode(i) >*/
	m = whichsnode[i__];
/*<         do j=m,l-1 >*/
	i__2 = l - 1;
	for (j = m; j <= i__2; ++j) {
/*<           proc = tainds(itainds+j) >*/
	    proc = tainds[itainds + j];
/*<           ptr_sendvals = wrkint(psdv+proc)+wrkint(pscv+proc) >*/
	    ptr_sendvals__ = wrkint[psdv + proc] + wrkint[pscv + proc];
/*<           sendvals(ptr_sendvals) = avals(k+tainds(itainds+l+j)) >*/
	    sendvals[ptr_sendvals__] = avals[k + tainds[itainds + l + j]];
/*<           wrkint(pscv+proc) = wrkint(pscv+proc)+1 >*/
	    ++wrkint[pscv + proc];
/*<         end do >*/
	}
/*<         itainds = itainds+2*l >*/
	itainds += l << 1;
/*<       end do  >*/
    }
/*     deallocate(whichsnode) */
/*     deallocate(tainds) */
/*<    >*/
    mpi_alltoall__(&wrkint[pscv], &c__1, &c__11, &wrkint[prcv], &c__1, &c__11,
	     comm, &ierr);
/*<       wrkint(prdv) = 0 >*/
    wrkint[prdv] = 0;
/*<       do proc=1,pp-1 >*/
    i__1 = *pp - 1;
    for (proc = 1; proc <= i__1; ++proc) {
/*<         wrkint(prdv+proc) = wrkint(prdv+proc-1)+wrkint(prcv+proc-1) >*/
	wrkint[prdv + proc] = wrkint[prdv + proc - 1] + wrkint[prcv + proc - 
		1];
/*<       end do >*/
    }
/*<    >*/
    mpi_alltoallv__(sendvals, &wrkint[pscv], &wrkint[psdv], &c__13, &pavals[1]
	    , &wrkint[prcv], &wrkint[prdv], &c__13, comm, &ierr);
/*     deallocate(sendvals) */
/*<       end >*/
    return 0;
} /* moveav_ */

