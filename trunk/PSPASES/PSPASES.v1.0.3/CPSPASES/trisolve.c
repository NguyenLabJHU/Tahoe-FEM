/* trisolve.f -- translated by f2c (version 20030320).
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
static integer c__112 = 112;
static integer c__3 = 3;
static integer c__113 = 113;
static integer c__111 = 111;
static integer c__0 = 0;

/* /+***************************************************************************+/ */
/* /+                                                                           +/ */
/* /+   (C) Copyright IBM Corporation, 1997                                     +/ */
/* /+   (C) Copyright Regents of the University of Minnesota, 1997              +/ */
/* /+                                                                           +/ */
/* /+   trisolve.f                                                              +/ */
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
/* /+ $Id: trisolve.c,v 1.1 2004-12-13 08:18:15 paklein Exp $ +/ */
/* /+***************************************************************************+/ */
/*<    >*/
/* Subroutine */ int trisolve_(integer *n, integer *rowdist, integer *order, 
	integer *lptrs, integer *linds, doublereal *lvals, integer *tptrs, 
	integer *tinds, integer *sup, integer *tsind, integer *lc, integer *
	iptrs, integer *ifopts, integer *nrhs, integer *options, doublereal *
	rhso, integer *ldo, doublereal *rhsc, integer *ldc, integer *ranmasks,
	 integer *comm, integer *hvbtemp, integer *lrud, integer *wrkord0, 
	integer *wrkord1, integer *wrkord2, doublereal *ty, doublereal *
	dworkmj, doublereal *ordvals, integer *mynnodes2, integer *bnrhs2, 
	integer *ordvalsiz2, integer *wsolvesize2, integer *trhsize2, integer 
	*wrkord1siz2, integer *wrkord2siz2, integer *ns2, integer *dd2, 
	integer *maxvsize2, integer *hvbsize2)
{
    /* System generated locals */
    integer rhso_dim1, rhso_offset, rhsc_dim1, rhsc_offset, ty_dim1, 
	    ty_offset, i__1, i__2;

    /* Builtin functions */
    integer lbit_shift(integer, integer), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle(void);

    /* Local variables */
    extern /* Subroutine */ int pbsolve1_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *), pfsolve1_(integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *), preordbc_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *), preordbe_(integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *), 
	    getmyhvb_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), preordxc_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, integer *), pbsolvem_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *);
    integer maxhsize, nsupnode, mynnodes;
    extern /* Subroutine */ int pfsolvem_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     doublereal *, integer *, integer *);
    integer maxvsize, i__;
    extern /* Subroutine */ int mpi_abort__(integer *, integer *, integer *);
    integer ordvalsiz, dd, pp, ns, sr, wrkord1siz, wrkord2siz, is1, 
	    supindsize, wsolvesize, rnr, psv, prv;
    extern /* Subroutine */ int getmysnodes_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer pdrc, pdsc, pdsd, pdrd, pirc, pisc, ierr, myid, pisd, pird, ppgr, 
	    nown, psvx, prvx, lgblk, myidh, bnrhs, pdrcx, pdscx, piscx, pisdx,
	     myidv, pircx, pirdx, pdsdx, pdrdx, piown, psloc;
    extern /* Subroutine */ int reordb_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *), reordx_(integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, integer *);
    integer pslocx, rhsptr, uvlptr, hvbsize, uindptr, recvptr, uvecptr, 
	    trhsize, supsize;

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___28 = { 0, 6, 0, 0, 0 };
    static cilist io___30 = { 0, 6, 0, 0, 0 };
    static cilist io___32 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };


/*<       implicit none >*/
/*<       include 'mpif.h' >*/
/*<       double precision zero >*/
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


/*  $Id: trisolve.c,v 1.1 2004-12-13 08:18:15 paklein Exp $ */

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
/*<       parameter(zero=0.d0) >*/
/*<       integer rowdist(0:*),order(0:*),lptrs(3,0:*),linds(*) >*/
/*<       integer tptrs(3,0:*),tinds(*),sup(*),tsind(*),lc(*),iptrs(2,0:*) >*/
/*<       integer ifopts(0:*),options(0:*),nrhs,ldo,ldc,ranmasks(5,0:*) >*/
/*<       double precision lvals(*),rhso(0:ldo-1,*),rhsc(0:ldc-1,*) >*/
/*<       integer N,dd,lgblk,wsolvesize,ns,supindsize,myid,myidh,myidv,pp >*/
/*<       integer comm,nsupnode,supsize,maxvsize,maxhsize,i,j,k,m,hvbsize >*/
/*<       integer uindptr,uvecptr,recvptr,rhsptr,trhsize,is1,nown,ierr >*/
/*<       integer mynnodes,uvlptr,bnrhs,sr,ir,rnr >*/
/*<       integer wrkord1siz,wrkord2siz,ordvalsiz >*/
/*<       integer pisc,pisd,pirc,pird,pdsc,pdsd,pdrc,pdrd >*/
/*<       integer piscx,pisdx,pircx,pirdx,pdscx,pdsdx,pdrcx,pdrdx,ppgr >*/
/*<       integer piown,psloc,pslocx >*/
/*<       integer psv,prv,psvx,prvx >*/
/*     integer, allocatable :: hvbtemp(:),lrud(:) */
/*<       integer hvbtemp, lrud >*/
/*<       integer ns2, dd2, maxvsize2, hvbsize2 >*/
/*<       dimension hvbtemp(hvbsize2 + maxvsize2) >*/
/*<       dimension lrud(ns2+ 4*dd2) >*/
/*     integer, allocatable :: wrkord0(:), wrkord1(:), wrkord2(:) */
/*<       integer wrkord0, wrkord1, wrkord2 >*/
/*<       integer mynnodes2, wrkord1siz2, wrkord2siz2 >*/
/*<       dimension wrkord0(0:mynnodes2-1) >*/
/*<       dimension wrkord1(0:wrkord1siz2-1) >*/
/*<       dimension wrkord2(0:wrkord2siz2-1) >*/
/*     double precision, allocatable :: ty(:,:),dworkmj(:) */
/*<       double precision ty, dworkmj >*/
/*<       integer bnrhs2, wsolvesize2, trhsize2 >*/
/*<       dimension ty(0:N-1,bnrhs2) >*/
/*<       dimension dworkmj(wsolvesize2 + 4*trhsize2) >*/
/*     double precision, allocatable :: ordvals(:) */
/*<       double precision ordvals >*/
/*<       integer ordvalsiz2 >*/
/*<       dimension ordvals(0:ordvalsiz2-1) >*/
/*<       N = ifopts(0) >*/
    /* Parameter adjustments */
    --lptrs;
    --linds;
    --lvals;
    --tptrs;
    --tinds;
    --sup;
    --tsind;
    --lc;
    --iptrs;
    rhso_dim1 = *ldo - 1 - 0 + 1;
    rhso_offset = 0 + rhso_dim1;
    rhso -= rhso_offset;
    rhsc_dim1 = *ldc - 1 - 0 + 1;
    rhsc_offset = 0 + rhsc_dim1;
    rhsc -= rhsc_offset;
    --ranmasks;
    ty_dim1 = *n - 1 - 0 + 1;
    ty_offset = 0 + ty_dim1;
    ty -= ty_offset;
    --dworkmj;
    --lrud;
    --hvbtemp;

    /* Function Body */
    *n = ifopts[0];
/*<       dd = ifopts(1) >*/
    dd = ifopts[1];
/*<       lgblk = ifopts(2) >*/
    lgblk = ifopts[2];
/*<       wsolvesize = ifopts(3)  >*/
    wsolvesize = ifopts[3];
/*<       ns = ifopts(4) >*/
    ns = ifopts[4];
/*<       supindsize = ifopts(5) >*/
    supindsize = ifopts[5];
/*<       myid = ifopts(6) >*/
    myid = ifopts[6];
/*<       myidh = ifopts(7) >*/
    myidh = ifopts[7];
/*<       myidv = ifopts(8) >*/
    myidv = ifopts[8];
/*<       supsize = ifopts(9) >*/
    supsize = ifopts[9];
/*<       bnrhs = min(options(0),nrhs) >*/
    bnrhs = min(options[0],*nrhs);
/*<       pp = ishft(1,dd) >*/
    pp = lbit_shift((ftnlen)1, dd);
/*     allocate(lrud(ns+4*dd),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'Allocate error' >*/
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, "Allocate error", (ftnlen)14);
	e_wsle();
/*<         call mpi_abort(comm,112,ierr) >*/
	mpi_abort__(comm, &c__112, &ierr);
/*<       end if >*/
    }
/*<    >*/
    i__1 = *n - 1;
    getmysnodes_(&i__1, &sup[1], &tinds[1], &tptrs[1], n, &supsize, &lrud[1], 
	    &nsupnode, &lrud[ns + 1], &dd, &maxhsize, &maxvsize, &ns, &myid);
/*<       i = sup(tptrs(3,lrud(nsupnode)))  >*/
    i__ = sup[tptrs[lrud[nsupnode] * 3 + 3]];
/*<       maxvsize = max(maxvsize,lptrs(2,i)) >*/
/* Computing MAX */
    i__1 = maxvsize, i__2 = lptrs[i__ * 3 + 2];
    maxvsize = max(i__1,i__2);
/*<       hvbsize=nsupnode*(10+maxvsize+3*((maxhsize+1)/2+(maxvsize+1)/2)) >*/
    hvbsize = nsupnode * (maxvsize + 10 + ((maxhsize + 1) / 2 + (maxvsize + 1)
	     / 2) * 3);
/*     allocate(hvbtemp(hvbsize+maxvsize),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,myid,': Cannot allocate memory for hvbtemp',hvbsize >*/
	s_wsle(&io___20);
	do_lio(&c__3, &c__1, (char *)&myid, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ": Cannot allocate memory for hvbtemp", (ftnlen)
		36);
	do_lio(&c__3, &c__1, (char *)&hvbsize, (ftnlen)sizeof(integer));
	e_wsle();
/*<         call mpi_abort(comm,113,ierr) >*/
	mpi_abort__(comm, &c__113, &ierr);
/*<       end if >*/
    }
/*<       uindptr = 1+hvbsize >*/
    uindptr = hvbsize + 1;
/*<    >*/
    getmyhvb_(&lrud[1], &nsupnode, &sup[1], &supsize, &tsind[1], &supindsize, 
	    &tptrs[1], &tinds[1], n, &dd, &lgblk, &hvbtemp[1], &hvbsize, &
	    myidh, &myidv, &myid);
/*<       wsolvesize = 4 *wsolvesize * nrhs >*/
    wsolvesize = (wsolvesize << 2) * *nrhs;
/*<       trhsize = maxvsize * nrhs >*/
    trhsize = maxvsize * *nrhs;
/*     allocate(dworkmj(wsolvesize+4*trhsize),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,myid,': Memory Allocation Failure' >*/
	s_wsle(&io___23);
	do_lio(&c__3, &c__1, (char *)&myid, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, ": Memory Allocation Failure", (ftnlen)27);
	e_wsle();
/*<         call mpi_abort(comm,111,ierr) >*/
	mpi_abort__(comm, &c__111, &ierr);
/*<       end if >*/
    }
/*<       uvecptr = 1+wsolvesize >*/
    uvecptr = wsolvesize + 1;
/*<       uvlptr  = uvecptr+trhsize >*/
    uvlptr = uvecptr + trhsize;
/*<       recvptr = uvlptr+trhsize >*/
    recvptr = uvlptr + trhsize;
/*<       rhsptr  = recvptr+trhsize >*/
    rhsptr = recvptr + trhsize;
/*     allocate(ty(0:N-1,bnrhs),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'memory allocation error' >*/
	s_wsle(&io___28);
	do_lio(&c__9, &c__1, "memory allocation error", (ftnlen)23);
	e_wsle();
/*<         call mpi_abort(comm,0,ierr) >*/
	mpi_abort__(comm, &c__0, &ierr);
/*<       end if >*/
    }
/*<       mynnodes = rowdist(myid+1)-rowdist(myid) >*/
    mynnodes = rowdist[myid + 1] - rowdist[myid];
/*     allocate(wrkord0(0:mynnodes-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'memory allocation error' >*/
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, "memory allocation error", (ftnlen)23);
	e_wsle();
/*<         call mpi_abort(comm,0,ierr) >*/
	mpi_abort__(comm, &c__0, &ierr);
/*<       end if >*/
    }
/*<       wrkord1siz = 19*pp >*/
    wrkord1siz = pp * 19;
/*     allocate(wrkord1(0:wrkord1siz-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'memory allocation error' >*/
	s_wsle(&io___32);
	do_lio(&c__9, &c__1, "memory allocation error", (ftnlen)23);
	e_wsle();
/*<         call mpi_abort(comm,0,ierr) >*/
	mpi_abort__(comm, &c__0, &ierr);
/*<       end if >*/
    }
/*<       pisc  = 0 >*/
    pisc = 0;
/*<       pisd  = pp >*/
    pisd = pp;
/*<       pirc  = 2*pp >*/
    pirc = pp << 1;
/*<       pird  = 2*pp+pp >*/
    pird = (pp << 1) + pp;
/*<       pdsc  = 4*pp >*/
    pdsc = pp << 2;
/*<       pdsd  = 4*pp+pp >*/
    pdsd = (pp << 2) + pp;
/*<       pdrc  = 4*pp+2*pp >*/
    pdrc = (pp << 2) + (pp << 1);
/*<       pdrd  = 4*pp+2*pp+pp >*/
    pdrd = (pp << 2) + (pp << 1) + pp;
/*<       piscx = 8*pp >*/
    piscx = pp << 3;
/*<       pisdx = 8*pp+pp >*/
    pisdx = (pp << 3) + pp;
/*<       pircx = 8*pp+2*pp >*/
    pircx = (pp << 3) + (pp << 1);
/*<       pirdx = 8*pp+2*pp+pp >*/
    pirdx = (pp << 3) + (pp << 1) + pp;
/*<       pdscx = 8*pp+4*pp >*/
    pdscx = (pp << 3) + (pp << 2);
/*<       pdsdx = 8*pp+4*pp+pp >*/
    pdsdx = (pp << 3) + (pp << 2) + pp;
/*<       pdrcx = 8*pp+4*pp+2*pp >*/
    pdrcx = (pp << 3) + (pp << 2) + (pp << 1);
/*<       pdrdx = 8*pp+4*pp+2*pp+pp >*/
    pdrdx = (pp << 3) + (pp << 2) + (pp << 1) + pp;
/*<       ppgr  = 16*pp >*/
    ppgr = pp << 4;
/*<    >*/
    preordbe_(n, order, &ranmasks[1], &nown, &wrkord1[pisc], &wrkord1[pisd], &
	    wrkord1[pirc], &wrkord1[pird], &mynnodes, &dd, &myid, &lgblk, &
	    wrkord1[ppgr], wrkord0, comm);
/*<       wrkord2siz = 3*nown+mynnodes >*/
    wrkord2siz = nown * 3 + mynnodes;
/*     allocate(wrkord2(0:wrkord2siz-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'memory allocation error' >*/
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, "memory allocation error", (ftnlen)23);
	e_wsle();
/*<         call mpi_abort(comm,0,ierr) >*/
	mpi_abort__(comm, &c__0, &ierr);
/*<       end if >*/
    }
/*<       piown  = 0 >*/
    piown = 0;
/*<       psloc  = 2*nown >*/
    psloc = nown << 1;
/*<       pslocx = 2*nown+mynnodes >*/
    pslocx = (nown << 1) + mynnodes;
/*<       ordvalsiz = 2*max(mynnodes,nown)*bnrhs >*/
    ordvalsiz = (max(mynnodes,nown) << 1) * bnrhs;
/*     allocate(ordvals(0:ordvalsiz-1),stat=is1) */
/*<       if(is1.ne.0) then >*/
    if (is1 != 0) {
/*<         print *,'memory allocation error' >*/
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, "memory allocation error", (ftnlen)23);
	e_wsle();
/*<         call mpi_abort(comm,0,ierr) >*/
	mpi_abort__(comm, &c__0, &ierr);
/*<       end if >*/
    }
/*<       psv  = 0 >*/
    psv = 0;
/*<       prv  = mynnodes*bnrhs >*/
    prv = mynnodes * bnrhs;
/*<       psvx = 0 >*/
    psvx = 0;
/*<       prvx = nown*bnrhs >*/
    prvx = nown * bnrhs;
/*<    >*/
    preordbc_(n, order, &ranmasks[1], &nown, &wrkord1[pisc], &wrkord1[pisd], &
	    wrkord1[pirc], &wrkord1[pird], &wrkord2[piown], rowdist, &
	    mynnodes, &wrkord2[psloc], &ordvals[psv], &dd, &myid, &lgblk, &
	    wrkord1[ppgr], wrkord0, comm);
/*<    >*/
    preordxc_(n, &wrkord2[piown], rowdist, &mynnodes, &wrkord1[piscx], &
	    wrkord1[pisdx], &wrkord1[pircx], &wrkord1[pirdx], &nown, &wrkord2[
	    pslocx], wrkord0, &ordvals[psvx], &dd, &myid, &ordvals[prvx], 
	    comm);
/*<       do i=0,pp-1 >*/
    i__1 = pp - 1;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<         wrkord1(pdsc+i) = ishft(wrkord1(pisc+i),-1)*bnrhs >*/
	wrkord1[pdsc + i__] = lbit_shift(wrkord1[pisc + i__], (ftnlen)-1) * 
		bnrhs;
/*<         wrkord1(pdsd+i) = ishft(wrkord1(pisd+i),-1)*bnrhs >*/
	wrkord1[pdsd + i__] = lbit_shift(wrkord1[pisd + i__], (ftnlen)-1) * 
		bnrhs;
/*<         wrkord1(pdrc+i) = ishft(wrkord1(pirc+i),-1)*bnrhs >*/
	wrkord1[pdrc + i__] = lbit_shift(wrkord1[pirc + i__], (ftnlen)-1) * 
		bnrhs;
/*<         wrkord1(pdrd+i) = ishft(wrkord1(pird+i),-1)*bnrhs >*/
	wrkord1[pdrd + i__] = lbit_shift(wrkord1[pird + i__], (ftnlen)-1) * 
		bnrhs;
/*<         wrkord1(pdscx+i) = wrkord1(piscx+i)*bnrhs >*/
	wrkord1[pdscx + i__] = wrkord1[piscx + i__] * bnrhs;
/*<         wrkord1(pdsdx+i) = wrkord1(pisdx+i)*bnrhs >*/
	wrkord1[pdsdx + i__] = wrkord1[pisdx + i__] * bnrhs;
/*<         wrkord1(pdrcx+i) = wrkord1(pircx+i)*bnrhs >*/
	wrkord1[pdrcx + i__] = wrkord1[pircx + i__] * bnrhs;
/*<         wrkord1(pdrdx+i) = wrkord1(pirdx+i)*bnrhs >*/
	wrkord1[pdrdx + i__] = wrkord1[pirdx + i__] * bnrhs;
/*<       end do >*/
    }
/*<       do sr=1,nrhs-bnrhs+1,bnrhs   >*/
    i__1 = *nrhs - bnrhs + 1;
    i__2 = bnrhs;
    for (sr = 1; i__2 < 0 ? sr >= i__1 : sr <= i__1; sr += i__2) {
/*<    >*/
	reordb_(n, &rhso[sr * rhso_dim1], ldo, &bnrhs, &ty[ty_offset], &nown, 
		&wrkord1[pdsc], &wrkord1[pdsd], &wrkord1[pdrc], &wrkord1[pdrd]
		, &wrkord2[piown], &mynnodes, &wrkord2[psloc], &ordvals[psv], 
		&ordvals[prv], comm);
/*<       if(bnrhs.ne.1) then >*/
	if (bnrhs != 1) {
/*<    >*/
	    pfsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &bnrhs, &dworkmj[rhsptr], &
		    dworkmj[uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &
		    hvbtemp[uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &
		    hvbsize, &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pfsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &dworkmj[rhsptr], &dworkmj[
		    uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &hvbtemp[
		    uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize,
		     &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<       if(bnrhs.ne.1) then >*/
	if (bnrhs != 1) {
/*<    >*/
	    pbsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &bnrhs, &dworkmj[uvecptr], 
		    &dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pbsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &dworkmj[uvecptr], &
		    dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<    >*/
	reordx_(n, &rhsc[sr * rhsc_dim1], ldc, &bnrhs, &ty[ty_offset], &nown, 
		&wrkord1[pdscx], &wrkord1[pdsdx], &wrkord1[pdrcx], &wrkord1[
		pdrdx], &wrkord2[piown], rowdist, &mynnodes, &wrkord2[pslocx],
		 wrkord0, &ordvals[psvx], &ordvals[prvx], &myid, comm);
/*<       end do  >*/
    }
/*<       rnr = nrhs-sr+1 >*/
    rnr = *nrhs - sr + 1;
/*<       if(rnr.ne.0) then >*/
    if (rnr != 0) {
/*<       do i=0,pp-1 >*/
	i__2 = pp - 1;
	for (i__ = 0; i__ <= i__2; ++i__) {
/*<         wrkord1(pdsc+i) = ishft(wrkord1(pisc+i),-1)*rnr >*/
	    wrkord1[pdsc + i__] = lbit_shift(wrkord1[pisc + i__], (ftnlen)-1) 
		    * rnr;
/*<         wrkord1(pdsd+i) = ishft(wrkord1(pisd+i),-1)*rnr >*/
	    wrkord1[pdsd + i__] = lbit_shift(wrkord1[pisd + i__], (ftnlen)-1) 
		    * rnr;
/*<         wrkord1(pdrc+i) = ishft(wrkord1(pirc+i),-1)*rnr >*/
	    wrkord1[pdrc + i__] = lbit_shift(wrkord1[pirc + i__], (ftnlen)-1) 
		    * rnr;
/*<         wrkord1(pdrd+i) = ishft(wrkord1(pird+i),-1)*rnr >*/
	    wrkord1[pdrd + i__] = lbit_shift(wrkord1[pird + i__], (ftnlen)-1) 
		    * rnr;
/*<         wrkord1(pdscx+i) = wrkord1(piscx+i)*rnr >*/
	    wrkord1[pdscx + i__] = wrkord1[piscx + i__] * rnr;
/*<         wrkord1(pdsdx+i) = wrkord1(pisdx+i)*rnr >*/
	    wrkord1[pdsdx + i__] = wrkord1[pisdx + i__] * rnr;
/*<         wrkord1(pdrcx+i) = wrkord1(pircx+i)*rnr >*/
	    wrkord1[pdrcx + i__] = wrkord1[pircx + i__] * rnr;
/*<         wrkord1(pdrdx+i) = wrkord1(pirdx+i)*rnr >*/
	    wrkord1[pdrdx + i__] = wrkord1[pirdx + i__] * rnr;
/*<       end do >*/
	}
/*<    >*/
	reordb_(n, &rhso[sr * rhso_dim1], ldo, &rnr, &ty[ty_offset], &nown, &
		wrkord1[pdsc], &wrkord1[pdsd], &wrkord1[pdrc], &wrkord1[pdrd],
		 &wrkord2[piown], &mynnodes, &wrkord2[psloc], &ordvals[psv], &
		ordvals[prv], comm);
/*<       if(rnr.ne.1) then >*/
	if (rnr != 1) {
/*<    >*/
	    pfsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &rnr, &dworkmj[rhsptr], &
		    dworkmj[uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &
		    hvbtemp[uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &
		    hvbsize, &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pfsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &dd, &
		    lgblk, n, &ty[ty_offset], &dworkmj[rhsptr], &dworkmj[
		    uvecptr], &dworkmj[uvlptr], &dworkmj[recvptr], &hvbtemp[
		    uindptr], &maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize,
		     &lc[1], &dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<       if(rnr.ne.1) then >*/
	if (rnr != 1) {
/*<    >*/
	    pbsolvem_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &rnr, &dworkmj[uvecptr], &
		    dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       else >*/
	} else {
/*<    >*/
	    pbsolve1_(&lrud[1], &nsupnode, &sup[1], &lptrs[1], &linds[1], &
		    lvals[1], &tptrs[1], &tinds[1], &myid, &myidh, &myidv, &
		    dd, &lgblk, n, &ty[ty_offset], &dworkmj[uvecptr], &
		    dworkmj[recvptr], &dworkmj[rhsptr], &hvbtemp[uindptr], &
		    maxvsize, &lrud[ns + 1], &hvbtemp[1], &hvbsize, &lc[1], &
		    dworkmj[1], &iptrs[1], comm);
/*<       end if >*/
	}
/*<    >*/
	reordx_(n, &rhsc[sr * rhsc_dim1], ldc, &rnr, &ty[ty_offset], &nown, &
		wrkord1[pdscx], &wrkord1[pdsdx], &wrkord1[pdrcx], &wrkord1[
		pdrdx], &wrkord2[piown], rowdist, &mynnodes, &wrkord2[pslocx],
		 wrkord0, &ordvals[psvx], &ordvals[prvx], &myid, comm);
/*<       end if  >*/
    }
/*      deallocate(lrud) */
/*      deallocate(hvbtemp) */
/*      deallocate(ty) */
/*      deallocate(dworkmj) */
/*      deallocate(wrkord0) */
/*      deallocate(wrkord1) */
/*      deallocate(wrkord2) */
/*      deallocate(ordvals) */
/*<       end >*/
    return 0;
} /* trisolve_ */

