
/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 * Modified:
 *     Feburary 7, 2001    use MPI_Isend/MPI_Irecv
 *     October 2, 2001     use MPI_Isend/MPI_Irecv with MPI_Test
 */

#include "superlu_zdefs.h"

#define ISEND_IRECV

/*
 * Function prototypes
 */
#ifdef _CRAY
fortran void CTRSM(_fcd, _fcd, _fcd, _fcd, int*, int*, doublecomplex*,
		   doublecomplex*, int*, doublecomplex*, int*);
fortran void CGEMM(_fcd, _fcd, int*, int*, int*, doublecomplex*, doublecomplex*, 
		   int*, doublecomplex*, int*, doublecomplex*, doublecomplex*, int*);
_fcd ftcs1;
_fcd ftcs2;
_fcd ftcs3;
#endif

/************************************************************************/
void zlsum_fmod
/************************************************************************/
(
 doublecomplex *lsum,    /* Sum of local modifications.                        */
 doublecomplex *x,       /* X array (local)                                    */
 doublecomplex *xk,      /* X[k].                                              */
 doublecomplex *rtemp,   /* Result of full matrix-vector multiply.             */
 int   nrhs,      /* Number of right-hand sides.                        */
 int   knsupc,    /* Size of supernode k.                               */
 int_t k,         /* The k-th component of X.                           */
 int_t *fmod,     /* Modification count for L-solve.                    */
 int_t nlb,       /* Number of L blocks.                                */
 int_t lptr,      /* Starting position in lsub[*].                      */
 int_t luptr,     /* Starting position in lusup[*].                     */
 int_t *xsup,
 gridinfo_t *grid,
 LocalLU_t *Llu,
 MPI_Request send_req[],
 SuperLUStat_t *stat
)
{
/*
 * Purpose
 * =======
 *   Perform local block modifications: lsum[i] -= L_i,k * X[k].
 */
    doublecomplex alpha = {1.0, 0.0}, beta = {0.0, 0.0};
    doublecomplex *lusup, *lusup1;
    doublecomplex *dest;
    int    iam, iknsupc, myrow, nbrow, nsupr, nsupr1, p, pi;
    int_t  i, ii, ik, il, ikcol, irow, j, lb, lk, rel;
    int_t  *lsub, *lsub1, nlb1, lptr1, luptr1;
    int_t  *ilsum = Llu->ilsum; /* Starting position of each supernode in lsum.   */
    int_t  *frecv = Llu->frecv;
    int_t  **fsendx_plist = Llu->fsendx_plist;
    MPI_Status status;
    int test_flag;

    iam = grid->iam;
    myrow = MYROW( iam, grid );
    lk = LBj( k, grid ); /* Local block number, column-wise. */
    lsub = Llu->Lrowind_bc_ptr[lk];
    lusup = Llu->Lnzval_bc_ptr[lk];
    nsupr = lsub[1];

    for (lb = 0; lb < nlb; ++lb) {
	ik = lsub[lptr]; /* Global block number, row-wise. */
	nbrow = lsub[lptr+1];
#ifdef _CRAY
	CGEMM( ftcs2, ftcs2, &nbrow, &nrhs, &knsupc,
	      &alpha, &lusup[luptr], &nsupr, xk,
	      &knsupc, &beta, rtemp, &nbrow );
#elif USE_VENDOR_BLAS
	zgemm_( "N", "N", &nbrow, &nrhs, &knsupc,
	       &alpha, &lusup[luptr], &nsupr, xk,
	       &knsupc, &beta, rtemp, &nbrow, 1, 1 );
#else
	zgemm_( "N", "N", &nbrow, &nrhs, &knsupc,
	       &alpha, &lusup[luptr], &nsupr, xk,
	       &knsupc, &beta, rtemp, &nbrow );
#endif
	stat->ops[SOLVE] += 8 * nbrow * nrhs * knsupc + 2 * nbrow * nrhs;
   
	lk = LBi( ik, grid ); /* Local block number, row-wise. */
	iknsupc = SuperSize( ik );
	il = LSUM_BLK( lk );
	dest = &lsum[il];
	lptr += LB_DESCRIPTOR;
	rel = xsup[ik]; /* Global row index of block ik. */
	for (i = 0; i < nbrow; ++i) {
	    irow = lsub[lptr++] - rel; /* Relative row. */
	    RHS_ITERATE(j)
		z_sub(&dest[irow + j*iknsupc],
		      &dest[irow + j*iknsupc],
		      &rtemp[i + j*nbrow]);
	}
	luptr += nbrow;
		    
	if ( (--fmod[lk])==0 ) { /* Local accumulation done. */
	    ikcol = PCOL( ik, grid );
	    p = PNUM( myrow, ikcol, grid );
	    if ( iam != p ) {
#ifdef ISEND_IRECV
		MPI_Isend( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			   SuperLU_MPI_DOUBLE_COMPLEX, p, LSUM, grid->comm,
                           &send_req[Llu->SolveMsgSent++] );
#else
#ifdef BSEND
		MPI_Bsend( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			   SuperLU_MPI_DOUBLE_COMPLEX, p, LSUM, grid->comm );
#else
		MPI_Send( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			 SuperLU_MPI_DOUBLE_COMPLEX, p, LSUM, grid->comm );
#endif
#endif
#if ( DEBUGlevel>=2 )
		printf("(%2d) Sent LSUM[%2.0f], size %2d, to P %2d\n",
		       iam, lsum[il-LSUM_H], iknsupc*nrhs+LSUM_H, p);
#endif
	    } else { /* Diagonal process: X[i] += lsum[i]. */
		ii = X_BLK( lk );
		RHS_ITERATE(j)
		    for (i = 0; i < iknsupc; ++i)
			z_add(&x[i + ii + j*iknsupc],
			      &x[i + ii + j*iknsupc],
			      &lsum[i + il + j*iknsupc]);
		if ( frecv[lk]==0 ) { /* Becomes a leaf node. */
		    fmod[lk] = -1; /* Do not solve X[k] in the future. */
		    lk = LBj( ik, grid );/* Local block number, column-wise. */
		    lsub1 = Llu->Lrowind_bc_ptr[lk];
		    lusup1 = Llu->Lnzval_bc_ptr[lk];
		    nsupr1 = lsub1[1];
#ifdef _CRAY
		    CTRSM(ftcs1, ftcs1, ftcs2, ftcs3, &iknsupc, &nrhs, &alpha,
			  lusup1, &nsupr1, &x[ii], &iknsupc);
#elif USE_VENDOR_BLAS
		    ztrsm_("L", "L", "N", "U", &iknsupc, &nrhs, &alpha, 
			   lusup1, &nsupr1, &x[ii], &iknsupc, 1, 1, 1, 1);
#else
		    ztrsm_("L", "L", "N", "U", &iknsupc, &nrhs, &alpha, 
			   lusup1, &nsupr1, &x[ii], &iknsupc);
#endif
		    stat->ops[SOLVE] += 4 * iknsupc * (iknsupc - 1) * nrhs
			+ 10 * knsupc * nrhs; /* complex division */
#if ( DEBUGlevel>=2 )
		    printf("(%2d) Solve X[%2d]\n", iam, ik);
#endif
		
		    /*
		     * Send Xk to process column Pc[k].
		     */
		    for (p = 0; p < grid->nprow; ++p) {
			if ( fsendx_plist[lk][p] != EMPTY ) {
			    pi = PNUM( p, ikcol, grid );
#ifdef ISEND_IRECV
			    MPI_Isend( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				       SuperLU_MPI_DOUBLE_COMPLEX, pi, Xk, grid->comm,
				       &send_req[Llu->SolveMsgSent++] );
#else
#ifdef BSEND
			    MPI_Bsend( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				       SuperLU_MPI_DOUBLE_COMPLEX, pi, Xk, grid->comm );
#else
			    MPI_Send( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				     SuperLU_MPI_DOUBLE_COMPLEX, pi, Xk, grid->comm );
#endif
#endif
#if ( DEBUGlevel>=2 )
			    printf("(%2d) Sent X[%2.0f] to P %2d\n",
				   iam, x[ii-XK_H], pi);
#endif
			}
                    }
		    /*
		     * Perform local block modifications.
		     */
		    nlb1 = lsub1[0] - 1;
		    lptr1 = BC_HEADER + LB_DESCRIPTOR + iknsupc;
		    luptr1 = iknsupc; /* Skip diagonal block L(I,I). */

		    zlsum_fmod(lsum, x, &x[ii], rtemp, nrhs, iknsupc, ik,
			       fmod, nlb1, lptr1, luptr1, xsup,
			       grid, Llu, send_req, stat);
		} /* if frecv[lk] == 0 */
	    } /* if iam == p */
	} /* if fmod[lk] == 0 */

    } /* for lb ... */

} /* zLSUM_FMOD */


/************************************************************************/
void zlsum_bmod
/************************************************************************/
(
 doublecomplex *lsum,        /* Sum of local modifications.                    */
 doublecomplex *x,           /* X array (local).                               */
 doublecomplex *xk,          /* X[k].                                          */
 int    nrhs,	      /* Number of right-hand sides.                    */
 int_t  k,            /* The k-th component of X.                       */
 int_t  *bmod,        /* Modification count for L-solve.                */
 int_t  *Urbs,        /* Number of row blocks in each block column of U.*/
 Ucb_indptr_t **Ucb_indptr,/* Vertical linked list pointing to Uindex[].*/
 int_t  **Ucb_valptr, /* Vertical linked list pointing to Unzval[].     */
 int_t  *xsup,
 gridinfo_t *grid,
 LocalLU_t *Llu,
 MPI_Request send_req[],
 SuperLUStat_t *stat
 )
{
/*
 * Purpose
 * =======
 *   Perform local block modifications: lsum[i] -= U_i,k * X[k].
 */
    doublecomplex alpha = {1.0, 0.0};
    int    iam, iknsupc, knsupc, myrow, nsupr, p, pi;
    int_t  fnz, gik, gikcol, i, ii, ik, ikfrow, iklrow, il, irow,
           j, jj, lk, lk1, nub, ub, uptr;
    int_t  *usub;
    doublecomplex *uval, *dest, *y;
    doublecomplex temp;
    int_t  *lsub;
    doublecomplex *lusup;
    int_t  *ilsum = Llu->ilsum; /* Starting position of each supernode in lsum.   */
    int_t  *brecv = Llu->brecv;
    int_t  **bsendx_plist = Llu->bsendx_plist;
    MPI_Status status;
    int test_flag;

    iam = grid->iam;
    myrow = MYROW( iam, grid );
    knsupc = SuperSize( k );
    lk = LBj( k, grid ); /* Local block number, column-wise. */
    nub = Urbs[lk];      /* Number of U blocks in block column lk */

    for (ub = 0; ub < nub; ++ub) {
	ik = Ucb_indptr[lk][ub].lbnum; /* Local block number, row-wise. */
	usub = Llu->Ufstnz_br_ptr[ik];
	uval = Llu->Unzval_br_ptr[ik];
	i = Ucb_indptr[lk][ub].indpos; /* Start of the block in usub[]. */
	i += UB_DESCRIPTOR;
	il = LSUM_BLK( ik );
	gik = ik * grid->nprow + myrow;/* Global block number, row-wise. */
	iknsupc = SuperSize( gik );
	ikfrow = FstBlockC( gik );
	iklrow = FstBlockC( gik+1 );

	RHS_ITERATE(j) {
	    dest = &lsum[il + j*iknsupc];
	    y = &xk[j*knsupc];
	    uptr = Ucb_valptr[lk][ub]; /* Start of the block in uval[]. */
	    for (jj = 0; jj < knsupc; ++jj) {
		fnz = usub[i + jj];
		if ( fnz < iklrow ) { /* Nonzero segment. */
		    /* AXPY */
		    for (irow = fnz; irow < iklrow; ++irow) {
			zz_mult(&temp, &uval[uptr], &y[jj]);
			z_sub(&dest[irow - ikfrow], &dest[irow - ikfrow],
			      &temp);
			++uptr;
		    }
		    stat->ops[SOLVE] += 8 * (iklrow - fnz);
		}
	    } /* for jj ... */
	}

	if ( (--bmod[ik]) == 0 ) { /* Local accumulation done. */
	    gikcol = PCOL( gik, grid );
	    p = PNUM( myrow, gikcol, grid );
	    if ( iam != p ) {
#ifdef ISEND_IRECV
		MPI_Isend( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			   SuperLU_MPI_DOUBLE_COMPLEX, p, LSUM, grid->comm,
                           &send_req[Llu->SolveMsgSent++] );
#else
#ifdef BSEND
		MPI_Bsend( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			   SuperLU_MPI_DOUBLE_COMPLEX, p, LSUM, grid->comm );
#else
		MPI_Send( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			  SuperLU_MPI_DOUBLE_COMPLEX, p, LSUM, grid->comm );
#endif
#endif
#if ( DEBUGlevel>=2 )
		printf("(%2d) Sent LSUM[%2.0f], size %2d, to P %2d\n",
		       iam, lsum[il-LSUM_H], iknsupc*nrhs+LSUM_H, p);
#endif
	    } else { /* Diagonal process: X[i] += lsum[i]. */
		ii = X_BLK( ik );
		dest = &x[ii];
		RHS_ITERATE(j)
		    for (i = 0; i < iknsupc; ++i)
			z_add(&dest[i + j*iknsupc], &dest[i + j*iknsupc],
			      &lsum[i + il + j*iknsupc]);
		if ( !brecv[ik] ) { /* Becomes a leaf node. */
		    bmod[ik] = -1; /* Do not solve X[k] in the future. */
		    lk1 = LBj( gik, grid ); /* Local block number. */
		    lsub = Llu->Lrowind_bc_ptr[lk1];
		    lusup = Llu->Lnzval_bc_ptr[lk1];
		    nsupr = lsub[1];
#ifdef _CRAY
		    CTRSM(ftcs1, ftcs3, ftcs2, ftcs2, &iknsupc, &nrhs, &alpha,
			  lusup, &nsupr, &x[ii], &iknsupc);
#elif USE_VENDOR_BLAS
		    ztrsm_("L", "U", "N", "N", &iknsupc, &nrhs, &alpha, 
			   lusup, &nsupr, &x[ii], &iknsupc, 1, 1, 1, 1);
#else
		    ztrsm_("L", "U", "N", "N", &iknsupc, &nrhs, &alpha, 
			   lusup, &nsupr, &x[ii], &iknsupc);
#endif
		    stat->ops[SOLVE] += 4 * iknsupc * (iknsupc + 1) * nrhs
			+ 10 * iknsupc * nrhs; /* complex division */
#if ( DEBUGlevel>=2 )
		    printf("(%2d) Solve X[%2d]\n", iam, gik);
#endif

		    /*
		     * Send Xk to process column Pc[k].
		     */
		    for (p = 0; p < grid->nprow; ++p) {
			if ( bsendx_plist[lk1][p] != EMPTY ) {
			    pi = PNUM( p, gikcol, grid );
#ifdef ISEND_IRECV
			    MPI_Isend( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				       SuperLU_MPI_DOUBLE_COMPLEX, pi, Xk, grid->comm,
				       &send_req[Llu->SolveMsgSent++] );
#else
#ifdef BSEND
			    MPI_Bsend( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				       SuperLU_MPI_DOUBLE_COMPLEX, pi, Xk, grid->comm );
#else
			    MPI_Send( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				     SuperLU_MPI_DOUBLE_COMPLEX, pi, Xk, grid->comm );
#endif
#endif
#if ( DEBUGlevel>=2 )
			    printf("(%2d) Sent X[%2.0f] to P %2d\n",
				   iam, x[ii-XK_H], pi);
#endif
			}
                     }
		    /*
		     * Perform local block modifications.
		     */
		    if ( Urbs[lk1] )
			zlsum_bmod(lsum, x, &x[ii], nrhs, gik, bmod, Urbs,
				   Ucb_indptr, Ucb_valptr, xsup, grid, Llu,
				   send_req, stat);
		} /* if brecv[ik] == 0 */
	    }
	} /* if bmod[ik] == 0 */

    } /* for ub ... */

} /* zlSUM_BMOD */

