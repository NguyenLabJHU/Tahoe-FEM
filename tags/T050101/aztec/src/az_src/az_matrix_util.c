/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_matrix_util.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:14 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_matrix_util.c,v 1.1.1.1 2001-01-30 20:59:14 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include "f2c.h"
#undef abs
#include "az_aztec.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef __MWERKS__
#include <malloc.h>
#endif /* __MWERKS__ - PAK (07/24/98) */
#include <float.h>

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmax_matrix_norm(double val[], int indx[],int bindx[], int rpntr[],
                           int cpntr[], int bpntr[], int proc_config[],
                           int data_org[])

/*******************************************************************************

  This routine returns maximum matrix norm for VBR matrix A: this is a parallel
  version for a distributed matrix.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, maximum matrix norm.
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  register int indx_ptr, irow, iblk_row, jcol, jblk_col, icol_blk, iblk = 0;
  int           num_blk_rows, num_col_blks, num_blk_cols;
  double        row_sum = 0.0, row_max = 0.0;
  int           k;
  int           j_last, bindx_row;

  /**************************** execution begins ******************************/

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

    for (irow = 0; irow < data_org[AZ_N_internal] + data_org[AZ_N_border];
         irow++) {

      /* compute diagonal contribution */

      row_sum = fabs(val[irow]);

      /* nonzero off diagonal contibution */

      j_last  = bindx[irow+1] - bindx[irow];
      bindx_row = bindx[irow];

      for (jcol = 0; jcol < j_last; jcol++) {
        k        = bindx_row + jcol;
        row_sum += fabs(val[k]);
      }
      row_max = max(row_sum, row_max);
    }

    row_max = AZ_gmax_double(row_max, proc_config);

    return row_max;
  }

  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {

    /* loop over the block rows */

    for (iblk_row = 0;
         iblk_row < data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
         iblk_row++) {

      /* find out how many rows are in this block row */

      num_blk_rows = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* find out how many block columns are in this block row */

      num_col_blks = bpntr[iblk_row+1] - bpntr[iblk_row];

      /* loop over all the rows in this block row */

      for (irow = 0; irow < num_blk_rows; irow++) {

        /* loop over all the column blocks in this block row */

        for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

          /* find out which column block this is */

          icol_blk = bindx[iblk];
          indx_ptr = indx[iblk++];

          /* find out how many columns are in this block */

          num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

          /* loop over all the columns in this block */

          for (jcol = 0; jcol < num_blk_cols; jcol++)
            row_sum += fabs(val[indx_ptr + jcol*num_blk_rows + irow]);
        }

        iblk   -= num_col_blks;
        row_max = max(row_sum, row_max);
        row_sum = 0.0;
      }

      iblk += num_col_blks;
    }

    row_max = AZ_gmax_double(row_max, proc_config);

    return row_max;
  }

  else {
    (void) fprintf(stderr, "ERROR: invalid matrix type %d\n",
                   data_org[AZ_matrix_type]);
    exit(1);
  }
  return(0.0);

} /* AZ_gmax_matrix_norm */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gvector_norm(int n, int p, double *x, int proc_config[])

/*******************************************************************************

  Function which returns the lp norm of the vector x, i.e., if p = 2, the
  standard Euclidean norm is returned. NOTE, to get the l-infinity norm,
  set p = -1.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, requested norm value.
  ============

  Parameter list:
  ===============

  n:               Order of vector x.

  p:               Order of the norm to perform, i.e., ||x||p.

  x:               Vector of length n (on this processor).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int    i, j;
  register double sum, power;
  int             index, one = 1, *n_ptr;
  double          norm;

  /**************************** execution begins ******************************/

  /* error checking */

  if (p <= 0 && p != -1) return -1.0;

  /* initialize */

  n_ptr = &n;

  switch (p) {

  case -1:                      /* infinity norm */

    index = idamax_(n_ptr, x, &one) - 1;
    if (index < 0) return -1.0;
    norm  = AZ_gmax_double(fabs(x[index]), proc_config);
    break;

  case 1:                       /* sum norm */
    sum  = dasum_(n_ptr, x, &one);
    norm = AZ_gsum_double(sum, proc_config);
    break;

  case 2:                       /* Euclidean norm */
    norm = sqrt(AZ_gdot(n, x, x, proc_config));
    break;

  default:                      /* any other p-norm */

    sum = 0.0;
    for (i = 0; i < n; i++) {
      power = *x;
      for (j = 0; j < p; j++)
        power *= *x;
      sum += fabs(power);
      x++;
    }

    norm = pow(AZ_gsum_double(sum, proc_config), 1.0 / (double) p);
  }

  return norm;

} /* vector norm */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_vbr_matrix(int matrix_flag, int Proc, int itotal_nodes,
                      int ext_nodes, double val[], int indx[], int bindx[],
                      int rpntr[], int bpntr[])

/*******************************************************************************

  Prints out the VBR matrix and pointers.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  matrix_flag:     = 0 no matrix output, = 1 output matrix.

  Proc:            Current processor number.

  itotal_nodes:    Number of internal + border nodes on this processor.

  ext_nodes:       Number of external nodes.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

*******************************************************************************/

{

  /* local variables */

  int iblk, i, iblk_row, m1, n1, ipoint;
  int ival, jblk, j, ib1, ib2, jpoint;

  /**************************** execution begins ******************************/

  /* print out the VBR indexing information for the matrix */

  (void) printf("\n----- Proc: %d indx -----\n\n", Proc);

  for (iblk = 0; iblk < itotal_nodes; iblk++) {
    for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
      (void) printf("%d ", *(indx+i));

    if (iblk == itotal_nodes - 1){
      (void) printf("%d\n", *(indx+i));
    }
    else
      (void) printf("\n");
  }

  (void) printf("\n----- Proc: %d bindx -----\n\n", Proc);

  for (iblk = 0; iblk < itotal_nodes; iblk++) {
    for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
      (void) printf("%d ", *(bindx+i));
    (void) printf("\n");
  }

  (void) printf("\n----- Proc: %d rpntr -----\n\n", Proc);

  for (i = 0; i < itotal_nodes + ext_nodes + 1; i++)
    (void) printf("%d ", *(rpntr+i));
  (void) printf("\n");

  (void) printf("\n----- Proc: %d bpntr -----\n\n", Proc);

  for (i = 0; i < itotal_nodes + 1; i++)
    (void) printf("%d ", *(bpntr+i));
  (void) printf("\n");

  /* dump of matrix in a block output format */

  if (matrix_flag) {

    /* loop over block rows */

    for (iblk_row = 0; iblk_row < itotal_nodes; iblk_row++) {

      /* number of rows in the current row block */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row]];

      /* loop over all the blocks in the current block-row */

      for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = rpntr[jblk];

        /* ending point column index of the current block */

        ib2 = rpntr[jblk+1];

        /* number of columns in the current block */

        n1 = ib2 - ib1;

        (void) printf("\nProc: %d Block Row: %d Block Column: %d "
                      "Row Pointer: %d Column Pointer: %d\n", Proc, iblk_row,
                      jblk, rpntr[iblk_row], rpntr[jblk]);

        (void) printf("---------------------------------------"
                      "---------------------------------------\n");

        for (ipoint = 0; ipoint < m1; ipoint++) {
          for (jpoint = 0; jpoint < n1; jpoint++)
            (void) printf("val[%d]: %e ", ival+jpoint*m1+ipoint,
                          val[ival+jpoint*m1+ipoint]);
          (void) printf("\n");
        }

        ival += m1*n1;
      }
    }
  }

} /* print_vbr_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_dtrans(int *m, int *n, double *A)

/*******************************************************************************

  Perform an in-place transpose of a general m x n matrix stored in "A".

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m, n:            Number of rows, columns in "A", respectively.

  A:               Matrix to be transposed.

*******************************************************************************/

{

  /* local variables */

  register int  i, j, w_index = 0, A_index;
  int           itemp;
  double       *work;

  /**************************** execution begins ******************************/

  /* calloc temporary array */

  work = (double *) calloc(*m * *n, sizeof(double));

  /* do the transpose */

  for (i = 0; i < *n; i++)
    for (j = 0; j < *m; j++) {
      A_index = i + j * *n;
      *(work + w_index++) = *(A + A_index);
    }

  /* copy from "work" back to "A" */

  for (i = 0; i < *m * *n; i++)
    *(A + i) = *(work + i);

  /* free "work" */

  free((void *) work);

  /* exchange m and n */

  itemp = *m;
  *m    = *n;
  *n    = itemp;

} /* dtrans */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_get_sym_indx(int iblk, int jblk, int *indx, int *bindx, int *bpntr)

/*******************************************************************************

  Given a block-row and block-column index, return the index to the symmetric
  starting point of the matrix.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     int, starting index into val of the symmetric portion of the
  ============          matrix.

  Parameter list:
  ===============

  iblk:            Current block row index.

  jblk:            Current block column index.

  indx,
  bindx,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

*******************************************************************************/

{

  register int i, icount = 0;
  int          itemp, pre_num_nz_blks, num_nz_blks;
  int         *bindx_ptr;

  itemp     = bpntr[jblk];
  bindx_ptr = bindx + itemp;

  /*
   * Count up how many nonzero block precede the current iblk' in the current
   * block row.
   */

  num_nz_blks = bpntr[jblk+1] - itemp;

  for (i = 0; i < num_nz_blks; i++) {
    if (*(bindx_ptr + icount) == iblk) {
      pre_num_nz_blks = icount;
      break;
    }
    else
      icount++;
  }

  return indx[itemp + pre_num_nz_blks];

} /* AZ_get_sym_indx */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

extern int AZ_sys_msg_type;

void AZ_print_out(int update_index[], int extern_index[], int update[], 
	int external[], double val[], int indx[], int bindx[], int rpntr[], 
	int cpntr[], int bpntr[], int proc_config[], int choice, int matrix, 
	int N_update, int N_external, int of)
{
/*******************************************************************************

  Print out the matrix in 1 of 3 formats described below.
  starting point of the matrix.

  Author:          Ray Tuminaro, SNL, 9222
  =======

  Return code:     none.
  ============ 

  Parameter list:
  ===============

  update_index,    On input, ordering of update and external locally on this
  extern_index:    processor. For example  'update_index[i]' gives the index
                   location of the block which has the global index 'update[i]'.
                   (AZ_global_mat only).
 
  update:          On input, blks updated on this node (AZ_global_mat only).

  external:        On input, list of external blocks (AZ_global_mat only).

  val,bindx        On input, matrix (MSR or VBR) arrays holding matrix values.
  indx, bnptr,     When using either AZ_input_form or AZ_explicit, these can
  rnptr, cnptr:    be either pre or post-AZ_transform() values depending on what
                   the user wants to see. When using AZ_global_form, these must
                   be post-AZ_transform() values.

  proc_config:     On input, processor information corresponding to:
                      proc_config[AZ_node] = name of this processor
                      proc_config[AZ_N_procs] = # of processors used
 
  choice:          On input, 'choice' determines the output to be printed
		   as described above.
 
  matrix:          On input, type of matrix (AZ_MSR_MATRIX or AZ_VBR_MATRIX).

  N_update:        On input, number of points/blks to be updated on this node.

  N_external:      On input, number of external points/blks needed by this node.

  of:              On input, an offset used with AZ_global_matrix and 
		   AZ_explicit. In particular, a(of,of) is printed for 
                   the matrix element stored as a(0,0).

*******************************************************************************/

   int type, neighbor, cflag;
   int ii,i,j,tj;
   int iblk, jblk, m1, n1, ival, new_iblk, new_jblk;
   MPI_Request request;  /* Message handle */

 
   type            = AZ_sys_msg_type;
   AZ_sys_msg_type =(AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +AZ_MSG_TYPE;

   /* Synchronize things so that processor 0 prints first and then */
   /* sends a message to processor 1 so that he prints, etc.      */

   neighbor = proc_config[AZ_node] - 1;
   if ( proc_config[AZ_node] != 0) {
      md_wrap_iread((void *) &i, 0, &neighbor, &type, &request);
      md_wrap_wait((void *) &i, 0, &neighbor, &type, &cflag, &request);
   }
   printf("proc %d:\n",proc_config[AZ_node]);

   if (choice == AZ_input_form ) {
     if ( update != (int *) NULL) {
        printf("  N_update: %5d\n",N_update); printf("  update:");
        AZ_list_print(update, N_update, (double *) NULL , 0);
     }

     if (matrix == AZ_MSR_MATRIX) {
        printf("  bindx: ");
        AZ_list_print(bindx, bindx[N_update], (double *) NULL , 0);

        printf("  val:   ");
        AZ_list_print((int *) NULL , N_update, val , bindx[N_update]);
     }
     else if (matrix == AZ_VBR_MATRIX) {
        printf("  rpntr: ");
        AZ_list_print(rpntr, N_update+1, (double *) NULL , 0);
        if ( cpntr != (int *) NULL ) {
           printf("  cpntr: ");
           AZ_list_print(cpntr, N_update+1+ N_external, (double *) NULL , 0);
        }
        printf("  bpntr: ");
        AZ_list_print(bpntr, N_update+1, (double *) NULL , 0);
        printf("  bindx: ");
        AZ_list_print(bindx, bpntr[N_update], (double *) NULL , 0);
        printf("  indx:  ");
        AZ_list_print(indx, bpntr[N_update]+1, (double *) NULL , 0);
        printf("  val:   ");
        AZ_list_print((int *) NULL, indx[bpntr[N_update]], val, 0);
     }
   }
   else if (choice == AZ_global_mat ) {
     if ( matrix == AZ_MSR_MATRIX) {
        for (i = 0 ; i < N_update; i++ ) {
          ii = update_index[i];
          printf("   a(%d,%d) = %20.13e;\n",update[i]+of,update[i]+of,val[ii]);
          for (j = bindx[ii] ; j < bindx[ii+1] ; j++ ) {
            tj = AZ_find_simple(bindx[j], update_index, N_update, extern_index,
                              N_external,update,external);
            if (tj != -1) 
               printf("   a(%d,%d) = %20.13e;\n",update[i]+of,tj+of,val[j]);
            else (void) fprintf(stderr,"col %d (= bindx[%d]) is undefined\n",
                                tj, j);
          }
        }
     }
     else if (matrix == AZ_VBR_MATRIX) {
        for (iblk= 0; iblk < N_update; iblk++) {
           new_iblk = update_index[iblk];

           m1 = rpntr[new_iblk+1] - rpntr[new_iblk];
 
           /* loop over blocks in the current block-row */
 
           for (ii = bpntr[new_iblk]; ii < bpntr[new_iblk+1]; ii++) {
              new_jblk = AZ_find_simple(bindx[ii], update_index, N_update, 
			 extern_index, N_external,update,external);
              if (new_jblk == -1) {
                 printf("local column %d not found\n",new_jblk);
                 exit(-1);
              }
              jblk = bindx[ii];
              ival =  indx[ii];

              n1 = cpntr[jblk+1] - cpntr[jblk];

              for (i = 0; i < m1; i++) {
                 for (j = 0; j < n1; j++)
                    (void) printf("   a(%d(%d),%d(%d)) = %20.13e;\n",update[iblk]+
				  of,i+of, new_jblk+of, j+of, val[ival+j*m1+i]);
              }
           }
        }
     }
   }
   else if (choice == AZ_explicit) {
     if ( matrix == AZ_MSR_MATRIX) {
        for (i = 0 ; i < N_update; i++ ) {
          if (update == NULL) tj = i+of;
          else tj = update[i] + of;
          printf("   a(%d,%d) = %20.13e;\n",tj,tj,val[i]);
          for (j = bindx[i] ; j < bindx[i+1] ; j++ ) {
               printf("   a(%d,%d) = %20.13e;\n",tj,bindx[j]+of,val[j]);
          }
        }
     }
     else if (matrix == AZ_VBR_MATRIX) {
        for (iblk = 0; iblk < N_update; iblk++) {
           if (update == NULL) tj = iblk+of;
           else tj = update[iblk] + of;

           m1 = rpntr[iblk+1] - rpntr[iblk];
 
           /* loop over blocks in the current block-row */
 
           for (ii = bpntr[iblk]; ii < bpntr[iblk+1]; ii++) {
              jblk = bindx[ii];
              ival =  indx[ii];
              n1 = (indx[ii+1]-ival)/m1;

              for (i = 0; i < m1; i++) {
                 for (j = 0; j < n1; j++)
                    (void) printf("   a(%d(%d),%d(%d)) = %20.13e;\n", tj, 
				  i+of, jblk+of, j+of, val[ival+j*m1+i]);
              }
 
           }
        }
     }
   }
   else (void) fprintf(stderr,"AZ_matrix_out: output choice unknown\n");

   neighbor = proc_config[AZ_node] + 1;
   if ( proc_config[AZ_node] != proc_config[AZ_N_procs] - 1) 
      md_wrap_write((char *) &i, 0, neighbor, type, &cflag);

   i = AZ_gsum_int(i,proc_config);


}
int AZ_find_simple(int k, int *update_index, int N_update, int *extern_index,
   int N_external, int *update, int *external)
{
   int i;

   if (k < N_update) {
      for (i = 0 ; i < N_update ; i++ ) 
         if ( update_index[i] == k) return(update[i]);
   }
   else {
      for (i = 0 ; i < N_external ; i++ ) 
         if ( extern_index[i] == k) return(external[i]);
   }

   return(-1);
}
void AZ_list_print(int ivec[] , int length1, double dvec[],int length2)

{
   int i, count = 0;

   if (ivec != (int *) NULL) {
      for (i = 0 ; i < length1 ; i++ ) {
         printf("%8d ",ivec[i]); count += 8;
         if (count > 50) { count = 0; printf("\n         "); }
      }
   }
   else if (dvec != (double *) NULL) {
      for (i = 0 ; i < length1 ; i++ ) {
         printf("%8.1e ",dvec[i]); count += 8;
         if (count > 50) { count = 0; printf("\n         "); }
      }
      if (length2 != 0) { 
         printf("      -- "); count += 8;
         if (count > 50) { count = 0; printf("\n         "); }
      }

      for (i = length1+1 ; i < length2 ; i++ ) {
         printf("%8.1e ",dvec[i]); count += 8;
         if (count > 50) { count = 0; printf("\n         "); }
      }
   }
   printf("\n");
}
