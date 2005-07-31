/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_domain_setup.c,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:13 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id: az_domain_setup.c,v 1.1.1.1 2001-01-30 20:59:13 paklein Exp $";
#endif


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


extern int AZ_getit(unsigned newbuf[], unsigned location);
extern void AZ_setit(unsigned newbuf[], unsigned location, unsigned value);
extern void AZ_grabit(int string[],int count,double *out,int bits);
extern void AZ_compress_it(double full[],int last_f, int compr[], int *last_c,
                           int row_bits, int proc_bits, int common_proc);
extern void AZ_uncompress_it(double full[], int *last_f, int compr[],
                             int last_c, int row_bits, int proc_bits,
                             int common_proc);
extern void AZ_set_integer(int exponent,int string[],int count,int exp_bits);

#include "f2c.h"
#undef abs
#include "az_aztec.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void AZ_msr2ilu_setup(double val[], int bindx[], double tval[], int tbindx[],
                      int new_N, int NN, int options[], int read_index,
                      double buffer[], int data_org[])

/*******************************************************************************

  Set up the submatrix to be factored by ilu as a domain decomposition
  preconditioner. Essentially, we copy the old submatrix and for the external
  points, we either setup Dirichlet conditions with 1's on the diagonal, setup
  Dirichlet conditions using the diagonal from the corresponding neighboring
  processor, or we use the rows from the neighboring processors responsible for
  updating the corresponding external point.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix on the
                   current processor (see file params.txt).

  bindx:           Pointer array used for DMSR and DVBR sparse matrix storage
                   (see file params.txt).

  tval,
  tbindx:          On output, new submatrix contains rows corresponding to the
                   external points.

  new_N:           Order of the enlargened system of equations on the current
                   processor.

  NN:              Order of the original system of equations defined on the
                   current processor.

  options:         Determines specific solution method and other parameters.

  read_index:      Length of storage buffer.

  buffer:          Storage buffer for exchanged rows (i.e. external rows).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

*******************************************************************************/

{

  /* local variables */

  int i, diff, tind;
  int temp;

  /**************************** execution begins ******************************/

  if (NN == 0) {
    tbindx[NN] = 0;
    return;
  }

  /* first copy old entries */

  diff = new_N - NN;
  for (i = bindx[NN]-1; i >= bindx[0]; i--) {
    tval[i+diff]   = val[i];
    tbindx[i+diff] = bindx[i];
  }

  for (i = 0; i < NN; i++) {
    tval[i]   = val[i];
    tbindx[i] = bindx[i] + diff;
  }

  tbindx[NN] = bindx[NN] + diff;
  temp       = tbindx[NN];

  if (options[AZ_overlap] == AZ_full) {

    /*
     * Add in the new nonzeros corresponding to the external rows. We do this by
     * first calculating the number of nonzeros in each row and recording this
     * in tbindx[]. We then use this number to compute the pointers to the
     * nonzeros and record these in tbindx[].
     */

    for (i = NN; i < new_N; i++) tbindx[i] = 0;

    i = 0;
    while (i < read_index) {
      while ( (i < read_index) && (buffer[i] == -1.0) ) i++;

      if (i < read_index) {
        tind = (int) buffer[i];
        i += 3;
        while (buffer[i] != -1.0) {
          i += 3;
          tbindx[tind]++;
        }
      }
    }

    for (i = NN; i < new_N; i++) {
      tind      = tbindx[i];
      tbindx[i] = temp;
      temp     += tind;
    }

    tbindx[new_N] = temp;

    /* Store the new nonzero entries in the matrix */

    i = 0;
    while (i < read_index) {
      while ( (i < read_index) && (buffer[i] == -1.0) ) i++;

      if (i < read_index) {
        tind       = (int) buffer[i];
        tval[tind] = buffer[i+1];
        tind       = tbindx[tind];
        i         += 3;

        while (buffer[i] != -1.0) {
          tbindx[tind] = (int) buffer[i];
          tval[tind++] = buffer[i+1];
          i           += 3;
        }
      }
    }

    if (buffer != 0) free((char *) buffer);
  }

  else {

    /* add Dirichlet conditions for the external points */

    for (i = NN; i < new_N; i++) {
      tval[i]     = 1.0;
      tbindx[i+1] = tbindx[i];
    }

    /* communicate the diagonals */

    if (options[AZ_overlap] == AZ_diag) {
      AZ_exchange_bdry(tval, data_org);
    }
  }

} /* AZ_msr2ilu_setup */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_extract(int row, double val[], int bindx[], int node, double buffer[],
                int *next, double mapper[], double mapper2[], int NN)

/*******************************************************************************

  Pull out the row denoted by 'row' and store it into 'buffer'. Rows are stored
  in the following fashion:
  a) each row is proceeded and terminated by the number -1.
  b) the first element of the row sent is the diagonal.
  c) each element is represented by a triplet of numbers (in double precision).
     The first number is the column index, the second number is the value of the
     matrix element, and the third number is the processor responsible for
     updating the value associated with the column.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  row:             Row to be stored into buffer.

  val:             Array containing the nonzero entries of the matrix on the
                   current processor (see file params.txt).

  bindx:           Pointer array used for DMSR and DVBR sparse matrix storage
                   (see file params.txt).

  node:            Current processor number.

  buffer:          Storage buffer for extracted rows.

  next:            On ouput, next available location in 'buffer' to store a new
                   element.

  mapper:          mapper[i] corresponds to the local index of the variable 'i'
                   on the processor responsible for updating this variable.

  mapper2:         mapper2[i] corresponds to the processor responsible for
                   updating the variable 'i'.

  NN:              Order of the original system of equations defined on the
                   current processor.

*******************************************************************************/

{

  /* local variables */

  int count = 0, i, start, end;

  /**************************** execution begins ******************************/

  buffer[count++] = (double) -1;
  buffer[count++] = (double) row;
  buffer[count++] = val[row];
  buffer[count++] = (double) node;

  start = bindx[row];
  end   = bindx[row+1];

  for (i = start; i < end; i++) {
    if (bindx[i] >= NN) {

      /* element corresponds to an external point */

      buffer[count++] = (double) mapper[bindx[i]];
      buffer[count++] = val[i];
      buffer[count++] = (double) mapper2[bindx[i]];
    }
    else {
      buffer[count++] = (double) bindx[i];
      buffer[count++] = val[i];
      buffer[count++] = (double) node;
    }
  }

  buffer[count] = (double) -1;
  *next         = count + 1;

} /* extract */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_exchange_rows(double val[], int indx[], int *nonzeros, int new_N,
                      double **tbuf, int *buf_length, int bindx[], 
                      int cpntr[], int bpntr[], int data_org[],
                      int proc_config[])

/*******************************************************************************

  Routine to exchange internal boundary rows between neighboring processors.
  This exchange is done by
  1) determining which processors are responsible for updating external points
     and what are the correcponding local indices for those points on that
     processor.
  2) determining the message length necessary to communicate the nonzeros in the
     rows corresponding to list_send_unknowns[].
  3) send the rows to the neighboring processors. Note: rows are sent in the
     following fashion:
     a) each row is proceeded and terminated by the number -1.
     b) the first element of the row sent is the diagonal.
     c) each element is represented by a triplet of numbers (in double
        precision). The first number is the column index, the second number is
        the value of the matrix element, and the third number is the processor
        responsible for updating the value associated with the column.
  4) receive the rows.
  5) replace column indecies that are defined on other processors by column
     indices on this processor (if possible).
  6) throw out any column indices that are not defined on this processor.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  nonzeros:        On output is equal to the number of nonzero matrix elements
                   received.

  new_N:           Order of the enlargened system of equations on the current
                   processor.

  tbuf:            Storage buffer for exchanged rows.

  buf_length:      Size of storage buffer for exchanged rows.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int     num_nonzeros_recv[AZ_MAX_NEIGHBORS];
  int     actual_send_length[AZ_MAX_NEIGHBORS];
  int     start_send_proc[AZ_MAX_NEIGHBORS];
  int     i, j, n;
  int     nz_count;
  int     mesg_from,rtype, size;
  int     total_send,total_recv;
  int     next, length;
  int     st, type;
  int     kk, k2;
  int     block_row;
  int    *node_name;
  int    *node_ptr;
  int     num_externs, otnode, tnode, *extern_loc_index, *extern_my_index, tind;
  int     newstart, *counts;
  int     Proc;
  int     num_neighbors, *proc_num_neighbor, *num_unknowns_send_neighbor;
  int    *list_send_unknowns, NN, blks;
  double *responsible_proc, *local_index;
  double *buffer;
  int    total_num_recv = 0;
  int    proc_bits, row_bits;
  MPI_Request request[AZ_MAX_NEIGHBORS];  /* Message handle */


  /* external variables */

  extern int AZ_sys_msg_type;

  extern void  AZ_splitup_big_msg(int num_neighbors, double  *buffer,
                                  int *start_send_proc,
                                  int *actual_send_length,
                                  int *actual_recv_length,
                                  int *proc_num_neighbor,
                                  int type, int *total_num_recv,
                                  int *proc_config);

  /**************************** execution begins ****************************/

  num_neighbors              = data_org[AZ_N_neigh];
  proc_num_neighbor          = &data_org[AZ_neighbors];
  num_unknowns_send_neighbor = &data_org[AZ_send_length];
  list_send_unknowns         = &data_org[AZ_send_list];
  NN                         = data_org[AZ_N_internal]+ data_org[AZ_N_border];
  blks                       = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
  Proc                       = proc_config[AZ_node];

  *tbuf       = 0;
  *buf_length = 0;


  /*
   * Setup 2 arrays: responsible_proc[] and local_index[]. These two arrays will
   * be used to determine which processors update which external point and what
   * is their local index on that processor.
   */

  responsible_proc = (double *) calloc(new_N+1, sizeof(double));
  local_index      = (double *) calloc(new_N+1, sizeof(double));
  if (local_index == NULL) {
    (void) fprintf(stderr, "no space in AZ_exchange_rows: local_index\n");
    exit(-1);
  }

  for (i = 0; i < new_N; i++) {
    responsible_proc[i] = (double) Proc;
    local_index[i]      = (double) i;
  }

  AZ_exchange_bdry(local_index, data_org);
  AZ_exchange_bdry(responsible_proc, data_org);

  /*
   * Compute the message length to be sent to each processor in order to sending
   * the rows corresponding to the array list_send_unknowns[].  Send the message
   * length to each processor.
   */

  block_row = 0;

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +
    AZ_MSG_TYPE;

  /* post receives */

  counts = (int *) calloc(num_neighbors+1,sizeof(int));
  if (counts == NULL) {
    (void) fprintf(stderr, "no space in AZ_exchange_rows: counts\n");
    exit(-1);
  }

  for (n = 0; n < num_neighbors; n++) {
    mesg_from = proc_num_neighbor[n];
    rtype     = type;
    (void) md_wrap_iread((void *) &(counts[n]), sizeof(int), &mesg_from,
                 &rtype, request+n);
  }


  /* send lengths */

  j = total_send = total_recv = 0;

  for (n = 0; n < num_neighbors; n++) {
    nz_count = 0;

    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      for (i = 0; i < num_unknowns_send_neighbor[n]; i++) {
        nz_count += 2;

        /* figure out the block row */

        AZ_which_block_row(&block_row, cpntr, list_send_unknowns[j]);

        for (kk = bpntr[block_row]; kk < bpntr[block_row+1]; kk++) {
          nz_count += 3 * (cpntr[bindx[kk]+1] - cpntr[bindx[kk]]);
        }

        j++;
      }
    }
    else {
      for (i = 0; i < num_unknowns_send_neighbor[n]; i++) {
        nz_count += 5;
        nz_count += 3 * (bindx[list_send_unknowns[j]+1] -
                         bindx[list_send_unknowns[j]]);

        j++;
      }
    }

    total_send += nz_count;
    (void) md_wrap_write((void *) &nz_count, sizeof(int), proc_num_neighbor[n],
                  type, &st);
  }

  /*
   * Receive the message length from each processor necessary to fill in the
   * rows corresponding to the external variables.
   */


  for (n = 0; n < num_neighbors; n++) {
    mesg_from = proc_num_neighbor[n];
    rtype     = type;
    (void) md_wrap_wait((void *) &(counts[n]), sizeof(int), &mesg_from,
                 &rtype, &st, request+n);
    total_recv += counts[n];
  }
  free(counts);

  /* allocate message buffer */

  size   = max(total_send, total_recv) + 10;
  buffer = (double *) calloc(size+1, sizeof(double));
  if ( buffer == NULL) {
    (void) fprintf(stderr, "no space in AZ_exchange_rows: buffer\n");
    exit(-1);
  }
  *tbuf = buffer;

  /* extract the rows (in an intermediate form) */

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +
    AZ_MSG_TYPE;

  /* figure out the number of bits needed to encode all the processors */

  proc_bits = 0;
  j         = proc_config[1]-1;
  while(j != 0) {
    proc_bits++;
    j /= 2;
  }

  /* figure out the number of bits needed to encode all the rows + 1 */

  row_bits = 0;
  j = AZ_gmax_int(data_org[AZ_N_internal] + data_org[AZ_N_border] +
                  data_org[AZ_N_external], proc_config);
  j++;
  while(j != 0) {
    row_bits++;
    j /= 2;
  }

  j = start_send_proc[0] = 0;

  for (n = 0; n < num_neighbors; n++) {
    next = actual_send_length[n] = 0;

    for (i = 0; i < num_unknowns_send_neighbor[n]; i++) {

      if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
        (void) AZ_newextract(list_send_unknowns[j], val, indx, Proc,
                             &(buffer[start_send_proc[n]+next]), &length,
                             local_index, responsible_proc, NN, blks, bindx,
                             cpntr, bpntr);

      }
      else
        (void) AZ_extract(list_send_unknowns[j], val, bindx, Proc,
                          &(buffer[start_send_proc[n]+next]), &length,
                          local_index, responsible_proc, NN);

      next += length;
      j++;
    }
    actual_send_length[n] = next;
    start_send_proc[n+1]  = start_send_proc[n] + actual_send_length[n];
  }

  /*
   * Send a message to neighbors indicating how many nonzeros will actually be
   * sent over.
   */

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +
    AZ_MSG_TYPE;

  for (n = 0; n < num_neighbors; n++) {
    mesg_from = proc_num_neighbor[n];
    rtype     = type;
    md_wrap_iread((void *) &num_nonzeros_recv[n], sizeof(int), &mesg_from,
                 &rtype, request+n);
  }

  for (n = 0; n < num_neighbors; n++) {
    md_wrap_write((void *) &actual_send_length[n], sizeof(int),
                  proc_num_neighbor[n], type, &st);
  }

  for (n = 0; n < num_neighbors; n++) {
    mesg_from = proc_num_neighbor[n];
    rtype     = type;
    (void) md_wrap_wait((void *) &num_nonzeros_recv[n], sizeof(int), &mesg_from,
                 &rtype, &st, request+n);
  }

  /* Splitup large messages so that communication buffer cannot overflow */

  AZ_splitup_big_msg(num_neighbors, buffer, start_send_proc,
                     actual_send_length, num_nonzeros_recv,
                     proc_num_neighbor, type, &total_num_recv,
                     proc_config);

  /* make a list of the external points */

  node_name = (int *) calloc(num_neighbors+1, sizeof(int));
  node_ptr  = (int *) calloc(num_neighbors+1, sizeof(int));
  if ( node_ptr== NULL) {
    (void) fprintf(stderr, "no space in AZ_exchange_rows: node_ptr\n");
    exit(-1);
  }

  for (i = 0; i < num_neighbors; i++) node_name[i] = proc_num_neighbor[i];

  AZ_sort(node_name, num_neighbors, NULL, NULL);

  num_externs = new_N - NN;
  otnode      = -1;
  for (i = 0; i < num_externs; i++) {
    tnode = (int) responsible_proc[i+NN];

    if (otnode != tnode) {
      kk     = AZ_find_index(tnode, node_name, num_neighbors);
      otnode = tnode;
    }
    node_ptr[kk]++;
  }

  tnode = otnode = 0;
  for (i = 0; i < num_neighbors; i++) {
    tnode       = node_ptr[i];
    node_ptr[i] = otnode;
    otnode     += tnode;
  }

  node_ptr[num_neighbors] = otnode;

  extern_loc_index = (int *) calloc(num_externs+1, sizeof(int));
  extern_my_index  = (int *) calloc(num_externs+1, sizeof(int));
  if ( extern_my_index == NULL) {
    (void) fprintf(stderr, "no space in AZ_exchange_rows: extern_my_index\n");
    exit(-1);
  }

  otnode = -1;
  for (i = 0; i < num_externs; i++) {
    tnode = (int) responsible_proc[i+NN];

    if (otnode != tnode) {
      kk     = AZ_find_index(tnode, node_name, num_neighbors);
      otnode = tnode;
    }
    tind = (int) local_index[i+NN];

    extern_my_index[node_ptr[kk]]  = i + NN;
    extern_loc_index[node_ptr[kk]] = tind;
    node_ptr[kk]++;
  }

  for (i = num_neighbors - 1; i > 0; i--) {
    node_ptr[i] = node_ptr[i-1];
  }
  node_ptr[0] = 0;

  for (i = 0; i < num_neighbors; i++) {
    AZ_sort(&(extern_loc_index[node_ptr[i]]), node_ptr[i+1] - node_ptr[i],
            &(extern_my_index[node_ptr[i]]), NULL);
  }

  /*
   * Throw out matrix elements that are not defined on this processor.
   * Additionally, replace matrix elements (column index) that are defined on
   * another processor, but also defined on this processor as an external point.
   */

  newstart = i = 0;
  while (i < total_num_recv) {
    while ( (i < total_num_recv) && (buffer[i] == -1.0)) {
      buffer[newstart++] = buffer[i];
      i++;
    }

    if (i == total_num_recv) break;

    tind  = (int) buffer[i];
    tnode = (int) buffer[i+2];

    if (tnode != Proc) {
      kk = AZ_find_index(tnode, node_name, num_neighbors);

      if (kk != -1) {
        k2 = AZ_find_index(tind, &(extern_loc_index[node_ptr[kk]]),
                           node_ptr[kk+1] - node_ptr[kk]);

        if (k2 != -1) {
          buffer[newstart++] = (double) extern_my_index[k2+node_ptr[kk]];
          buffer[newstart++] = buffer[i+1];
          buffer[newstart++] = (double) Proc;
        }
      }
    }
    else {
      buffer[newstart++] = buffer[i];
      buffer[newstart++] = buffer[i+1];
      buffer[newstart++] = buffer[i+2];
    }
    i += 3;
  }

  if (num_neighbors != 0) {
   if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX)
    *nonzeros = (total_num_recv - 2*new_N + 2*NN) / 3 + bindx[NN];
   else
    *nonzeros = (total_num_recv - 2*new_N + 2*NN) / 3 + indx[bpntr[blks]] + 1;
  }

  free((char *) local_index);
  free((char *) responsible_proc);
  free(node_name);
  free(node_ptr);
  free(extern_loc_index);
  free(extern_my_index);

  *buf_length = newstart;

} /* AZ_exchange_rows */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void AZ_newextract(int row, double val[], int indx[], int node, double buffer[],
                   int *next, double mapper[], double mapper2[], int NN,
                   int blks, int bindx[], int cpntr[], int bpntr[])

/*******************************************************************************

  Pull out the row denoted by 'row' and store it into 'buffer'. Rows are stored
  in the following fashion:
  a) each row is proceeded and terminated by the number -1.
  b) the first element of the row sent is the diagonal.
  c) each element is represented by a triplet of numbers (in double precision).
     The first number is the column index, the second number is the value of the
     matrix element, and the third number is the processor responsible for
     updating the value associated with the column.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  row:             Row to be stored into buffer.

  val:             Array containing the nonzero entries of the matrix (see file
                   params.txt).

  indx,
  bindx,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file params.txt).

  node:            Current processor number.

  buffer:          Buffer array where extracted rows are stored.

  next:            On ouput, next available location in 'buffer' to store a new
                   element.

  mapper:          mapper[i] corresponds to the local index of the variable 'i'
                   on the processor responsible for updating this variable.

  mapper2:         mapper2[i] corresponds to the processor responsible for
                   updating the variable 'i'.

  NN:              Order of the submatrix system.

  blks:            The total number of block rows in the matrix stored on this
                   processor.


*******************************************************************************/

{

  /* local variables */

  int    count = 0, i, start, end;
  int    diag_count;
  int    block_row;
  int    incr;
  int    icol_blk, indx_ptr, num_blk_cols, jcol, realcol, num_blk_rows;
  double the_element;

  /**************************** execution begins ******************************/

  buffer[count++] = -1.0;
  buffer[count++] = (double) row;
  diag_count      = count;
  buffer[count++] = 0.0;
  buffer[count++] = (double) node;

  /* figure out which block row we are in */

  for (block_row = 0; block_row < blks; block_row++) {
    if (cpntr[block_row+1] > row) break;
  }

  num_blk_rows = cpntr[block_row+1] - cpntr[block_row];
  incr         = row - cpntr[block_row];
  start        = bpntr[block_row];
  end          = bpntr[block_row+1];

  for (i = start; i < end; i++) {
    icol_blk     = bindx[i];
    indx_ptr     = indx[i];
    num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

    for (jcol = 0; jcol < num_blk_cols; jcol++) {
      realcol     = jcol + cpntr[icol_blk];
      the_element = val[indx_ptr+jcol*num_blk_rows+incr];

      if (the_element != 0.0) {
        if (realcol == row)
          buffer[diag_count] = the_element;
        else if (realcol >= NN) {

          /* element corresponds to an external point */

          buffer[count++] = (double) mapper[realcol];
          buffer[count++] = the_element;
          buffer[count++] = (double) mapper2[realcol];
        }
        else {
          buffer[count++] = (double) realcol;
          buffer[count++] = the_element;
          buffer[count++] = (double) node;
        }
      }
    }
  }

  buffer[count] = (double) -1;
  *next         = count + 1;

} /* newextract */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_which_block_row(int *block_row, int cpntr[], int index)

{
  *block_row = 0;
  while (cpntr[*block_row+1] <= index) (*block_row)++;

} /* which_block_row */
