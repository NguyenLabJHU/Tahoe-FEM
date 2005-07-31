/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_aztec.h,v $
 *
 * $Author: paklein $
 *
 * $Date: 2001-01-30 20:59:14 $
 *
 * $Revision: 1.1.1.1 $
 *
 * $Name: not supported by cvs2svn $
 *====================================================================*/


/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


/*
 * Include file for inclusion in any routine which will call the solver
 * library. Contains necessary constants and prototypes.
 *
 * Author:  Scott A. Hutchinson, SNL
 *          John  N. Shadid,     SNL
 *          Ray   S. Tuminaro,   SNL
 */

#ifndef __AZ_AZTEC_H_
#define __AZ_AZTEC_H_

/* PAK (04/27/98) */
#ifdef __MPI__
#include "mpi.h"
#else
#define MPI_Request int
#endif

#include <stdio.h>

/* constants */
#include "az_aztec_defs.h"

/* function definitions */
#if 0
#define max(x,y) (( x > y ) ? x : y)     /* max function  */
#define min(x,y) (( x < y ) ? x : y)     /* min function  */
#endif
#define sgn(x) ((x < 0.0) ? -1.0 : 1.0)  /* sign function */

/*
 * There are different conventions for external names for fortran subroutines.
 * In addition, different compilers return differing caluse for a fortran
 * subroutine call. In this section we take these into account.
 */

#define matched

#if defined(caps)
#   define dgemm_                     DGEMM
#   define dgemv_                     DGEMV
#   define dgetrf_                    DGETRF
#   define dgetri_                    DGETRI
#   define dgetrs_                    DGETRS
#   define dgeco_                     DGECO
#   define dgedi_                     DGEDI
#   define y12mbf_                    Y12MBF
#   define y12mcf_                    Y12MCF
#   define y12mdf_                    Y12MDF
#   define idamax_                    IDAMAX
#   define dasum_                     DASUM
#   define daxpy_                     DAXPY
#   define ddot_                      DDOT
#   define dcopy_                     DCOPY
#   define dscal_                     DSCAL
#   define dtrsm_                     DTRSM
#   define dpotrf_                    DPOTRF
#   define dlaic1_                    DLAIC1
#   define az_broadcast_              AZ_BROADCAST
#   define az_check_input_            AZ_CHECK_INPUT
#   define az_check_msr_              AZ_CHECK_MSR
#   define az_check_vbr_              AZ_CHECK_VBR
#   define az_defaults_               AZ_DEFAULTS
#   define az_exchange_bdry_          AZ_EXCHANGE_BDRY
#   define az_find_index_             AZ_FIND_INDEX
#   define az_find_local_indices_     AZ_FIND_LOCAL_INDICES
#   define az_find_procs_for_externs_ AZ_FIND_PROCS_FOR_EXTERNS
#   define az_free_memory_            AZ_FREE_MEMORY
#   define az_gavg_double_            AZ_GAVG_DOUBLE
#   define az_gdot_                   AZ_GDOT
#   define az_gmax_double_            AZ_GMAX_DOUBLE
#   define az_gmax_int_               AZ_GMAX_INT
#   define az_gmax_matrix_norm_       AZ_GMAX_MATRIX_NORM
#   define az_gmax_vec_               AZ_GMAX_VEC
#   define az_gmin_double_            AZ_GMIN_DOUBLE
#   define az_gmin_int_               AZ_GMIN_INT
#   define az_gsum_double_            AZ_GSUM_DOUBLE
#   define az_gsum_int_               AZ_GSUM_INT
#   define az_gsum_vec_               AZ_GSUM_VEC
#   define az_gvector_norm_           AZ_GVECTOR_NORM
#   define az_init_quick_find_        AZ_INIT_QUICK_FIND
#   define az_matvec_mult_            AZ_MATVEC_MULT
#   define az_msr2vbr_                AZ_MSR2VBR
#   define az_order_ele_              AZ_ORDER_ELE
#   define az_pr_error_               AZ_PR_ERROR
#   define az_print_out_              AZ_PRINT_OUT
#   define az_processor_info_         AZ_PROCESSOR_INFO
#   define az_quick_find_             AZ_QUICK_FIND
#   define az_read_msr_matrix_        AZ_READ_MSR_MATRIX
#   define az_read_update_            AZ_READ_UPDATE
#   define az_reorder_matrix_         AZ_REORDER_MATRIX
#   define az_reorder_vec_            AZ_REORDER_VEC
#   define az_set_message_info_       AZ_SET_MESSAGE_INFO
#   define az_sort_                   AZ_SORT
#   define az_solve_                  AZ_SOLVE
#   define az_transform_              AZ_TRANSFORM
#elif defined(matched)

/* do not map BLAS function names PAK (08/01/2000) */
#if 0
#   define dgemm_                     dgemm
#   define dgemv_                     dgemv
#   define idamax_                    idamax
#   define dasum_                     dasum
#   define daxpy_                     daxpy
#   define ddot_                      ddot
#   define dcopy_                     dcopy
#   define dscal_                     dscal
#   define dtrsm_                     dtrsm
#endif 

#   define dlaic1_                    dlaic1
#   define dgetrf_                    dgetrf
#   define dgetri_                    dgetri
#   define dgetrs_                    dgetrs
#   define dgeco_                     dgeco
#   define dgedi_                     dgedi
#   define dpotrf_                    dpotrf
#   define y12mbf_                    y12mbf
#   define y12mcf_                    y12mcf
#   define y12mdf_                    y12mdf
#   define az_broadcast_              az_broadcast
#   define az_check_input_            az_check_input
#   define az_check_msr_              az_check_msr
#   define az_check_vbr_              az_check_vbr
#   define az_defaults_               az_defaults
#   define az_exchange_bdry_          az_exchange_bdry
#   define az_find_index_             az_find_index
#   define az_find_local_indices_     az_find_local_indices
#   define az_find_procs_for_externs_ az_find_procs_for_externs
#   define az_free_memory_            az_free_memory
#   define az_gavg_double_            az_gavg_double
#   define az_gdot_                   az_gdot
#   define az_gmax_double_            az_gmax_double
#   define az_gmax_int_               az_gmax_int
#   define az_gmax_matrix_norm_       az_gmax_matrix_norm
#   define az_gmax_vec_               az_gmax_vec
#   define az_gmin_double_            az_gmin_double
#   define az_gmin_int_               az_gmin_int
#   define az_gsum_double_            az_gsum_double
#   define az_gsum_int_               az_gsum_int
#   define az_gsum_vec_               az_gsum_vec
#   define az_gvector_norm_           az_gvector_norm
#   define az_init_quick_find_        az_init_quick_find
#   define az_matvec_mult_            az_matvec_mult
#   define az_msr2vbr_                az_msr2vbr
#   define az_order_ele_              az_order_ele
#   define az_pr_error_               az_pr_error
#   define az_print_out_              az_print_out
#   define az_processor_info_         az_processor_info
#   define az_quick_find_             az_quick_find
#   define az_read_msr_matrix_        az_read_msr_matrix
#   define az_read_update_            az_read_update
#   define az_reorder_matrix_         az_reorder_matrix
#   define az_reorder_vec_            az_reorder_vec
#   define az_set_message_info_       az_set_message_info
#   define az_sort_                   az_sort
#   define az_solve_                  az_solve
#   define az_transform_              az_transform
#endif

#ifndef FSUB_TYPE
#  if defined(ncube)
#     define  FSUB_TYPE void
#  elif defined(paragon)
#     define  FSUB_TYPE void
#  elif defined(hp)
#     define  FSUB_TYPE void
#  else
#     define  FSUB_TYPE int
#  endif
#endif

#ifdef __cplusplus
#include <stdio.h>
extern "C" {
#endif

extern double     ddot_(int *n1, double *v1, int *dum11, double *v2,
                        int *dum21);
extern FSUB_TYPE  dscal_(int *n1, double *a1, double *r, int *stride1);
extern FSUB_TYPE  daxpy_(int *n1, double *a1, double *x, int *skipx1, double *y,
                         int *skipy1);
extern FSUB_TYPE  dcopy_(int *n1, double *src, int *src1, double *dst,
                         int *dst1);
extern FSUB_TYPE  dgetrf_(int *, int *, double *, int *, int *, int *);

extern FSUB_TYPE  dgetri_(int *, double *, int *, int *, double *, int *,
                          int *);

extern FSUB_TYPE  dgemm_(char *, char *, int *, int *, int *, double *,
                         double *, int *, double *, int *, double *, double *,
                         int *, int, int);

extern FSUB_TYPE  dgeco_(double *, int *, int *, int *, double *, double *);

extern FSUB_TYPE  dgedi_(double *, int *, int * , int *, double *, double *,
                         int *);

extern FSUB_TYPE  dgetrs_(char *, int *, int *, double *, int *, int *,
                          double *, int *, int *, int);

extern FSUB_TYPE  dlaic1_(int * , int *, double *, double *, double *, double *, double *,
                          double *, double *);

extern FSUB_TYPE  dtrsm_(char *, char *, char *, char *, int *, int *,
                         double *, double *, int *, double *, int *, int, int,
                         int, int);

extern FSUB_TYPE  dpotrf_(char *, int *, double *,int *, int *, int);

/* Aztec function prototypes that can be called by the user */

extern void AZ_async_gop(int N, double v1[], double v2[], int proc_config[],
                         double *result, int operation);

extern void AZ_solve(
        double x[],     /* On input 'x' contains the initial guess. On output*/
                        /* 'x' contains the solution to our linear system.   */
                        /* NOTE: THis vector must be of size >= N + NExt     */
        double b[],     /* right hand side of linear system.                 */
                        /* NOTE: This vector must be of size >= N            */
        int options[],
        double params[],
        int indx[],     /* The ith element of indx points to the location in */
                        /* val of the (0,0) entry of the ith block entry. The*/
                        /* last element is the number of nonzero entries of  */
                        /* matrix A plus one.                                */
        int bindx[],    /* Contains the block column indices of the non-zero */
                        /* block entries.                                    */
        int rpntr[],    /* The ith element of rpntr indicates the first point*/
                        /* row in the ith block row. The last element is the */
                        /* number of block rows plus one.                    */
        int cpntr[],    /* The jth element of cpntr indicates the first point*/
                        /* column in the jth block column. The last element  */
                        /* is the number of block columns plus one.          */
        int bpntr[],    /* The ith element of bpntr points to the first block*/
                        /* entry of the ith row in bindx. The last element is*/
                        /* the number of nonzero blocks of matrix A plus one.*/
        double val[],   /* matrix A in sparse format (VBR)  .                */
                        /* Indicates current level of factorization          */
                        /* factor_flag =                                     */
                        /*      1: indicates first call to precond. routine  */
                        /*      that performs some type of factorization     */
                        /*      preprocessing such as an incomplete LU.      */
                        /*                                                   */
                        /*      2: use preprocessing info. from a previous   */
                        /*      call. Implies some further change in the     */
                        /*      the numerical entries rather than the sparse */
                        /*      pattern.                                     */
                        /*                                                   */
                        /*      3: use precondtioner from last level 1 or 2  */
                        /*      call to precond. (see specific precondioner  */
                        /*      for more info)                               */
        int data_org[], double status[], int proc_config[]);


extern void AZ_broadcast(char *ptr, int length, int proc_config[], int action);

extern int  AZ_check_input(int data_org[], int options[], double params[],
                           int proc_config[]);

extern void AZ_check_msr(int *bindx, int N_update, int N_external,
                         int option, int *proc_config);

extern void AZ_check_update(int update[], int N_update,
                            int proc_config[]);

extern void AZ_check_vbr(int N_update, int N_external, int option,
                         int bindx[], int bnptr[], int cnptr[], int rnptr[],
                         int proc_config[]);

extern void AZ_defaults(int options[], double params[]);

extern void AZ_free_memory(int label);

extern void AZ_exchange_bdry(double x[], int data_org[]);

extern int AZ_find_simple(int, int *, int, int *, int, int *, int *);

extern int  AZ_find_index(int key, int list[], int length);

extern void AZ_find_local_indices(int N_update, int bindx[], int update[],
                                  int **external, int *N_external, int mat_type,
                                  int bpntr[]);

extern void AZ_find_procs_for_externs(int N_update, int update[],
                                      int external[], int N_external,
                                      int proc_config[], int **extern_proc);

extern double AZ_gavg_double(double var, int proc_config[]);

extern double AZ_gdot(int N, double r[], double z[], int proc_config[]);

extern double AZ_gmax_double(double, int proc_config[]);

extern int    AZ_gmax_int(int val, int proc_config[]);

extern double AZ_gmax_vec(int N, double vec[], int proc_config[]);

extern int AZ_gmin_int(int val, int proc_config[]);

extern double AZ_gmin_double(double var, int proc_config[]);

extern double AZ_gsum_double(double , int proc_config[]);

extern int    AZ_gsum_int(int totals, int proc_config[]);

extern void   AZ_gsum_vec_int(int vals[], int vals2[], int length,
                              int proc_config[]);

extern double AZ_gmax_matrix_norm(double val[], int indx[], int bindx[],
                                  int rpntr[], int cpntr[], int bpntr[],
                                  int proc_config[], int data_org[]);

extern double AZ_gvector_norm(int n, int p, double *x, int *);

extern void   AZ_init_quick_find(int list[], int length, int *shift, int *bins);

extern void AZ_list_print(int ivec[] , int length, double dvec[], int length2);

extern void   AZ_msr2vbr(double val[], int indx[], int rnptr[], int cnptr[],
                         int bnptr[], int bindx[], int msr_bindx[],
                         double msr_val[], int total_blk_rows,
                         int total_blk_cols, int blk_space, int nz_space,
                         int blk_type);

extern void   AZ_order_ele(int update_index[], int extern_index[],
                           int *internal, int *border, int N_update,
                           int msr_bindx[], int bindx[], int extern_proc[],
                           int N_external, int option, int m_type);

extern void   AZ_output_matrix(double val[], int indx[], int bindx[],
                               int rpntr[], int cpntr[], int bpntr[],
                               int proc_config[], int data_org[]);

extern void   AZ_print_error(int error_code);

extern void AZ_print_out(int update_index[], int extern_index[], int update[],
                        int external[],
                        double val[], int indx[],  int
                        bindx[], int rpntr[], int cpntr[], int bpntr[], int
                        proc_config[], int choice, int matrix, int N_update,
                        int N_external, int off_set );

extern void   AZ_processor_info(int proc_config[]);

extern int    AZ_quick_find(int key, int list[],int length, int shift,
                            int bins[]);

extern void   AZ_read_msr_matrix(int update[], double **val, int **bindx,
                                 int N_update, int proc_config[]);

extern void   AZ_read_update(int *N_update_blks, int *update_blks[],
                                 int proc_config[], int bigN, int chunk,
                                 int input_option);

/* missing function prototypes: PAK (07/24/98) */
extern int    AZ_read_external(int N_external, int external[],
                     int **extern_proc, FILE *fp, int proc_config[]);

extern void   AZ_reorder_vec(double vec[], int data_org[], int update_index[],
			     int rpntr[]);
extern void   AZ_invorder_vec(double vec[], int data_org[], int update_index[],
			     int rpntr[],double newvec[]);

extern void   AZ_reorder_matrix(int N_update, int bindx[], double val[],
                                int update_index[], int extern_index[],
                                int indx[], int rnptr[], int bnptr[],
                                int N_external, int cnptr[], int option,
                                int);

extern void   AZ_set_message_info(int N_external, int extern_index[],
                                  int N_update, int external[],
                                  int extern_proc[], int update[],
                                  int update_index[], int proc_config[],
                                  int cnptr[], int *data_org[], int);

extern void   AZ_sort(int list[], int N, int list2[], double list3[]);

extern void   AZ_matvec_mult(double *val, int *indx, int *bindx, int *rpntr,
                             int *cpntr, int *bpntr, double *b, double *c,
                             int exchange_flag, int *data_org);

extern int    md_read(char *, int , int *, int *, int *);

extern int    md_write(char *, int , int , int , int *);

extern int    md_wrap_iread(void *, int , int *, int *, MPI_Request *);
extern int    md_wrap_write(void *, int , int , int , int *);
extern int    md_wrap_wait(void *, int , int *, int *, int *, MPI_Request *);

/* added to work around using MPI_ANY_SOURCE for Irecv */
extern int md_max(int, int*);
extern int md_allgather(char*, char*, int);
/* PAK (10/07/2000) */

extern double AZ_second(void);

extern void   AZ_gdot_vec(int N, double dots[], double dots2[],
                          int proc_config[]);

extern void   AZ_transform(int proc_config[], int *external[], int bindx[],
                           double val[], int update[], int *update_index[],
                           int *extern_index[], int *data_org[], int N_update,
                           int indx[], int bnptr[], int rnptr[], int *cnptr[],
                           int mat_type);

extern void   AZ_vb2msr(int m, double val[], int indx[], int bindx[],
                        int rpntr[], int cpntr[], int bpntr[], double msr_val[],
                        int msr_bindx[]);

extern void AZ_print_vbr_matrix(
        int matrix_flag, /* = 0 no matrix output, = 1 output matrix */
        int Proc,        /* Processor number                  */
        int itotal_nodes,/* Number of internal + border nodes */
        int ext_nodes,   /* Number of external nodes          */
        double  val[],   /* matrix A in sparse format (VBR)   */
        int  indx[],     /* The ith element of indx points to the location in */
                         /* val of the (0,0) entry of the ith block entry. The*/
                         /* last element is the number of nonzero entries of  */
                         /* matrix A plus one.                                */
        int bindx[],     /* Contains the block column indices of the non-zero */
                         /* block entries.                                    */
        int rpntr[],     /* The ith element of rpntr indicates the first point*/
                         /* row in the ith block row. The last element is the */
                         /* number of block rows plus one.                    */
        int bpntr[]      /* The ith element of bpntr points to the first block*/
                         /* entry of the ith row in bindx. The last element is*/
                         /* the number of nonzero blocks of matrix A plus one.*/
        );

extern void   AZ_add_new_ele(int cnptr[], int col, int blk_row, int bindx[],
                             int bnptr[], int indx[], double val[], int therow,
                             double new_ele, int maxcols, int blk_space,
                             int nz_space, int blk_type);

extern int    AZ_find_block_col(int cnptr[], int column, int maxcols,
                                int blk_type);
extern int    AZ_find_block_in_row(int bindx[], int bnptr[], int i, int blk_col,
                                   int indx[], int, double val[], int blk_space,
                                   int nz_space);

extern void   AZ_transpose(int N, double l[], int ijl[], double lt[],
                           int ijlt[], int row_counter[]);

extern void   AZ_sortqlists(char a[], int b[], int lists[], int length,
                            int type_length, int ind_length);

extern void   AZ_mysleep(int i);

extern void   AZ_convert_values_to_ptrs(int array[], int length, int start);

extern void   AZ_convert_ptrs_to_values(int array[], int length);

extern int    AZ_find_closest_not_larger(int key, int list[], int length);

extern void   AZ_add_new_row(int therow, int *nz_ptr, int *current, double
                             **val, int **bindx, char *input,FILE *dfp,
                             int *msr_len, int *column0);

extern void   AZ_sort_ints(char a[], int indx[], int start, int end, int b[],
                           int *mid, int real_lists, char buffer[], int buf_len,
                           int afirst, int );

extern void   AZ_sort_dble(char a[], int indx[], int start, int end, int b[],
                           int *mid, int real_lists, char buffer[], int buf_len,
                           int afirst, int );

extern void   AZ_change_it(int indx[], int length, int *first, int *total,
                           int b[]);

extern void   AZ_reverse_it(int indx[], int length, int first, int total,
                            int b[]);

extern void   AZ_direct_sort(int b[], int indx[], char buffer[], char a[],
                             int *start, int buf_len, int *ind_index,
                             int *the_first, int *real_lists, int *pre_mid);

extern void   AZ_get_x_incr(double val[], int indx[], int bindx[], int rpntr[],
                            int cpntr[], int bpntr[], int options[],
                            int data_org[], int proc_config[], double params[],
                            int i, double **hh, double *rs, double *trash,
                            double **ss);

#ifdef TIME_VB
extern void   AZ_time_kernals(int , int , double , double *, int *, int *,
                              int *, int *, int *, double *, double *, int);
#endif

extern double *AZ_manage_memory(int size, int action, int type, char *name,
                                int *status);

extern void   AZ_block_diagonal_scaling(double val[], int indx[], int bindx[],
                                        int rpntr[], int cpntr[], int bpntr[],
                                        double b[], int option_i[],
                                        int data_org[], int proc_config[]);

extern void   AZ_sym_block_diagonal_scaling(double val[], int indx[],
                                            int bindx[], int rpntr[],
                                            int cpntr[], int bpntr[],
                                            double b[], int options[],
                                            int data_org[], int proc_config[]);

extern void   AZ_precondition(double val[], int indx[], int bindx[],
                              int rpntr[], int cpntr[], int bpntr[], double x[],
                              int options[], int data_org[], int proc_config[],
                              double params[]);

extern void   AZ_row_sum_scaling(double val[], int indx[], int bindx[],
                                 int rpntr[], int cpntr[], int bpntr[],
                                 double b[], int data_org[], int option_i[]);


extern void   AZ_x_scale(int proc_config[]);

extern void   AZ_sym_diagonal_scaling(double val[], int bindx[], double b[],
                                      int data_org[], int option_i[],
                                      double x[], int indx[], int bpntr[],
				      int rpntr[], int cpntr[]);

extern void   AZ_sym_rescale_sl(double x[], int data_org[]);

#ifdef next_version
extern void   AZ_sym_rescale_vbr(double x[], int data_org[]);
#endif

extern void   AZ_sym_row_sum_scaling_sl(double val[], int bindx[], double b[],
                                        int data_org[], int option_i[],
                                        double x[]);

extern void   AZ_x_scale_sl(double x[], int data_org[]);

extern void   AZ_compute_global_scalars(double val[], int indx[], int bindx[],
                                        int rpntr[], int cpntr[], int bpntr[],
                                        double x[], double b[], double r[],
                                        double w[], double *r_norm,
                                        double *scaled_r_norm, int option_i[],
                                        int data_org[], int proc_config[],
                                        int *use_r, double v1[], double v2[],
                                        double *value, int first_time);

extern void   AZ_scale_true_residual(double val[], int indx[], int bindx[],
                                     int rpntr[], int cpntr[], int bpntr[],
                                     double x[], double b[], double v[],
                                     double w[], double *actual_residual,
                                     double *scaled_r_norm, int options[],
                                     int data_org[], int proc_config[]);

extern void   AZ_compute_residual(double val[], int indx[], int bindx[],
                                  int rpntr[], int cpntr[], int bpntr[],
                                  double b[], double u[], double r[],
                                  int data_org[]);

extern void   AZ_dvbr_sparax_overlap(int m, double *val, int *indx, int *bindx,
                                     int *rpntr, int *cpntr, int *bpntr,
                                     double *b, double *c, int data_org[]);

extern void   AZ_dvbr_sparax_overlap2(int m, double *val, int *indx, int *bindx,
                                      int *rpntr, int *cpntr, int *bpntr,
                                      double *b, double *c);

extern void   AZ_dgemv2(int m, int n, double *a, double *x, double *y);

extern void   dgemv_(char *, int *, int *, double *, double *, int *, double *,
                     int *, double *, double *, int *, int);

/* When calling this fortran routine from C we need to include an extra     */
/* parameter on the end indicating the string length of the first parameter */

#ifdef hp
extern void   dgemvnsqr_(int *, double *, double *, double *);
#endif

extern void   AZ_domain_decomp(double x[], double val[], int indx[],
                               int rpntr[], int cpntr[], int bindx[],
                               int bpntr[], int options[], int data_org[],
                               int proc_config[], double params[]);

extern void   AZ_polynomial_expansion(double val[], int indx[], int bindx[],
                                      int rpntr[], int cpntr[], int bpntr[],
                                      double z[], int options[], int data_org[],
                                      int proc_config[]);

extern void   AZ_dvbr_diag_sparax(int m, double *val, int *rpntr, int *bpntr,
				  double *b, double *c);

#if defined (hp)
extern void   vec_$dcopy(double *, double *, int *);
extern void   blas_$dgemm(char *, char *, int *, int *, int *, double *,
                          double *, int *, double *, int *, double *, double *,
                          int *, int, int);
#endif

extern void   AZ_calc_blk_diag_inv(double *val, int *indx, int *bindx,
                                   int *rpntr, int *cpntr, int *bpntr,
                                   double *d_inv, int *d_indx, int *d_bindx,
                                   int *d_rpntr, int *d_bpntr, int data_org[]);

extern void   AZ_get_poly_coefficients(int power, double b, double c[],
                                       int param_flag);

extern void   AZ_reorganize_matrix(int bindx[], double val[], int N,
                                   int extra[]);

extern void   AZ_lower(double rhs[], double val[],int bindx[], int extra[],
                       int n);

extern void   AZ_upper(double rhs[],double val[], int bindx[], int extra[],
                       int n);

extern void   AZ_msr2ilu_setup(double val[], int bindx[], double val2[],
                               int bindx2[], int new_N, int NN, int options[],
                               int read_index, double buffer[], int data_org[]);

extern double AZ_sum_comp(int first1, int last1, int first2, int last2,
                          double newval[], int newbindx[]);

extern int    idamax_(int *, double *, int *);

extern double dasum_(int *, double *, int *);

extern void   AZ_factandsolve(double newa[], double aflag[], double pivot[],
                              double x[], int snr[], int rnr[], int ha[],
                              int iflag[], int *z, int *ifail, int *nn, int *n,
                              int *iha, int *nn1);

extern void   AZ_backsolve(double newa[], double pivot[], double x[], int snr[],
                           int ha[], int iflag[], int *ifail, int *nn, int *n,
                           int *iha);

extern void   AZ_msr2lu(int oldN, double val[], double newa[], int snr[],
                        int rnr[], int *z, int indx[]);

extern void   AZ_vb2lu(int m, double val[], int indx[], int bindx[],
                       int rpntr[], int cpntr[], int bpntr[], double newa[],
                       int snr[], int rnr[], int *nz_ptr);

extern void   y12mbf_(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, int ha[], int *iha, double aflag[],
                      int iflag[], int *ifail);

extern void   y12mcf_(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, double pivot[], double b[], int ha[],
                      int *iha, double aflag[], int iflag[], int *ifail);

extern void   y12mdf_(int *n, double val[], int *nn, double b[], double pivot[],
                      int snr[], int ha[], int *iha, int iflag[], int *ifail);

extern void   AZ_print_call_iter_solve(int * , int );

extern int    AZ_check_options(int * , int ,int data_org[], int,double *);

extern int    AZ_get_block(int j, int k, int bindx[], int bpntr[], int *ptr_j);

extern void   AZ_update_block(int i, int k, int j, double val[], int indx[],
                              int bindx[], int cpntr[]);

extern void   AZ_divide_block(int i, int j, double val[], int indx[],
                              int bindx[], int cpntr[], double *z,
                              double *blockj, double *blocki, int *ipvt);

/* prototype typo? - PAK (07/24/98) */
extern void   AZ_divide_block0(int i, int j, double val[], int indx[],
                                  int bindx[], int cpntr[], int *ipvt);

extern void   AZ_AZ_divide_block0(int i, int j, double val[], int indx[],
                                  int bindx[], int cpntr[], int *ipvt);

extern void   AZ_which_block_row(int *block_row, int cpntr[], int row_index);

extern void   AZ_vb2bilu_setup(int cpntr2[], int bpntr2[],
                               int bindx2[], int indx2[], double val2[],
                               int rpntr[], int cpntr[], int bpntr[],
                               int bindx[], int indx[], double val[], int newN,
                               int old_blks, int new_blks, int length,
                               double buffer[], int options[], int data_org[]);

extern void   AZ_order(int M, double *val_old, double *val_new, int *bindx,
                       int *indx_old, int *indx_new, int *bpntr,
                       int *diag_bloc);

extern void   AZ_lower_triang_vbr_solve(int N, int M, double *val, int *indx,
                                        int *bindx, int *rpntr, int *cpntr,
                                        int *bpntr, int *diag_block, double *y,
                                        double *b);

extern void   AZ_upper_triang_vbr_solve(int N, int M, double *val, int *indx,
                                        int *bindx, int *rpntr, int *cpntr,
                                        int *bpntr, int *diag_block, double *y,
                                        double *b,int *);

extern void   AZ_dgemv3(int m, int n, double *a, double *x, double *y);

extern void   AZ_gather_ptrs(int i, int N, double val[], int bindx[],
                             int *nptrs, int j_ptrs[], double a_j_ptrs[]);

extern void   AZ_sort_ptrs(int n, int arr[], double brr[]);

double        AZ_get_elm(int i, int j, double val[], int bindx[]);

extern void   AZ_extract(int row, double val[], int bindx[], int node,
                         double buffer[], int *next, double mapper[],
                         double mapper2[], int NN);

extern void   AZ_newextract(int row, double val[], int msr_bindx[], int node,
                            double buffer[], int *next, double mapper[],
                            double mapper2[], int NN, int blks, int bindx[],
                            int cpntr[], int bpntr[]);

extern void   AZ_ilu_routine(int Nexp, int Mexp, double x[], int length,
                             double buffer[], double val[], int indx[],
                             int bindx[], int rpntr[], int cpntr[], int bpntr[],
                             int option_i[], int data_org[]);

extern void   AZ_block_ilu(double val2[], int indx2[], int bindx2[],
                           int cpntr2[], int bpntr2[],
                           int new_blks, int new_N, int length, double buffer[],
                           double x[], int option_i[], int data_org[]);

extern void   AZ_exchange_rows(double val[], int indx[], int *Mexp, int Nexp,
                               double **buffer, int *length, int bindx[],
                               int cpntr[], int bpntr[], int data_org[],
			       int proc_config[]);

extern void   AZ_sym_gauss_seidel(void);

extern void   AZ_sym_gauss_seidel_sl(double val[], int bindx[], double x[],
                                     int data_org[], int options[]);

extern void   AZ_scale_f(double val[], int indx[], int bindx[], int rpntr[],
                         int cpntr[], int bpntr[], double rhs[], double x[],
                         int options[], int data_org[], int proc_config[],
                         int action);

extern int    AZ_breakdown_f(int N, double v[], double w[], double inner,
                             int proc_config[]);

extern void   AZ_change_sign(double *lambda_max, double val[], int indx[],
                             int bindx[], int rpntr[], int cpntr[], int bpntr[],
                             int data_org[]);

/* potentially available to the user */

extern void   AZ_exchange_local_info(int N_neighbors, int proc_num_neighbor[],
                                     char *message_send_add[],
                                     int message_send_length[],
                                     char *message_recv_add[],
                                     int message_recv_length[], int type);

extern double AZ_sync_timer(int proc_config[]);

extern void   AZ_sync(int, int);

extern void   AZ_gappend(int vals[], int *cur_length, int total_length,
                         int proc_config[]);

extern void   AZ_pgmres(double *, int *, int *, int *, int *, int *, double *,
                        double *, double *, int *, double * , int *, int *,
                        double *);

extern void   AZ_pcg_f(double *, int *, int *, int *, int *, int *, double *,
                       double *, double *, int *, double * , int *, int * ,
                       double *);

extern void   AZ_pcgs(double *, int *, int *, int *, int *, int *, double *,
                      double *, double *, int *, double * , int *, int * ,
                      double *);

extern void   AZ_psymmlq(double *, int *, int *, int *, int *, int *, double *,
                        double *, double *, int *, double * , int *, int * ,
                        double *);

extern void   AZ_pqmrs(double *, int *, int *, int *, int *, int *, double *,
                       double *, double *, int *, double *, int *, int *,
                       double *);

extern void   AZ_pbicgstab(double *, int *, int *, int *, int *, int *,
                           double *, double *, double *, int *, double * ,
                           int *, int * , double *);

extern void   AZ_lu_y12m(int, int, double [], int, double [], double [], int [],
                         int [], int [],int [],int [], int [], int [], double);

extern int    AZ_broadcast_info(char buffer[], int proc_config[], int length);

extern void   AZ_gather_mesg_info(double x[],int data_org[],char **, char **,
                                  int *, int *);

extern void   AZ_read_local_info(int data_org[], char *message_recv_add[],
                                 int message_recv_length[]);

extern void   AZ_write_local_info(int data_org[], char *message_recv_add[],
                                  char *message_send_add[],
                                  int message_recv_length[],
                                  int message_send_length[]);

extern double AZ_calc_solve_flops(int options[], int, double , int , double,
                                  int data_org[], int proc_config[]);

extern double AZ_calc_iter_flops(int solver_flag, double inner_flops,
                                 double daxpy_flops, double matvec_flops,
                                 int total_its, double gnnz, double K);

extern double AZ_calc_precond_flops(int solver_flag, int options[],
                                    double daxpy_flops,
                                    double matvec_flops, int total_its, int gn,
                                    double gnnz, int data_org[],
                                    int proc_config[]);

extern void   AZ_print_sync_start(int proc, int do_print_line);
extern void   AZ_print_sync_end(int proc, int nprocs, int do_print_line);

extern void   AZ_p_error(char *str, int proc);

extern int    AZ_get_new_eps(double *epsilon, double, double,
                             int proc_config[]);

extern void   AZ_terminate_status_print(int situation, int iter,
                                        double status[], double rec_residual,
                                        double params[], double scaled_r_norm,
                                        double actual_residual, int options[],
                                        int proc_config[]);

extern void   AZ_random_vector(double u[], int data_org[], int proc_config[]);

extern void   AZ_splitup_big_msg(int num_neighbors, double *buffer,
                                 int *start_send_proc, int *actual_send_length,
                                 int *num_nonzeros_recv, int *proc_num_neighbor,
                                 int type, int *total_num_recv,
                                 int *proc_config);

extern void   AZ_dtrans(int *, int *, double *);

extern int    AZ_get_sym_indx(int, int, int *, int *, int *);

extern double AZ_srandom1(int *seed);

#ifdef eigen
extern void AZ_lu_ng(int n, int nonzeros, double x[], double val[], int 
		     bindx[], int options[], int data_org[], double 
                     drop_tolerance);
#endif

#ifdef __cplusplus
}
#endif

#endif /* _AZ_AZTEC_H_ */
