/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: az_aztec_defs.h,v $
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


/*******************************************************************************
 *
 *              message types used in communication
 *
 ******************************************************************************/


/* message space used corresponds to AZ_MSG_TYPE --> AZ_MSG_TYPE + NUM_MSGS */

/* In general, AZTEC uses message types that lie between the values: AZ_MSG_TYPE
   and AZ_MSG_TYPE + AZ_NUM_MSGS. Each time we reach a code segment which every
   processor executes, the current message type to be used (the global variable
   AZ_sys_msg_type) is incremeted (modulo NUM_MSGS) to determine the new message
   type to use. */

#ifndef _AZ_AZTEC_DEF_H_
#define _AZ_AZTEC_DEF_H_

#define AZ_MSG_TYPE      1234
#define AZ_NUM_MSGS        20

/*******************************************************************************
 *
 *              various constants
 *
 ******************************************************************************/

#ifndef AZ_MAX_MEMORY_SIZE
#define AZ_MAX_MEMORY_SIZE   16000000  /* maximum memory size used for the LU */
                                       /* within domain decomposition.        */
#endif
#ifndef AZ_MAX_MSG_BUFF_SIZE
#define AZ_MAX_MSG_BUFF_SIZE 100000    /* max usable message buffer size      */
#endif
#define AZ_MAX_NEIGHBORS     250
#define AZ_MAX_MESSAGE_SIZE  (AZ_MAX_MSG_BUFF_SIZE / (2*AZ_MAX_NEIGHBORS))
#define AZ_FALSE               0
#define AZ_TRUE                1
#define AZ_MAX_POLY_ORDER     10 /* max order for polynomial preconditioners */
#define AZ_default           -10 /* options[i] = AZ_default ==>
                                        AZ_check_input() sets options[i] to
                                        its default value.
                                  */

/*******************************************************************************
 *
 *              Array sizes for declarations
 *
 ******************************************************************************/

#define AZ_OPTIONS_SIZE       15
#define AZ_PARAMS_SIZE         5
#define AZ_PROC_SIZE           3
#define AZ_STATUS_SIZE         7
#define AZ_COMM_SIZE          (10 + 3*AZ_MAX_NEIGHBORS)

/*******************************************************************************
 *
 *              constants for solver types
 *
 ******************************************************************************/

#define AZ_cg               0 /* preconditioned conjugate gradient method     */
#define AZ_gmres            1 /* preconditioned gmres method                  */
#define AZ_cgs              2 /* preconditioned cg squared method             */
#define AZ_tfqmr            3 /* preconditioned transpose-free qmr method     */
#define AZ_bicgstab         4 /* preconditioned stabilized bi-cg method       */
#define AZ_slu              5 /* super LU direct method.                      */
#define AZ_symmlq           6 /* indefinite symmetric like symmlq             */
#define AZ_lu               8 /* sparse LU direct method. Also used for a     */
                              /* preconditioning option.  NOTE: this should   */
                              /* be the last solver so that AZ_check_input()  */
                              /* works properly.                              */

/*******************************************************************************
 *
 *              constants for scaling types
 *
 ******************************************************************************/

/* #define AZ_none          0    no scaling                                   */
/* #define AZ_Jacobi        1    Jacobi scaling                               */
#define AZ_BJacobi          2 /* block Jacobi scaling                         */
#define AZ_row_sum          3 /* point row-sum scaling                        */
#define AZ_sym_diag         4 /* symmetric diagonal scaling                   */
#define AZ_sym_row_sum      5 /* symmetric diagonal scaling                   */
#define AZ_sym_BJacobi      6 /* symmetric block Jacobi scaling. NOTE: this   */
                              /* should be last so that AZ_check_input()      */
                              /* works properly.                              */

/*******************************************************************************
 *
 *              constants for preconditioner types
 *
 ******************************************************************************/

#define AZ_none             0 /* no preconditioning. Note: also used for      */
                              /* scaling, output, overlap options options     */
#define AZ_Jacobi           1 /* Jacobi preconditioning. Note: also used for  */
                              /* scaling options                              */
#define AZ_sym_GS           2 /* symmetric Gauss-Siedel preconditioning       */
#define AZ_Neumann          3 /* Neumann series polynomial preconditioning    */
#define AZ_ls               4 /* least-squares polynomial preconditioning     */
#define AZ_ilu              6 /* domain decomp with  ilu in subdomains        */
#define AZ_bilu             7 /* domain decomp with block ilu in subdomains   */
/* #define AZ_lu            8    domain decomp with   lu in subdomains        */
#define AZ_icc              9 /* domain decomp with incomp Choleski in domains*/
                              /* NOTE: AZ_icc should be last so that          */
                              /* AZ_check_input() works properly.             */

/*******************************************************************************
 *
 *              constants for convergence types
 *
 ******************************************************************************/

#define AZ_r0               0 /* ||r||_2 / ||r^{(0)}||_2                      */
#define AZ_rhs              1 /* ||r||_2 / ||b||_2                            */
#define AZ_Anorm            2 /* ||r||_2 / ||A||_infty                        */
#define AZ_sol              3 /* ||r||_infty/(||A||_infty ||x||_1+||b||_infty)*/
#define AZ_weighted         4 /* ||r||_WRMS                                   */
                              /* NOTE: AZ_weighted should be last so that     */
                              /* AZ_check_input() works properly.             */

/*******************************************************************************
 *
 *              constants for output types
 *
 ******************************************************************************/

#define AZ_all             -3 /* Print out everything including matrix        */
                              /* Must be lowest value so that AZ_check_input()*/
                              /* works properly.                              */
/* #define AZ_none          0    Print out no results (not even warnings)     */
#define AZ_last            -1 /* Print out final residual and warnings        */
#define AZ_warnings        -2 /* Print out only warning messages              */

/*******************************************************************************
 *
 *              constants for matrix output
 *
 ******************************************************************************/

#define AZ_input_form       0 /* Print out the matrix arrays as they appear   */
                              /* along with some additional information. The  */
                              /* idea here is to print out the information    */
                              /* that the user must supply as input to the    */
                              /* function AZ_transform()                      */
#define AZ_global_mat       1 /* Print out the matrix as a(i,j) where i and j */
                              /* are the global indices. This option must     */
                              /* be invoked only after AZ_transform() as the  */
                              /* array update_index[] is used.                */
                              /* NOTE: for VBR matrices the matrix is printed */
                              /* as a(I(i),J(j)) where I is the global block  */
                              /* row and J is the global block column and i   */
                              /* and j are the row and column indices within  */
                              /* the block.                                   */
#define AZ_explicit         2 /* Print out the matrix as a(i,j) where i and j */
                              /* are the local indices.                       */
                              /* NOTE: for VBR matrices the matrix is printed */
                              /* as a(I(i),J(j)) where I is the global block  */
                              /* row and J is the global block column and i   */
                              /* and j are the row and column indices within  */
                              /* the block.                                   */

/*******************************************************************************
 *
 *              constants for using factorization information
 *
 ******************************************************************************/

#define AZ_calc             1 /* use no previous information                  */
#define AZ_recalc           2 /* use last symbolic information                */
#define AZ_reuse            3 /* use a previous factorization to precondition */
#define AZ_sys_reuse        4 /* use last factorization to precondition       */
                              /* NOTE: AZ_sys_reuse should be last so that    */
                              /* AZ_check_input() works properly.             */

/*******************************************************************************
 *
 *              constants for domain decompositon overlap
 *
 ******************************************************************************/

/* #define AZ_none          0    No overlap                                   */
#define AZ_diag             1 /* Use diagonal blocks for overlapping          */
#define AZ_sym_full         2 /* Use external rows (symmetric) for overlapping*/
#define AZ_full             3 /* Use external rows   for overlapping          */
                              /* Note: must be highest value so that          */
                              /*       AZ_check_input() works properly.       */

/*******************************************************************************
 *
 *              constants for GMRES orthogonalization procedure
 *
 ******************************************************************************/

#define AZ_classic          0
#define AZ_modified         1

/*******************************************************************************
 *
 *              constants for determining rtilda (used in bicgstab, cgs, tfqmr)
 *
 ******************************************************************************/

#define AZ_resid            0
#define AZ_rand             1

/*******************************************************************************
 *
 *              constants indicating reason for iterative method termination
 *
 ******************************************************************************/

#define AZ_normal           0 /* normal termination                           */
#define AZ_param            1 /* requested option not implemented             */
#define AZ_breakdown        2 /* numerical breakdown during the computation   */
#define AZ_maxits           3 /* maximum iterations exceeded                  */
#define AZ_loss             4 /* loss of precision                            */
#define AZ_ill_cond         5 /* GMRES hessenberg is ill-conditioned          */

/*******************************************************************************
 *
 *              array indices into options array
 *
 ******************************************************************************/

#define AZ_solver              0
#define AZ_scaling             1
#define AZ_precond             2
#define AZ_conv                3
#define AZ_output              4
#define AZ_pre_calc            5
#define AZ_max_iter            6
#define AZ_poly_ord            7
#define AZ_overlap             8
#define AZ_kspace              9
#define AZ_orthog              10
#define AZ_aux_vec             11
#define AZ_print_freq          12

/*******************************************************************************
 *
 *              array indices into params array
 *
 ******************************************************************************/

#define AZ_tol                 0
#define AZ_drop                2
#define AZ_weights             3

/*******************************************************************************
 *
 *              array indices into data_org array
 *
 ******************************************************************************/

#define AZ_matrix_type         0
#define AZ_N_internal          1
#define AZ_N_border            2
#define AZ_N_external          3
#define AZ_N_int_blk           4
#define AZ_N_bord_blk          5
#define AZ_N_ext_blk           6
#define AZ_N_neigh             7
#define AZ_total_send          8
#define AZ_name                9
#define AZ_neighbors           10
#define AZ_rec_length          (10 +   AZ_MAX_NEIGHBORS)
#define AZ_send_length         (10 + 2*AZ_MAX_NEIGHBORS)
#define AZ_send_list           (10 + 3*AZ_MAX_NEIGHBORS)

/*******************************************************************************
 *
 *              array indices into status array
 *
 ******************************************************************************/

#define AZ_its                 0
#define AZ_why                 1
#define AZ_r                   2
#define AZ_rec_r               3
#define AZ_scaled_r            4
#define AZ_first_precond       5     /* This is used to record the time for */
                                     /* the first preconditioning step. The */
                                     /* intention is time factorization     */
                                     /* routines. Note: not mentioned in    */
                                     /* manual                              */

/*******************************************************************************
 *
 *              array indices into proc_config array
 *
 ******************************************************************************/

#define AZ_node                0
#define AZ_N_procs             1
#define AZ_dim                 2

/*******************************************************************************
 *
 *              partitioning option choices
 *
 ******************************************************************************/

#define AZ_linear              0
#define AZ_file                1
#define AZ_box                 2

/*******************************************************************************
 *
 *              constants for memory management
 *
 ******************************************************************************/

#define AZ_ALLOC               0
#define AZ_CLEAR               1
#define AZ_REALLOC             2
#define AZ_SYS                 -14901
#define AZ_OLD_ADDRESS         0
#define AZ_NEW_ADDRESS         1

/*******************************************************************************
 *
 *              constants for matrix types
 *
 ******************************************************************************/

#define AZ_MSR_MATRIX          0
#define AZ_VBR_MATRIX          1
#define AZ_USER_MATRIX         2

/*******************************************************************************
 *
 *              constants for scaling action
 *
 ******************************************************************************/

#define AZ_SCALE_MAT           0
#define AZ_RESCALE_SOL         1

/*******************************************************************************
 *
 *              constants used for residual expresion calculations
 *              (performed by AZ_compute_global_scalars) within iterative methods
 *
 ******************************************************************************/

#define AZ_NOT_FIRST           0   /* not the first residual expression       */
                                   /* request. Information should be available*/
                                   /* from a previous request and the residual*/
                                   /* may not be available.                   */
#define AZ_FIRST_TIME          1   /* first time that a residual expression   */
                                   /* is requested for a particular iterative */
                                   /* solve. This means that the true residual*/
                                   /* is available and that certain invariant */
                                   /* information (e.g. r_0, ||A||) must be   */
                                   /* computed.                               */

/*******************************************************************************
 *
 *              constants (see AZ_get_new_eps) used to determine whether to
 *              continue or to quit the iterative method when the real
 *              residual does not match the updated residual
 *
 ******************************************************************************/

#define AZ_QUIT             5
#define AZ_CONTINUE         6

/*******************************************************************************
 *
 *              constants (see AZ_broadcast) used to determine whether to
 *              concatenate information to be broadcast or to send information
 *              already stored in an internal buffer
 *
 ******************************************************************************/

#define AZ_PACK             0
#define AZ_SEND             1


/*******************************************************************************
 *
 *              software tool constants
 *
 ******************************************************************************/

#define AZ_TEST_ELE         3
#define AZ_ALL              1 /* All elements are reordered.                  */
#define AZ_EXTERNS          2 /* Only external elements are reordered.        */
#define AZ_GLOBAL           1 /* MSR entries correspond to global columns     */
#define AZ_LOCAL            2 /* MSR entries correspond to local columns      */

#endif /* _AZ_AZTEC_DEF_H_ */
