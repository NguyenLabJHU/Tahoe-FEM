/* $Id: debug.c,v 1.1 2005-01-05 07:37:08 paklein Exp $ */
#include <stdio.h>
#include "pspases_f2c.h"

void dsyrk_(char *UL, char *NT, integer *N, integer *K,
	doublereal *alpha, doublereal *A, integer *lda, 
	doublereal *beta, doublereal *C, integer *ldc, 
	ftnlen UL_len, ftnlen NT_len);

#if 0
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
#endif
   
void mydsyrk_(char *UL, char *NT, integer *N, integer *K,
	doublereal *alpha, doublereal *A, integer *lda, 
	doublereal *beta, doublereal *C, integer *ldc, 
	ftnlen UL_len, ftnlen NT_len)
{
	int i;

	/* report values */
	printf("dsyrk:\n");
	printf("UT = \"");
	for (i = 0; i < UL_len; i++)
		printf("%c", UL[i]);
	printf("\"\n");
	printf("  N = %d\n", *N);
	printf("  K = %d\n", *K);
	printf("lda = %d\n", *lda);
	printf("ldc = %d\n", *ldc);
	printf("NT = \"");
	for (i = 0; i < NT_len; i++)
		printf("%c", NT[i]);
	printf("\"\n");

	/* call */
	dsyrk_(UL, NT, N, K, alpha, A, lda, beta, C, ldc, UL_len, NT_len);
}
