/*****************************************************************************/
/*                                                                           */
/*   (C) Copyright IBM Corporation, 1997                                     */
/*   (C) Copyright Regents of the University of Minnesota, 1997              */
/*                                                                           */
/*   parmetisf.c                                                             */
/*                                                                           */
/*   Written by Mahesh Joshi, U of MN.                                       */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/* This code is meant to be used solely for educational, research, and       */
/* benchmarking purposes by non-profit institutions and US government        */
/* agencies only.  Use by any other organization requires prior written      */
/* permission from both IBM Corporation and the University of Minnesota.     */
/* The software may not be sold or redistributed.  One may make copies       */
/* of the software or modify it for their use provided that the copies,      */
/* modified or otherwise, are not sold or distributed, are used under the    */
/* same terms and conditions, and this notice and any part of the source     */
/* code that follows this notice are not separated.                          */
/*                                                                           */
/* As unestablished research software, this code is provided on an           */
/* ``as is'' basis without warranty of any kind, either expressed or         */
/* implied, including but not limited to implied warranties of               */
/* merchantability and fitness for a particular purpose.  IBM does not       */
/* warrant that the functions contained in this software will meet the       */
/* user's requirements or that the operation of its routines will be         */
/* uninterrupted or error-free.  Acceptance and use of this program          */
/* constitutes the user's understanding that he/she will have no recourse    */
/* to IBM for any actual or consequential damages, including, but not        */
/* limited to, lost profits or savings, arising out of the use or inability  */
/* to use these libraries.  Even if the user informs IBM of the possibility  */
/* of such damages, IBM expects the user to accept the risk of any such      */
/* harm, or the user shall not attempt to use these libraries for any        */
/* purpose.                                                                  */
/*                                                                           */
/* The downloading, compiling, or executing any part of this software        */
/* constitutes an implicit agreement to these terms.  These terms and        */
/* conditions are subject to change at any time without prior notice.        */
/*                                                                           */
/*****************************************************************************/
/* $Id: parmetisf.c,v 1.2 2004-12-11 09:27:22 paklein Exp $ */
/*****************************************************************************/

#include <pspaces.h>
#include "mpi.h"

/* #define __DO_DEBUG__ 1 */
#undef __DO_DEBUG__

void parometisf(idxtype *vtxdist,idxtype *xadj,idxtype *adjncy,
		idxtype *order,idxtype *sizes,int *options,
		int *serialorder,MPI_Comm *commin) 
{
	int i, rank;
	int numflag;

#if __DO_DEBUG__
	MPI_Comm_rank(*commin, &rank);
	printf("%5d: vtxdist = ", rank);
	for (i = 0; i < 10; i++)
		printf("%8d", vtxdist[i]);
	printf("\n");

	printf("%5d: xadj = ", rank);
	for (i = 0; i < 10; i++)
		printf("%8d", xadj[i]);
	printf("\n");

	printf("%5d: sizes = ", rank);
	for (i = 0; i < 10; i++)
		printf("%12d", sizes[i]);
	printf("\n");

	printf("%5d: options = ", rank);
	for (i = 0; i < 5; i++)
		printf("%8d", options[i]);
	printf("\n");
#endif

  numflag = 0;
  if(*serialorder) {
    ParMETIS_SerialNodeND(
                vtxdist,xadj,adjncy,&numflag,options,order,sizes,commin);  
  }
  else {
    ParMETIS_NodeND(
                vtxdist,xadj,adjncy,&numflag,options,order,sizes,commin);  
  }

}

void ikeysortf(int *pn,idxtype *inds,idxtype *perm) {

  int n,i;
  KeyValueType *base;

  n = *pn;
  base = (KeyValueType *)malloc(n*sizeof(KeyValueType));

  for(i=0;i<n;i++) {
    base[i].key = inds[i];
    base[i].val = i;
  }

  ikeysort(n,base);

  for(i=0;i<n;i++) {
    inds[i] = base[i].key;
    perm[i] = base[i].val;
  }

  free(base);

}
